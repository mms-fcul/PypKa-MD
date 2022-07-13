import io
import logging
import os
import subprocess as sb
import time
from contextlib import redirect_stdout
from datetime import datetime, timedelta
from typing import Dict, Tuple, List

import sys

sys.path.insert(0, "/home/pedror/MMS@FCUL/pypka/")

from pypka import Titration

from pdbmender.formats import read_gro_line

from pypkamd.configs import Config
from pypkamd.misc import create_link, get_curtime, remove_comments


def create_cphdm_directory() -> None:
    os.system("rm -rf {}".format(Config.md_configs.tmpDIR))
    os.system("mkdir -p {}".format(Config.md_configs.tmpDIR))

    os.system(
        "cp -f {} {}/{}".format(
            Config.md_configs.TOPin, Config.md_configs.tmpDIR, Config.md_configs.TOP
        )
    )
    os.system(
        "cp -f {} {}/{}".format(
            Config.md_configs.GROin, Config.md_configs.tmpDIR, Config.md_configs.GRO
        )
    )
    os.system("cp -f {} {}/.".format(Config.md_configs.MDPin, Config.md_configs.tmpDIR))

    os.chdir("{}".format(Config.md_configs.tmpDIR))

    create_link("{}".format(Config.md_configs.ffDIR))
    # create_link("{}/residuetypes.dat".format(Config.md_configs.ffDIR))

    relax_mdp_new_params = [
        "nsteps",
        "Pcoupl",
        "freezegrps",
        "freezedim",
    ]

    with open(Config.md_configs.MDPin) as f, open(
        Config.md_configs.MDP, "w"
    ) as fmdp_new, open(Config.md_configs.MDP_relax, "w") as frelax_new:
        new_mdp_content = ""
        relax_mdp_content = ""
        cutoff_scheme = ""
        for line in f:
            cleaned_line = remove_comments(line)
            parts = cleaned_line.split("=")
            if len(parts) >= 2:
                param = parts[0].strip()
                if param not in relax_mdp_new_params:
                    relax_mdp_content += line
                if param == "cutoff-scheme":
                    cutoff_scheme = parts[1].split(";")[0].strip().lower()

            new_mdp_content += line

        preffix_energygrp_excl = ""
        if cutoff_scheme == "verlet":
            preffix_energygrp_excl = ";"

        relax_mdp_content += """
;Solvent relaxation parameters

nsteps = {0}
Pcoupl = No
;constraints = none
;energygrps = {1}
freezegrps = {1}
freezedim = Y Y Y
{2}energygrp_excl = {1} {1}
        """.format(
            Config.md_configs.RelaxSteps,
            Config.md_configs.titrating_group,
            preffix_energygrp_excl,
        )

        fmdp_new.write(new_mdp_content)
        frelax_new.write(relax_mdp_content)
    os.remove(Config.md_configs.MDPin)

    # Create a Index File
    if Config.md_configs.NDXin:
        os.system("cp ../{} {}".format(Config.md_configs.NDXin, Config.md_configs.NDX))
        print("cp ../{} {}".format(Config.md_configs.NDXin, Config.md_configs.NDX))
    else:
        new_ndx_cmd = "{}/gmx make_ndx -f {} -o {} " "-quiet".format(
            Config.md_configs.GroDIR, Config.md_configs.GRO, Config.md_configs.NDX
        )
        sb.run(
            new_ndx_cmd,
            shell=True,
            input=b"q\n",
            stdout=Config.md_configs.LOG,
            stderr=sb.STDOUT,
        )

    # Create a .tpr File to be used in the centering procedure
    new_tpr_cmd = (
        "{}/gmx grompp -f {} -c {} -p {} -n {} -o fixgro.tpr "
        "-maxwarn 1000 -quiet".format(
            Config.md_configs.GroDIR,
            Config.md_configs.MDP,
            Config.md_configs.GRO,
            Config.md_configs.TOP,
            Config.md_configs.NDX,
        )
    )
    sb.run(new_tpr_cmd, shell=True, stdout=Config.md_configs.LOG, stderr=sb.STDOUT)

    # Create a .dat file for fixbox
    if Config.md_configs.DATin:
        os.system("cp ../{} {}".format(Config.md_configs.DATin, Config.md_configs.DAT))
    else:
        with open(Config.md_configs.NDX) as f:
            save_trigger = False
            first_index = None
            last_index = None
            for line in f:
                if line.startswith("[ Protein ]"):
                    save_trigger = True
                elif line.startswith("["):
                    save_trigger = False
                elif save_trigger:
                    if not first_index:
                        first_index = line.split()[0]
                    last_index = line.split()[-1]

        with open(Config.md_configs.DAT, "w") as f:
            content = """# Molecular definitions file.
# All lines not starting with the following characters are ignored:
#   G : group name definition, followed by its molecule definitions.
#       a : molecule defined by atom (ordinal) index range.
#       n : molecule(s) defined by residue name.
#       g : molecules(s) defined through previously defined group.
#   A : group for assemblage stage, one line per (sequential) stage.
#   C : groups to be centered along the three box vectors (a,b,c or v1,v2,v3).

G Protein
#1st chain goes from 1st atom to atom 5743 (ordinal, not atom numbers)
a    {} {}

# Groups to be (sequentially) assembled:
A Protein

# Groups to be centered along each of the three box vectors:
C Protein Protein Protein W W W

# Use PBC (1) or not (0) along each of the three box vectors:
P Protein Protein Protein
    """.format(
                first_index, last_index
            )
            f.write(content)


def center_titrable_molecule(titrating_group: str, gro: str, ndx: str) -> str:
    # Centering procedure
    # effective.gro -> fixgro_input.gro -> centered.gro

    f_out_name = "centered.gro"

    fixgro_cmd = (
        "{}/gmx trjconv -f {} "
        "-o fixgro_input.gro -s fixgro.tpr -n {} "
        "-pbc mol -quiet; ".format(Config.md_configs.GroDIR, gro, ndx)
    )
    fixgro_cmd += "{}/fixbox fixgro_input.gro {} " "> {}".format(
        Config.md_configs.fixboxDIR, Config.md_configs.DAT, f_out_name
    )

    input_str = "{}\n".format(titrating_group)
    sb.run(
        fixgro_cmd,
        shell=True,
        input=str.encode(input_str),
        stdout=Config.md_configs.LOG,
        stderr=sb.STDOUT,
    )

    if os.path.isfile("fixgro_input.gro"):
        rm_cmd = "rm fixgro_input.gro"
        sb.run(rm_cmd, shell=True)
    else:
        raise Exception(
            f"Centering procedure failed. Check {Config.md_configs.LOG_fname} for more info."
        )

    return f_out_name


def run_dynamics(mdp: str, gro: str, top: str, ndx: str, sysname: str) -> None:
    tpr_cmd = (
        "{0}/gmx grompp -f {1} -c {2} -p {3} -n {4} -o {5}.tpr -po {5}_out.mdp "
        "-maxwarn 1000 -quiet".format(
            Config.md_configs.GroDIR, mdp, gro, top, ndx, sysname
        )
    )
    sb.run(tpr_cmd, shell=True, stdout=Config.md_configs.LOG, stderr=sb.STDOUT)

    mdrun_cmd = (
        "{0}/gmx mdrun -nt {1} -s {2}.tpr -x {2}.xtc -c {2}.gro -e {2}.edr "
        "-g {2}.log -o {2}.trr -rcon {3}".format(
            Config.md_configs.GroDIR,
            Config.md_configs.nCPUs,
            sysname,
            Config.md_configs.rcon,
        )
    )
    sb.run(mdrun_cmd, shell=True, stdout=Config.md_configs.LOG, stderr=sb.STDOUT)


def remove_freezed_velocities(top: str, relaxgro: str, effectivegro: str) -> None:
    freezed_atoms = top.index_atoms[Config.md_configs.titrating_group]

    final_gro = ""
    final_gro_footer = ""
    final_gro_atoms = {}

    # Get titrating atoms coords and velocities from EffectiveGro
    with open(effectivegro) as f_effective, open(relaxgro) as f_relax:
        line_counter = 0
        natoms_left = 0
        for line_eff, line_relax in zip(f_effective, f_relax):
            line_counter += 1
            if natoms_left > 0:
                natoms_left -= 1
                (aname, anumb, resname, resnumb, x, y, z) = read_gro_line(line_eff)
                if anumb in freezed_atoms:
                    final_gro_atoms[anumb] = line_eff
                else:
                    final_gro_atoms[anumb] = line_relax
            elif line_counter == 2:
                natoms_left = int(line_eff.strip())
                final_gro += line_eff
            elif line_counter == 1:
                final_gro += line_eff
            else:
                final_gro_footer += line_eff

    for atom in sorted(final_gro_atoms.keys()):
        final_gro += final_gro_atoms[atom]
    final_gro += final_gro_footer

    with open(relaxgro, "w") as f_new:
        f_new.write(final_gro)


def concat_results(
    runname: str,
    cycle: int,
    effective_name: str,
    init_time: float,
    simtime_begin: float,
) -> None:
    if cycle == 0:
        edr_in = "{}.edr".format(effective_name)
        edr_out = "{}.edr".format(runname)
        edr_stdin = "{}\nc".format(simtime_begin)

        xtc_in = "{}.xtc".format(effective_name)
        xtc_out = "{}.xtc".format(runname)
        xtc_cmd = "{}/gmx trjconv -f {} -o {} -t0 {} -s {}.tpr".format(
            Config.md_configs.GroDIR, xtc_in, xtc_out, init_time, effective_name
        )
        xtc_stdin = "System"
    else:
        edr_cur = "{}.edr".format(runname)
        edr_in = "{} {}.edr".format(edr_cur, effective_name)
        edr_out = "aux.edr"
        edr_stdin = "{}\n{}".format(init_time, simtime_begin)

        xtc_cur = "{}.xtc".format(runname)
        xtc_in = "{} {}.xtc".format(xtc_cur, effective_name)
        xtc_out = "aux.xtc"
        xtc_cmd = "{}/gmx trjcat -f {} -o {} -settime".format(
            Config.md_configs.GroDIR, xtc_in, xtc_out
        )
        xtc_stdin = "{}\nc\n".format(init_time)

    edr_cmd = "{}/gmx eneconv -o {} -f {} -settime".format(
        Config.md_configs.GroDIR, edr_out, edr_in
    )
    sb.run(
        edr_cmd,
        shell=True,
        input=str.encode(edr_stdin),
        stdout=Config.md_configs.LOG,
        stderr=sb.STDOUT,
    )

    sb.run(
        xtc_cmd,
        shell=True,
        input=str.encode(xtc_stdin),
        stdout=Config.md_configs.LOG,
        stderr=sb.STDOUT,
    )

    if cycle != 0:
        os.rename(edr_out, edr_cur)
        os.rename(xtc_out, xtc_cur)

    save_logs_cmd = "cat {}.log >> {}.log".format(effective_name, runname)
    sb.run(save_logs_cmd, shell=True)


def clean_dir(effective_name: str, relax_name: str, pypka_input: str) -> None:
    clean_cmd = (
        "rm -f {0}.edr {0}.log {0}.xtc {0}.tpr {1}.edr"
        " {1}.gro {1}.log {1}.tpr {1}.xtc state.cpt"
        " state_prev.cpt mdout.mdp {2} *#".format(
            effective_name, relax_name, pypka_input
        )
    )
    sb.run(clean_cmd, shell=True)


def final_cleanup(sysname: str, tmpDIR: str, effective_name: str) -> None:
    os.system("cp {0}.gro ../{1}.gro".format(effective_name, sysname))
    os.system("cp {0}.edr {0}.log {0}.mocc {0}.occ {0}.xtc {0}.top ../".format(sysname))
    os.chdir("../")
    os.system("rm -r {}".format(tmpDIR))


def write_pypka_inputgro(f_centered_name: str, nontitrating_res: str) -> None:
    new_lines = []
    with open(f_centered_name) as f:
        line_counter = 0
        natoms_left = 0
        for line in f:
            line_counter += 1
            if natoms_left > 0:
                natoms_left -= 1
                resnumb = read_gro_line(line)[3]

                if resnumb in nontitrating_res:
                    resname = nontitrating_res[resnumb]
                    line = "{}{}{}".format(line[:5], resname, line[8:])

            elif line_counter == 2:
                natoms_left = int(line.strip())

            new_lines.append(line)

    with open(Config.md_configs.pypka_input, "w") as f_new:
        f_new.write("".join(new_lines))


def run_pbmc(
    params: Dict, sites: List, pH: float, offset: int, fixed_sites: Dict
) -> Tuple[Dict[int, int], Dict[int, float], Dict[int, list]]:
    prot_states = {}
    prot_avgs = {}
    taut_probs = {}
    if len(sites) > 0:
        with redirect_stdout(io.StringIO()) as f:
            tit = Titration(params, sites={"A": sites}, fixed_sites={"A": fixed_sites})

        Config.md_configs.LOG.write(f.getvalue())
        print("\r" + " " * 90, end="\r")

        for site in tit:
            resnumb = site.getResNumber()
            resname = site.getName()
            if resname in ("NTR", "CTR"):
                resnumb += offset
            prot_states[resnumb] = site.getFinalState(pH)
            prot_avgs[resnumb] = site.getTitrationCurve()[pH]
            taut_probs[resnumb] = site.getTautsProb(pH)

    return (prot_states, prot_avgs, taut_probs)


def run_pkai(
    fname, sites, pH, offset
) -> Tuple[Dict[int, int], Dict[int, float], Dict[int, list]]:
    try:
        Pege
    except NameError:
        from pege import Pege

    prot_states = {}
    prot_avgs = {}
    taut_probs = {}

    protein = Pege(fname, save_final_pdb=False, fix_pdb=False)

    for site in sites:
        resnumb = site
        if isinstance(resnumb, str):
            resnumb = int(site[:-1]) + offset

        (
            prot_states[resnumb],
            prot_avgs[resnumb],
            taut_probs[resnumb],
        ) = protein.get_residue_taut_probs("A", site, pH)

    return (prot_states, prot_avgs, taut_probs)


def run_cphmd(top: str) -> None:
    cur_time, prev_time = None, None
    cycle_times = []
    for cycle in range(Config.md_configs.InitCycle, Config.md_configs.EndCycle):
        cur_time = time.time()

        simtime_begin, simtime_end = Config.md_configs.get_simtime(cycle)

        info = "{} | Starting CpHMD Cycle #{}".format(get_curtime(), cycle)
        logging.info(info)

        if prev_time:
            cycle_time = cur_time - prev_time
            cycle_times.append(cycle_time)
            cycles_left = Config.md_configs.EndCycle - cycle
            avg_time = sum(cycle_times) / len(cycle_times) * cycles_left
            end_guess = datetime.now() + timedelta(seconds=avg_time)
            info = "{} | Estimated Time of Completion: {}".format(
                get_curtime(), end_guess.strftime("%d %B %Y %H:%M")
            )
            logging.info(info)

        ### PB/MC ###
        f_centered_name = center_titrable_molecule(
            Config.md_configs.titrating_group,
            Config.md_configs.GRO,
            Config.md_configs.NDX,
        )

        if not Config.md_configs.pkai:
            info = "{} | PB/MC - titrating {} sites".format(
                get_curtime(), len(top.titrating_sites)
            )
            if top.fixed_sites:
                info += " and {} fixed states".format(len(top.fixed_sites))
            logging.info(info)

            write_pypka_inputgro(f_centered_name, top.nontit_tautomers)
            os.system(f"cp {Config.md_configs.pypka_input} pypka_input_{cycle}.gro")

            pb_prot_states, pb_prot_avgs, pb_taut_probs = run_pbmc(
                Config.md_configs.pypka_params,
                top.titrating_sites[:],
                Config.md_configs.pH,
                top.offset,
                {k: v["state"] for k, v in top.fixed_sites.items()},
            )
        else:
            info = "{} | AI/MC - titrating {} sites".format(
                get_curtime(), len(top.titrating_sites)
            )
            logging.info(info)

            pb_prot_states, pb_prot_avgs, pb_taut_probs = run_pkai(
                f_centered_name,
                top.titrating_sites[:],
                Config.md_configs.pH,
                top.offset,
            )

        top.update(pb_prot_states, pb_prot_avgs, pb_taut_probs)
        top.write_top_file()

        ### MD ###
        info = "{} | MD    - solvent relaxation for {} steps".format(
            get_curtime(), Config.md_configs.RelaxSteps
        )
        logging.info(info)
        os.system(f"cp {Config.md_configs.TOP} topol_{cycle}.top")

        run_dynamics(
            Config.md_configs.MDP_relax,
            Config.md_configs.GRO,
            Config.md_configs.TOP,
            Config.md_configs.NDX,
            Config.md_configs.relax_name,
        )

        remove_freezed_velocities(
            top, Config.md_configs.GRO_relax, Config.md_configs.GRO
        )

        info = "{} | MD    - production of {} steps".format(
            get_curtime(), Config.md_configs.EffectiveSteps
        )
        logging.info(info)

        run_dynamics(
            Config.md_configs.MDP,
            Config.md_configs.GRO_relax,
            Config.md_configs.TOP,
            Config.md_configs.NDX,
            Config.md_configs.effective_name,
        )

        ### Concatenate Results ###
        concat_results(
            Config.md_configs.sysname,
            cycle,
            Config.md_configs.effective_name,
            Config.md_configs.InitTime,
            simtime_begin,
        )
        clean_dir(
            Config.md_configs.effective_name,
            Config.md_configs.relax_name,
            Config.md_configs.pypka_input,
        )

        prev_time = cur_time
