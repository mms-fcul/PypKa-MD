import io
import logging
import os
import subprocess as sb
import time
from contextlib import redirect_stdout
from datetime import datetime, timedelta

from pypka import Titration
from pdbmender.formats import read_gro_line

from pypkamd.configs import Config
from pypkamd.misc import create_link, get_curtime, remove_comments


def create_cphdm_directory():
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
    os.system("cp -f {} {}/.".format(Config.md_configs.DATin, Config.md_configs.tmpDIR))

    os.chdir("{}".format(Config.md_configs.tmpDIR))

    create_link("{}".format(Config.md_configs.ffDIR))
    create_link("{}/residuetypes.dat".format(Config.md_configs.ffDIR))

    # convert ionic strength from molar to molecule/nm^3
    # ionicstrMolecule = ionicstr * 0.6022

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
        for line in f:
            cleaned_line = remove_comments(line)
            parts = cleaned_line.split("=")
            if len(parts) >= 2:
                param = parts[0].strip()
                if param not in relax_mdp_new_params:
                    relax_mdp_content += line
            new_mdp_content += line

        relax_mdp_content += """
;Solvent relaxation parameters

nsteps = {0}
Pcoupl = No
;constraints = none
;energygrps = {1}
freezegrps = {1}
freezedim = Y Y Y
energygrp_excl = {1} {1}
        """.format(
            Config.md_configs.RelaxSteps, Config.md_configs.titrating_group
        )

        fmdp_new.write(new_mdp_content)
        frelax_new.write(relax_mdp_content)
    os.remove(Config.md_configs.MDPin)

    # Create a Index File
    if Config.md_configs.NDXin:
        os.system("cp ../{} {}".format(Config.md_configs.NDXin, Config.md_configs.NDX))
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

    # Create a Tpr File to be used in the centering procedure
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


def center_titrable_molecule(titrating_group, gro, ndx, outfile):
    # Centering procedure
    # effective.gro -> fixgro_input.gro -> pypka_input.gro

    # TODO change input variable fixgr

    fixgro_cmd = (
        "{}/gmx trjconv -f {} "
        "-o fixgro_input.gro -s fixgro.tpr -n {} "
        "-pbc mol -quiet; ".format(Config.md_configs.GroDIR, gro, ndx)
    )
    fixgro_cmd += "{}/fixbox fixgro_input.gro {} " "> {}".format(
        Config.md_configs.fixboxDIR, Config.md_configs.DATin, outfile
    )

    input_str = "{}\n".format(titrating_group)
    sb.run(
        fixgro_cmd,
        shell=True,
        input=str.encode(input_str),
        stdout=Config.md_configs.LOG,
        stderr=sb.STDOUT,
    )

    rm_cmd = "rm fixgro_input.gro"
    sb.run(rm_cmd, shell=True)


def run_dynamics(mdp, gro, top, ndx, sysname):
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


def remove_freezed_velocities(top, relaxgro, effectivegro):
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


def concat_results(runname, cycle, effective_name, init_time, simtime_begin):
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


def clean_dir(effective_name, relax_name, pypka_input):
    clean_cmd = (
        "rm -f {0}.edr {0}.log {0}.xtc {0}.tpr {1}.edr"
        " {1}.gro {1}.log {1}.tpr {1}.xtc state.cpt"
        " state_prev.cpt mdout.mdp {2} *#".format(
            effective_name, relax_name, pypka_input
        )
    )
    sb.run(clean_cmd, shell=True)


def final_cleanup(sysname, tmpDIR, effective_name):
    os.system("cp {0}.gro ../{1}.gro".format(effective_name, sysname))
    os.system("cp {0}.edr {0}.log {0}.mocc {0}.occ {0}.xtc {0}.top ../".format(sysname))
    os.chdir("../")
    os.system("rm -r {}".format(tmpDIR))


def run_pbmc(params, sites, pH, offset, fixed_sites):
    # print(sites, fixed_sites)
    prot_states = {}
    prot_avgs = {}
    taut_probs = {}
    if len(sites) > 0:

        tit = Titration(params, sites={"A": sites}, fixed_sites={"A": fixed_sites})

        for site in tit:
            resnumb = site.getResNumber()
            resname = site.getName()
            if resname in ("NTR", "CTR"):
                resnumb += offset
            prot_states[resnumb] = site.getFinalState(pH)
            prot_avgs[resnumb] = site.getTitrationCurve()[pH]
            taut_probs[resnumb] = site.getTautsProb(pH)

            # print(site.final_states, site.states_prob)
            # print(resname, resnumb, pb_prot_states[resnumb], pb_prot_avgs[resnumb])
    # print(prot_avgs)

    return prot_states, prot_avgs, taut_probs


def run_cphmd(top):
    cur_time, prev_time = None, None
    cycle_times = []
    for cycle in range(Config.md_configs.InitCycle, Config.md_configs.EndCycle):
        cur_time = time.time()

        simtime_begin, simtime_end = Config.md_configs.get_simtime(cycle)

        info = "\n{} | Starting CpHMD Cycle #{}".format(get_curtime(), cycle)
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
        info = "{} | PB/MC - titrating {} sites".format(
            get_curtime(), len(top.titrating_sites)
        )
        if top.fixed_sites:
            info += " and {} fixed states".format(len(top.fixed_sites))
        logging.info(info)

        center_titrable_molecule(
            Config.md_configs.titrating_group,
            Config.md_configs.GRO,
            Config.md_configs.NDX,
            Config.md_configs.pypka_input,
        )

        with redirect_stdout(io.StringIO()) as f:
            pb_prot_states, pb_prot_avgs, pb_taut_probs = run_pbmc(
                Config.md_configs.pypka_params,
                top.titrating_sites[:],
                Config.md_configs.pH,
                top.offset,
                {k: v[1] for k, v in top.fixed_sites.items()},
            )
        print("\r" + " " * 90, end="\r")
        Config.md_configs.LOG.write(f.getvalue())

        top.update(pb_prot_states, pb_prot_avgs, pb_taut_probs)
        top.write_top_file()

        ### MD ###
        info = "{} | MD    - solvent relaxation for {} steps".format(
            get_curtime(), Config.md_configs.RelaxSteps
        )
        logging.info(info)

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
