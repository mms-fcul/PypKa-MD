import logging
import random
from typing import Optional, List, Tuple, Dict

from pypkamd.misc import (
    check_convert_termini,
    get_titrable_sites,
    remove_comments,
    read_ff_dict_section,
)
from pypkamd.constants import TITRABLE_AAS
from pypkamd.configs import Config


class Topology:
    def __init__(self, configs: Config):

        self.sections = ["atoms", "bonds", "pairs", "angles", "dihedrals", "impropers"]
        self.section = {section: [] for section in self.sections}
        self.section["others"] = {}

        self.tit_atoms = {}
        self.tit_atypes = {}
        self.tit_charges = {}
        self.tit_bonds = {}
        self.tit_angles = {}
        self.tit_dihedrals = {}
        self.tit_impropers = {}

        self.prop_types = {
            "t": (1, self.tit_atypes),
            "q": (1, self.tit_charges),
            "b": (2, self.tit_bonds, self.section["bonds"]),
            "a": (3, self.tit_angles, self.section["angles"]),
            "d": (4, self.tit_dihedrals, self.section["dihedrals"]),
            "i": (4, self.tit_impropers, self.section["impropers"]),
        }

        self.offset = None
        self.ters = {}

        self.config_vars = configs

        if configs.sites == ["all"]:
            configs.sites = get_titrable_sites(configs.GROin)

            logging.info(
                " {} titrable sites have been identified".format(len(configs.sites))
            )

        self.all_sites = configs.sites[:]
        self.titrating_sites = configs.sites[:]
        self.rt_cur_cycle = 0
        self.fixed_sites = {}

        self.read_input_top(
            configs.TOPin, configs.sites, configs.ffDIR.split("/")[-1], configs.ff_dict
        )
        self.read_ff_dict(configs.ff_dict)

        self.nontit_tautomers = {}
        self.get_non_tit_tauts(configs.ff_tautomers)

        self.index_atoms = {configs.titrating_group: []}

        self.titres_states = {}
        self.mocc = {}
        self.occ = {}
        for resnumb in configs.sites:
            if type(resnumb) == str:
                resnumb = self.offset + int(resnumb[:-1])

            self.titres_states[resnumb] = []
            self.mocc[resnumb] = []
            self.occ[resnumb] = []

        self.f_output = configs.TOP
        self.f_occ = configs.sysname + ".occ"
        self.f_mocc = configs.sysname + ".mocc"

    def get_nres(self) -> int:
        return int(self.section["atoms"][-1][2])

    def get_atom(self, anumb: str) -> Optional[int]:
        for atom in self.section["atoms"]:
            if atom[0] == anumb:
                return atom
        return None

    def get_titrable_resnames(self) -> List[str]:
        resnames = []
        for resnumb, atom_i_list in self.tit_atoms.items():
            if resnumb < self.offset:
                atom_i = atom_i_list[0]
                atom = self.section["atoms"][atom_i]
                resname = atom[3]
                resnames.append(resname)
            else:
                ter_numb = resnumb - self.offset
                resname = self.ters[ter_numb]
                resnames.append(resname)

        return list(set(resnames))

    def atom_is_titrable(self, anumb: int) -> Optional[int]:
        anumb = int(anumb) - 1
        for residue, atoms in self.tit_atoms.items():
            if anumb in atoms:
                return residue
        return False

    def get_atoms(self, resname: str, aname: str) -> Tuple[int, str, str]:
        # print("###################################", resname, aname)
        for atom in self.section["atoms"]:
            anumb, resnumb, resname_, aname_ = atom[0], int(atom[2]), atom[3], atom[4]
            tit_residue = self.atom_is_titrable(anumb)
            if (
                resnumb in self.ters.keys()
                and int(anumb) in self.tit_atoms[resnumb + self.offset]
            ):
                resname_ = self.ters[resnumb]
            if aname == aname_ and tit_residue and resname == resname_:
                yield (anumb, resname, aname)

    def search_by_anames(
        self, search_item: str, resname: str, anames: str
    ) -> Tuple[int, list]:
        anumbs = []
        for aname in anames:
            for i in self.get_atoms(resname, aname):
                anumbs.append(i[0])

        for i, item in enumerate(search_item):
            item_anames = []

            for anumb in item[: len(anames)]:
                atom = self.get_atom(anumb)
                anumb, resname_, aname_ = atom[0], atom[3], atom[4]
                if anumb in anumbs:
                    item_anames.append(aname_)

            if item_anames == anames:
                yield (i, item)

    def get_resnumb_from_anumb(self, anumb: str) -> int:
        for atom in self.section["atoms"]:
            anumb_ = atom[0]
            if anumb_ == anumb:
                return atom[2]

    def read_input_top(
        self, f_top: str, titrating_resnumbs: List[int], ffDIR: str, ff_dict: str
    ):
        def save_other_lines(last_section: str, line: str):
            if last_section == None:
                last_section = "begin"
            elif last_section == False:
                last_section = "end"
            if last_section not in self.section["others"]:
                self.section["others"][last_section] = []
            self.section["others"][last_section].append(line)

        ter_resnumbs = {}
        tit_resnumbs = []
        for i in titrating_resnumbs:
            if type(i) == str:
                numb = int(i.strip("NC"))
                ter_resnumbs[numb] = i[-1]
                self.ters[numb] = "NTR"
            else:
                tit_resnumbs.append(i)

        cph_ready_resnames = {i: ii for i, ii in TITRABLE_AAS.items()}

        termini_atoms = {}
        for line in read_ff_dict_section(ff_dict, "replace"):
            cols = line.split()
            resname = cols[0]
            atom_type = cols[2]

            if resname[:2] not in ("NT", "CT"):
                continue

            if resname not in termini_atoms:
                termini_atoms[resname] = []
            if atom_type not in termini_atoms[resname]:
                termini_atoms[resname].append(atom_type)

        last_section = None
        section = "begin"
        tmp_ter_atoms = {}
        with open(f_top) as f:
            for line in f:
                clean_line = remove_comments(line, comment_delimiter=";").strip()
                if not clean_line:
                    if not line.strip() and last_section != None:
                        last_section, section = section, False
                    continue

                if "#include" in line and "./" in line:
                    ff = line.split("./")[1].split("/")[0]
                    assert (
                        ff == ffDIR
                    ), "Topology file error\nLine {}\nShould be replaced by: {}".format(
                        line, line.replace(ff, ffDIR)
                    )

                if clean_line.startswith("[ "):
                    tmp_section = clean_line.strip("[] \n")
                    if self.section["dihedrals"] and "dihedrals" == tmp_section:
                        tmp_section = "impropers"

                    if tmp_section in self.sections:
                        last_section, section = section, tmp_section

                    if section not in self.sections:
                        save_other_lines(last_section, clean_line)

                elif section in self.sections:
                    cols = clean_line.split()
                    self.section[section].append(cols)

                    resnumb = int(cols[2])
                    if section == "atoms":
                        resname = cols[3]
                        atype = cols[4]
                        atom_i = len(self.section["atoms"]) - 1
                        add_ter_atom = False

                        if resnumb in tit_resnumbs:
                            if resname[:2] in cph_ready_resnames.keys():
                                self.section["atoms"][-1][3] = cph_ready_resnames[
                                    resname[:2]
                                ]

                        if resnumb in ter_resnumbs:
                            ter_type = ter_resnumbs[resnumb]
                            if ter_type == "N":
                                if (
                                    resname not in ("GLY", "PRO")
                                    and atype in termini_atoms["NT"]
                                ):
                                    self.ters[resnumb] = "NT"
                                    add_ter_atom = True
                                    # self.section["atoms"][-1][3] = "NT"

                                elif (
                                    resname in ("GLY", "PRO")
                                    and atype
                                    in termini_atoms["NTGLY"] + termini_atoms["NTPRO"]
                                ):
                                    self.ters[resnumb] = "NT" + resname
                                    add_ter_atom = True
                            elif ter_type == "C":
                                if atype in termini_atoms["CT"]:
                                    self.ters[resnumb] = "CT"
                                    add_ter_atom = True
                                    # self.section["atoms"][-1][3] = "CT"

                            if add_ter_atom:
                                if resnumb not in tmp_ter_atoms:
                                    tmp_ter_atoms[resnumb] = []
                                tmp_ter_atoms[resnumb].append(atom_i)

                        if (
                            resnumb in tit_resnumbs
                            and atype not in ("N", "H", "CA", "C", "O")
                            and not add_ter_atom
                        ):
                            if resnumb not in self.tit_atoms:
                                self.tit_atoms[resnumb] = []
                            self.tit_atoms[resnumb].append(atom_i)

                else:
                    save_other_lines(last_section, clean_line)

        self.offset = 10 ** (len(str(self.get_nres())) + 1)

        for resnumb, atoms in tmp_ter_atoms.items():
            resnumb += self.offset
            self.tit_atoms[resnumb] = atoms

        # from pprint import pprint

        # pprint(self.tit_atoms)
        # pprint(self.tit_angles)
        # pprint(self.section["atoms"])
        # pprint(self.ters)
        # exit()

    def read_ff_dict(self, ff_dict: dict) -> None:

        tit_res_cph = self.get_titrable_resnames()

        for line in read_ff_dict_section(ff_dict, "replace"):
            cols = line.split()
            resname = cols[0]
            line_type = cols[1]
            natoms = self.prop_types[line_type][0]
            n_tauts = int(cols[2 + natoms])
            assert n_tauts + natoms + 3 == len(
                cols
            ), "There is an inconsistency in the ff_dict file"

            if resname not in tit_res_cph:
                # print(resname, tit_res)
                continue

            properties = cols[-n_tauts:]
            prop_obj = self.prop_types[line_type][1]

            if line_type in ("t", "q"):
                aname = cols[2]
                for atom in self.get_atoms(resname, aname):
                    prop_obj[atom[0]] = properties

            elif line_type in ("b", "a", "d", "i"):
                anames = cols[2 : 2 + natoms]
                search_in = self.prop_types[line_type][2]

                for i, item in self.search_by_anames(search_in, resname, anames):
                    resnumb = self.get_resnumb_from_anumb(item[0])
                    prop_obj[i] = (resnumb, item[:natoms], properties)

        # from pprint import pprint
        # pprint(self.tit_angles)
        # pprint(self.prop_types["a"])
        # pprint(self.tit_charges)
        # exit()

    def get_non_tit_tauts(self, titratable_residues) -> None:
        anames = {
            resname: [i[1] for i in cols]
            for resname, cols in titratable_residues.items()
        }
        for atom in self.section["atoms"]:
            _, _, resnumb, resname, aname, _, charge, _ = atom

            if (
                resname in titratable_residues
                and aname in anames[resname]
                and int(resnumb) not in self.all_sites
            ):
                taut_numb = None
                for taut in titratable_residues[resname]:
                    taut_numb, taut_aname, taut_acharge = taut
                    if taut_aname == aname and taut_acharge == charge:
                        break
                self.nontit_tautomers[int(resnumb)] = taut_numb

    def read_index(self, ndx: str) -> None:
        index_groups = self.index_atoms.keys()
        with open(ndx) as f:
            save_group = False
            for line in f:
                if "[" in line and "]" in line:
                    index_group = line.strip("[ ]\n")
                    if index_group in index_groups:
                        save_group = True
                    else:
                        save_group = False
                elif save_group:
                    try:
                        indexes = [int(i) for i in line.strip().split()]
                        self.index_atoms[index_group] += indexes
                    except:
                        msg = "Index group {} in file {} " "is not valid.".format(
                            index_group, ndx
                        )
                        raise IOError(msg)

    def sample_fixed(self, prob_states: Dict[int, int]) -> Dict[int, int]:
        for site, fixed_details in self.fixed_sites.items():
            tauts = list(range(len(fixed_details["taut_probs"])))
            new_state = random.choices(tauts, fixed_details["taut_probs"])[0]

            if isinstance(site, str) and site[-1] in "NC":
                s = int(site[:-1]) + self.offset
            else:
                s = site

            prob_states[s] = new_state
            self.fixed_sites[site]["state"] = new_state
            self.mocc[s].append(fixed_details["prob"])

        return prob_states

    def update(self, new_states, new_protavg, new_tautprobs):
        if self.config_vars["reduced_titration"]:
            self.rt_cur_cycle += 1

            if self.fixed_sites:
                new_states = self.sample_fixed(new_states)

            if self.config_vars["rt_cycles"] == self.rt_cur_cycle:
                self.rt_cur_cycle = 0
                self.titrating_sites = self.config_vars["sites"][:]
                self.fixed_sites = {}

        # Save the current average protonation
        for site, avg in new_protavg.items():
            self.mocc[site].append(avg)

            if self.config_vars["reduced_titration"] and self.rt_cur_cycle == 1:
                if (
                    avg < self.config_vars["rt_limit"]
                    or avg > 1 - self.config_vars["rt_limit"]
                ):
                    if site > self.offset:
                        s = site - self.offset
                        if "NT" in self.ters[s]:
                            s = "{}N".format(s)
                        elif "CT" in self.ters[s]:
                            s = "{}C".format(s)
                        else:
                            raise Exception("Something is wrong!")
                    else:
                        s = site

                    site_i = self.titrating_sites.index(s)
                    del self.titrating_sites[site_i]

                    self.fixed_sites[s] = {
                        "prob": avg,
                        "state": new_states[site],
                        "taut_probs": new_tautprobs[site],
                    }

        # print(self.titres_states[site])
        # print("rt cycle:", self.rt_cur_cycle)
        # print("fixed:", self.fixed_sites)
        # print("titrating next:", self.titrating_sites)

        # Deal with the new protonation states
        for site, state in new_states.items():
            self.titres_states[site] = state
            self.occ[site].append(state)

            for atom_id in self.tit_atoms[site]:
                atom = self.section["atoms"][atom_id]
                anumb = atom[0]

                if anumb in self.tit_atypes:
                    atom[1] = self.tit_atypes[anumb][state]

                if anumb in self.tit_charges:
                    atom[6] = self.tit_charges[anumb][state]

        for prop_code, prop_objs in self.prop_types.items():
            if prop_code in ("q", "t"):
                continue
            natoms, titprop_list, top_list = prop_objs
            for index, item in titprop_list.items():
                resnumb, anumbs = int(item[0]), item[1]

                assert top_list[index][:natoms] == anumbs

                if resnumb in self.ters:
                    resnumb += self.offset

                new_state = self.titres_states[resnumb]
                new_prop = item[2][new_state]

                top_list[index][-1] = new_prop

        self.write_prot_file(self.occ, self.f_occ, "occ")
        self.write_prot_file(self.mocc, self.f_mocc, "mocc")

    def get_ordered_sites(self):
        def get_resnumber(res):
            if isinstance(res, str):
                res = int(res[:-1])
            return res

        sorted_sites = []
        termini = []
        for res in self.all_sites:
            if isinstance(res, str):
                termini.append(res)
            else:
                sorted_sites.append(res)
        sorted_sites.sort()
        for res in termini:
            res_numb = int(res[:-1])
            if not sorted_sites:
                sorted_sites.append(res)
                continue
            sorted_res = get_resnumber(sorted_sites[0])

            i = 0
            while sorted_res < res_numb and i + 2 <= len(sorted_sites):
                i += 1
                sorted_res = get_resnumber(sorted_sites[i])
                # print(sorted_res, res_numb, i, len(sorted_sites))

            if res_numb > sorted_res:
                i += 1

            if res[-1] == "C" and sorted_res <= res_numb:
                i += 1

            sorted_sites.insert(i, res)
        return sorted_sites

    def write_prot_file(self, prot_obj, fprot_name, mode):
        content = "# "
        nframes = len(list(prot_obj.values())[0])

        for site in self.get_ordered_sites():
            if mode == "mocc":
                content += "{:>8} ".format(site)
            elif mode == "occ":
                content += "{:>4} ".format(site)

        content += "\n  "
        for i in range(nframes):
            for site in self.get_ordered_sites():
                site = check_convert_termini(site, self.offset)
                site_occ = prot_obj[site][i]
                if mode == "mocc":
                    content += "{:8.6f} ".format(site_occ)
                elif mode == "occ":
                    content += "{:4} ".format(site_occ)
            content += "\n  "
        with open(fprot_name, "w") as f_new:
            f_new.write(content)

    def write_top_file(self):
        def seclist_to_str(seclist, extra_spaces=5):
            sec_str = ""
            maxsize = len(self.section["atoms"][-1][0]) + extra_spaces
            for line in seclist:
                tmp_line = ""
                for item in line[:-1]:
                    col = "{}".format(item)
                    while len(col) < maxsize:
                        col = " " + col
                    tmp_line += col
                tmp_line += " " * maxsize
                tmp_line += line[-1]
                sec_str += tmp_line + "\n"
            return sec_str

        top_content = ""
        if "begin" in self.section["others"]:
            top_content += "\n".join(self.section["others"]["begin"]) + "\n"

        top_content += "[ atoms ]\n"
        top_content += seclist_to_str(self.section["atoms"], extra_spaces=8) + "\n"

        # print(seclist_to_str(self.section["atoms"], extra_spaces=8))

        for section in self.sections[1:]:
            section_name = section
            if section == "impropers":
                section_name = "dihedrals"

            top_content += "[ {} ]\n".format(section_name)
            top_content += seclist_to_str(self.section[section])
            if section in self.section["others"]:
                top_content += "\n".join(self.section["others"][section]) + "\n"

        if "end" in self.section["others"]:
            top_content += "\n".join(self.section["others"]["end"]) + "\n"

        final_content = ""
        for line in top_content.splitlines():
            if line.startswith("[ ") or line.startswith("#ifdef"):
                line = "\n" + line
            final_content += line + "\n"

        with open(self.f_output, "w") as ftop_new:
            ftop_new.write(final_content)
