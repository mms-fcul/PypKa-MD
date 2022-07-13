import os
from datetime import datetime
from typing import List

from pypka import getTitrableSites
from pdbmender.formats import gro2pdb, read_pdb_line


def get_curtime() -> str:
    return datetime.today().strftime("%Y-%m-%d %H:%M")


def remove_comments(line: str, comment_delimiter: str = ";") -> str:
    if len(comment_delimiter) == 1:
        return line.split(comment_delimiter)[0]
    return remove_comments(line.split(comment_delimiter[1:])[0])


def create_link(path: str) -> None:
    os.system("ln -s -f {} .".format(path))


def get_titrable_sites(groname: str) -> List[int]:
    gro2pdb(groname, "tmp.pdb")
    sites, _ = getTitrableSites("tmp.pdb", ser_thr_titration=False)

    potential = []
    confirmed = []
    with open("tmp.pdb") as f:
        for line in f:
            if line.startswith("ATOM "):
                atom, _, res, chain, resnumb, *_ = read_pdb_line(line)
                if resnumb in sites["A"] and res == "CYS":
                    if resnumb not in potential:
                        potential.append(resnumb)
                    if atom == "HG" and resnumb not in confirmed:
                        confirmed.append(resnumb)
    to_del = set(potential) - set(confirmed)
    sites = set(sites["A"]) - to_del

    return list(sites)


def check_convert_termini(site: str, offset: int) -> int:
    if isinstance(site, str) and site[-1] in "NC":
        site = offset + int(site[:-1])
    return site

def read_ff_dict_section(fname, section):
    trigger_read = False
    with open(fname) as f_dict:
        for line in f_dict:
            if "[" in line and "]" in line:
                if section in line:
                    trigger_read = True
                else:
                    trigger_read = False
            elif trigger_read:
                yield line