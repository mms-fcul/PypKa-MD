import os
from datetime import datetime

from pypka import getTitrableSites
from pdbmender.formats import gro2pdb, read_pdb_line


def get_curtime():
    return datetime.today().strftime("%Y-%m-%d %H:%M")


def remove_comments(line, comment_delimiter=";"):
    if len(comment_delimiter) == 1:
        return line.split(comment_delimiter)[0]
    return remove_comments(line.split(comment_delimiter[1:])[0])


def create_link(path):
    os.system("ln -s -f {} .".format(path))


def get_titrable_sites(groname):
    gro2pdb(groname, "tmp.pdb", save_box=False)
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


def check_convert_termini(site, offset):
    if isinstance(site, str) and site[-1] in "NC":
        site = offset + int(site[:-1])
    return site
