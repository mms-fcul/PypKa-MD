import pytest
import sys

sys.path.insert(1, "../")
from pypkamd.cphmd import run_pbmc
from pypkamd.topology import Topology


def read_mocc(fname):
    mocc = []
    with open(fname) as f:
        for line in f:
            value = float(line.strip())
            mocc.append(value)
    return mocc


class TestPypKa(object):
    default_params = {
        "pH": "3.7",
        "epsin": 2,
        "ionicstr": 0.1,
        "pbc_dimensions": 0,
        "ncpus": 4,
        "convergence": 0.01,
        "clean_pdb": False,
        "CpHMD_mode": True,
    }

    def test_asp(self):
        offset = 100
        params = self.default_params
        sites = [4]
        pH = 3.7

        mocc = read_mocc("penta_control/ASP/old_mocc")
        for i in range(50):
            params["structure"] = f"penta_control/ASP/protein_{i:0>3}.gro"
            print(params["structure"])
            _, prot_avgs = run_pbmc(params, sites, pH, offset)
            diff = prot_avgs[4] - mocc[i]
            assert diff < 0.01, f"{i} Difference is not small enough."


class TestTop(object):
    def test_get_ordered_sites(self):
        from collections import namedtuple

        sites = (1, "1N", 2, 10, "30C", 25)
        configs = namedtuple(
            "configs",
            ("sites", "TOPin", "ffID", "ff_dict", "titrating_group", "TOP", "sysname"),
        )
        config = configs(
            sites,
            "lyso/MYLYSO_001.top",
            "G54a7pHt.ff",
            "../pypkamd/G54a7pHt.ff/protstates.dic",
            None,
            None,
            "test",
        )
        top = Topology(config)
        assert top.get_ordered_sites() == ["1N", 1, 2, 10, 25, "30C"]

        sites = (1, 2, 10, "30C", 25, 45, "10N", 55, 11, 88)
        config = configs(
            sites,
            "lyso/MYLYSO_001.top",
            "G54a7pHt.ff",
            "../pypkamd/G54a7pHt.ff/protstates.dic",
            None,
            None,
            "test",
        )
        top = Topology(config)
        assert top.get_ordered_sites() == [1, 2, "10N", 10, 11, 25, "30C", 45, 55, 88]

        sites = ("88C", 1, 2, 10, "30C", 25, 45, "10N", 55, 11, 88, "2N", "55C", "55N", "1N")
        config = configs(
            sites,
            "lyso/MYLYSO_001.top",
            "G54a7pHt.ff",
            "../pypkamd/G54a7pHt.ff/protstates.dic",
            None,
            None,
            "test",
        )
        top = Topology(config)
        assert top.get_ordered_sites() == [
            "1N",
            1,
            "2N",
            2,
            "10N",
            10,
            11,
            25,
            "30C",
            45,
            "55N",
            55,
            "55C",
            88,
            "88C",
        ]
