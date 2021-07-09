import pytest
import sys

sys.path.insert(1, "../")
from pypkamd.cphmd import run_pbmc


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
