import pytest
import sys

sys.path.insert(1, "../")
from pypkamd.cphmd import run_pbmc
from pypkamd.topology import Topology


def read_mocc(fname):
    mocc = []
    with open(fname) as f:
        for line in f:
            values = line.split()
            line_mocc = []
            for value in values:
                value = float(value)
                line_mocc.append(value)
            mocc.append(line_mocc)
    return mocc


class TestPypKa(object):
    default_params = {
        "epsin": 2,
        "ionicstr": 0.1,
        "pbc_dimensions": 0,
        "ncpus": -1,
        "convergence": 0.01,
        "clean_pdb": False,
        "CpHMD_mode": True,
        "nlit": 100,
        "nonit": 5,
    }

    def test_ntr(self):
        offset = 100
        params = self.default_params
        sites = ["1N"]
        pH = 7.2
        params["pH"] = str(pH)
        base_dir = "penta_control/NTR"

        mocc = read_mocc(f"{base_dir}/old_mocc")
        mocc_line = -1
        diffs = []
        for rep in range(1, 4):
            for i in range(50):
                mocc_line += 1
                params["structure"] = f"{base_dir}/{rep}_protein_{i:0>3}.gro"
                print(params["structure"])
                (_, prot_avgs, _) = run_pbmc(params, sites, pH, offset, {})
                diff = prot_avgs[101] - mocc[mocc_line]
                assert (
                    abs(diff) < 0.005
                ), f"{params['structure']} Difference is not small enough: {diff}"
                diffs.append(str(diff[0]))

        with open("tmp_NTR.out", "w") as f:
            f.write("\n".join(diffs))

    def test_tetra(self):
        offset = 100
        params = self.default_params
        sites = ["1N", 2, 3, "4C"]
        pH = 4.0
        params["pH"] = str(pH)

        mocc = read_mocc("tetra/AHHA/old_mocc")
        diffs = []
        mocc_line = -1
        for rep in range(1, 4):
            for i in range(50):
                mocc_line += 1
                params["structure"] = f"tetra/AHHA/{rep}_protein_{i:0>3}.gro"
                print(params["structure"])
                _, prot_avgs, _ = run_pbmc(params, sites, pH, offset, {})

                for resi, res in enumerate((101, 2, 3, 104)):
                    mocc_prot = mocc[mocc_line][resi]
                    diff = prot_avgs[res] - mocc_prot
                    diffs.append(str(diff))
                    assert diff < 0.06, (
                        f"Difference is not small enough for residue {res}: {round(diff, 3)}"
                        f"; Expected: {mocc_prot}; Obtained: {prot_avgs[res]}"
                    )
                    if diff > 0.03:
                        with open("worst_AHHA.out", "a") as f:
                            f.write(
                                f"AHHA {rep}_protein_{i:0>3}.gro {mocc_prot} {prot_avgs[res]} {round(diff, 3)} \n"
                            )
        with open("tmp_AHHA.out", "w") as f:
            f.write("\n".join(diffs))


class TestTop(object):
    def test_get_ordered_sites(self):
        from collections import namedtuple

        def perf_test(sites, result):
            configs = namedtuple(
                "configs",
                (
                    "sites",
                    "TOPin",
                    "ffID",
                    "ff_dict",
                    "titrating_group",
                    "TOP",
                    "sysname",
                ),
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
            assert top.get_ordered_sites() == result

        perf_test(["1N"], ["1N"])
        perf_test(["1N", 3], ["1N", 3])
        perf_test(["1N", 0], [0, "1N"])
        perf_test([2, "1C", "3N", "3C", "2N"], ["1C", "2N", 2, "3N", "3C"])
        perf_test(["1C", 1, "3N", "3C", "2N"], [1, "1C", "2N", "3N", "3C"])
        perf_test(["1C", "3N", "3C", "2N"], ["1C", "2N", "3N", "3C"])
        perf_test([1, "1N", 2, 10, "30C", 25], ["1N", 1, 2, 10, 25, "30C"])

        perf_test(
            [1, 2, 10, "30C", 25, 45, "10N", 55, 11, 88],
            [1, 2, "10N", 10, 11, 25, "30C", 45, 55, 88],
        )
        perf_test(
            [
                "88C",
                1,
                2,
                10,
                "30C",
                25,
                45,
                "10N",
                55,
                11,
                88,
                "2N",
                "55C",
                "55N",
                "1N",
            ],
            [
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
            ],
        )
