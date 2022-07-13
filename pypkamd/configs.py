import os
import pwd
from builtins import bool
import logging
from pypka import __version__ as pypka_version
from pypka.config import ParametersDict as ParametersDictPypKa
from pypkamd.misc import remove_comments, read_ff_dict_section


def get_username():
    return pwd.getpwuid(os.getuid())[0]


class Config:
    @classmethod
    def storeParams(cls, mdp: str):
        if pypka_version < "2.3.0":
            raise Exception(
                "Please update your PypKa installation. \nCurrent version: {}\nRequired version: >= 2.3.0".format(
                    pypka_version
                )
            )

        cls.md_configs = MDConfig(mdp)
        cls.md_configs.read_input_mdp()


class ParametersDict(ParametersDictPypKa):
    input_conversion = {}
    input_conversion["dt"] = "TimeStep"
    input_conversion["tinit"] = "InitTime"
    input_conversion["ref_t"] = "temp"
    input_conversion["nsteps"] = "EffectiveSteps"


class MDConfig(ParametersDict):
    def __init__(self, fmdp: str):

        self.TimeStep = None
        self.InitTime = None
        self.EffectiveSteps = None
        self.RelaxSteps = 100

        self.nCycles = 0

        self.InitCycle = 0

        self.temp = None

        self.tmpDIR = "./tmp"
        self.script_dir = os.path.dirname(os.path.abspath(__file__))

        self.ff_options = {"CHARMM36mpH": "CHARMM36m", "G54A7pH": "G54A7"}
        self.ffID = None

        self.USER = get_username()
        self.HOST = os.uname()[1]

        self.fixboxDIR = "{0}/fixbox/".format(self.script_dir)

        self.effective_name = "effective"
        self.relax_name = "relax"

        self.MDP = "{}.mdp".format(self.effective_name)
        self.MDP_relax = "{}.mdp".format(self.relax_name)
        self.GRO = "{}.gro".format(self.effective_name)
        self.GRO_relax = "{}.gro".format(self.relax_name)

        self.sysname = None

        self.TOP = None
        self.NDX = None

        self.GROin = None
        self.TOPin = None
        self.MDPin = fmdp
        self.DATin = ""
        self.NDXin = ""

        self.LOG_fname = "LOG_CpHMD"
        self.LOG = open(self.LOG_fname, "a")

        self.GroDIR = None

        self.nCPUs = 1
        self.rcon = 0

        self.pH = None
        self.ionicstr = None
        self.nlit = 100
        self.nonit = 5

        self.reduced_titration = True
        self.rt_cycles = 10
        self.rt_limit = 0.001

        self.pkai = False

        self.titrating_group = "Protein"

        self.pypka_input = "pypka_input.gro"

        self.sites = None
        self.ff_tautomers = {}

        self.pypka_ffs_dir = "default"
        self.pypka_ffID = None
        self.pypka_nlit = "default"
        self.pypka_nonit = "default"

        self.input_type = {
            "TimeStep": float,
            "InitTime": float,
            "temp": float,
            "EffectiveSteps": int,
            "RelaxSteps": int,
            "InitCycle": int,
            "EndCycle": int,
            "nCycles": int,
            "GROin": str,
            "TOPin": str,
            "DATin": str,
            "sysname": str,
            "titrating_group": str,
            "nCPUs": int,
            "pH": float,
            "ionicstr": float,
            "sites": list,
            "GroDIR": str,
            "NDXin": str,
            "rcon": float,
            "rt_cycles": int,
            "rt_limit": float,
            "reduced_titration": bool,
            "pkai": bool,
            "pypka_ffs_dir": str,
            "pypka_ffID": str,
            "pypka_nlit": int,
            "pypka_nonit": int,
            "ffID": str,
        }

        # TODO implement '>=0' condition
        #                newEndCyle condition for InitCycle nCycles
        #                newEffectiveTime condition for EffectiveSteps TimeStep
        self.input_special_conditions = {}
        #    'TimeStep': '>=0',
        #    'InitTime': '>=0',
        #    'temp':     '>=0',
        #    'EffectiveSteps': '>=0'
        # }

    def get_simtime(self, cycle: int):
        simtime_begin = self.InitTime + self.EffectiveTime * cycle
        simtime_end = simtime_begin + self.EffectiveTime
        return simtime_begin, simtime_end

    def calcEndCyle(self):
        return self.InitCycle + self.nCycles

    def calcEffectiveTime(self):
        return self.EffectiveSteps * self.TimeStep

    def setFF(self):
        if self.ffID not in self.ff_options:
            raise IOError(
                "ffID parameter not valid. \nPlease choose one of available options: {}".format(
                    self.ffID, " ".join(self.ff_options.keys())
                )
            )

        self.ffDIR = "{0}/{1}.ff".format(self.script_dir, self.ffID)
        self.ff_dict = "{0}/protstates.dic".format(self.ffDIR)

        if not self.pypka_ffID:
            self.pypka_ffID = self.ff_options[self.ffID]

        for line in read_ff_dict_section(self.ff_dict, "tautomers"):
            cols = line.split()
            if cols[0] not in self.ff_tautomers:
                self.ff_tautomers[cols[0]] = []
            self.ff_tautomers[cols[0]] += [(cols[1:])]

    def setFileNames(self):
        self.TOP = "{}.top".format(self.sysname)
        self.NDX = "{}.ndx".format(self.sysname)
        self.DAT = "{}.dat".format(self.sysname)

    def setPypKaParams(self):
        self.pypka_params = {
            "structure": self.pypka_input,
            "pH": str(self.pH),
            "epsin": 2,
            "ionicstr": self.ionicstr,
            "pbc_dimensions": 0,
            "ncpus": self.nCPUs,
            "convergence": 0.01,
            "clean_pdb": False,
            "CpHMD_mode": True,
            "sts": "sts_cphmd",
            "nlit": self.nlit,
            "nonit": self.nonit,
            "ffID": self.pypka_ffID,
            "save_pdb": "delphi.pdb",
        }
        if self.pypka_ffs_dir != "default":
            self.pypka_params["ffs_dir"] = self.pypka_ffs_dir
        if self.pypka_nlit != "default":
            self.pypka_params["nlit"] = self.pypka_nlit
        if self.pypka_nonit != "default":
            self.pypka_params["nonit"] = self.pypka_nonit

        if self.pkai and self.reduced_titration:
            self.reduced_titration = False
            logging.warn("Reduced titration has been turned off")

    def check_mandatory_files(self):
        for f in (self.GROin, self.TOPin):
            if not os.path.isfile(f):
                raise IOError("Missing input file {}".format(f))
        for f in (self.NDXin, self.DATin):
            os.path.isfile(f)

    def read_input_mdp(self):
        def add_sites(mdp_sites):
            sites = []
            if mdp_sites == "all":
                return ["all"]

            for site in mdp_sites.split():
                if "N" in site or "C" in site:
                    sites.append(site)
                else:
                    try:
                        sites.append(int(site))
                    except:
                        raise IOError(
                            "sites has not been correctly defined. Example:\n"
                            "; sites = 1N 4 168 243C"
                        )
            return sites

        with open(self.MDPin) as f:
            for line in f:
                if line.startswith(";"):
                    line = line[1:]
                if line.startswith("SLURM"):
                    continue

                cleaned_line = remove_comments(line)

                parts = cleaned_line.split("=")
                param = parts[0].strip()

                if param in self:
                    param_value = parts[1].strip()
                    if param == "ref_t":
                        param_value = param_value.split()[0]

                    if param == "sites":
                        param_value = add_sites(param_value)

                    self[param] = param_value

        self.setFF()
        self.setFileNames()
        self.EndCycle = self.calcEndCyle()
        self.EffectiveTime = self.calcEffectiveTime()
        self.setPypKaParams()

        for i in self.__dict__:
            if self[i] == None:
                raise IOError("{} has not been defined.".format(i))

        self.check_mandatory_files()
