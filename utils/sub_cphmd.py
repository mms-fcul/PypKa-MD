#! /usr/bin/python3

import argparse
import os
import getpass

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""
Slurm Submission

Example:
python3 sub_cphmd.py <mdp file>

""",
)

# Mandatory Arguments
parser.add_argument("mdp", help=" settings file name", action="store")

args = parser.parse_args()

OPTIONAL_PARAMS = ["NDXin"]


class slurm_params:
    def __init__(self):
        self.segments = None
        self.partition = None
        self.nCPUs = None

        self.sysname = None

        self.GROin = None
        self.TOPin = None
        self.DATin = None
        self.NDXin = None

        self.nCycles = None
        self.dt = None
        self.nsteps = None
        self.effective_time = None

        self.user = getpass.getuser()
        self.cur_dir = os.getcwd()
        self.runfile = os.path.abspath(__file__)

        self.requeue = 1

        self.cphdm_params = (
            "sysname",
            "nCPUs",
            "GROin",
            "TOPin",
            "DATin",
            "NDXin",
            "nCycles",
        )
        self.md_params = ("dt", "nsteps")

    def set_eff_time(self):
        self.effective_time = float(self.dt) * int(self.nsteps) * int(self.nCycles)

    def __getitem__(self, name):
        if name not in self.__dict__:
            raise Exception("{} is not a parameter of slurm".format(name))
        return self.__dict__[name]

    def __setitem__(self, name, key):
        self.__dict__[name] = key

    def __contains__(self, item):
        return item in self.__dict__

    def check_params(self):
        for param, param_value in self.__dict__.items():
            if not param_value and param not in OPTIONAL_PARAMS:
                raise IOError(f"slurm parameter {param} has not been set.")


def get_param_from_line(line, delimiter=";SLURM"):
    line = line.split(delimiter)[1].split(";")[0]
    cols = line.split("=")
    param = cols[0].strip()
    param_value = None
    if len(cols) > 1:
        param_value = cols[1].strip().replace(" ", "")
    return param, param_value


sparams = slurm_params()

if not os.path.isfile(args.mdp):
    raise IOError(f"{args.mdp} does not exist.")

with open(args.mdp) as f:
    for line in f:
        if line.startswith(";SLURM"):
            param, param_value = get_param_from_line(line)
            sparams[param] = param_value

        if line.startswith(";") and "=" in line:
            param, param_value = get_param_from_line(line, delimiter=";")
            if param in sparams.cphdm_params:
                sparams[param] = param_value

        elif "=" in line:
            param, param_value = get_param_from_line(";" + line, delimiter=";")
            if param in sparams.md_params:
                sparams[param] = param_value

sparams.set_eff_time()
sparams.check_params()

if not os.path.isfile(f"{sparams.sysname}_000.gro"):
    os.system(f"ln -s {sparams.GROin} {sparams.sysname}_000.gro")
if not os.path.isfile(f"{sparams.sysname}_000.top"):
    os.system(f"ln -s {sparams.TOPin} {sparams.sysname}_000.top")

find_current_segment = """\
#!/bin/bash -e

# Finds current block
i=1
j=001

while [ -f {sysname}_${{j}}.occ ] ; do
    i=$((i+1))
    j=`printf "%03d\\n" ${{i}}`
done

if [ $i -gt {max_block} ]
then
    exit
fi

k=$((i-1))
l=`printf "%03d\\n" ${{k}}`
""".format(
    sysname=sparams.sysname, max_block=sparams.segments
)

create_blockinfo = """
# Info about the Job:
echo "Job executed in Machine $HOSTNAME" > {rundir}/{sysname}_${{j}}.blockinfo
echo "Job executed with {ncpus} processors" >> {rundir}/{sysname}_${{j}}.blockinfo
echo "Job executed in DIR: /tmp/{user}_CpHMD$$ " >> {rundir}/{sysname}_${{j}}.blockinfo

echo "" >> {rundir}/{sysname}_${{j}}.blockinfo
echo -e "Job started on: " >> {rundir}/{sysname}_${{j}}.blockinfo
date >> {rundir}/{sysname}_${{j}}.blockinfo
""".format(
    rundir=sparams.cur_dir,
    sysname=sparams.sysname,
    ncpus=sparams.nCPUs,
    user=sparams.user,
)

edit_mdp = """
# Edit mdp file
sed "s/;*.sysname.*=.*/; sysname = {sysname}_${{j}}/" {rundir}/{mdp} > {sysname}_${{j}}.mdp
sed -i "s/;*.GROin.*=.*/; GROin = {sysname}_${{l}}.gro/" {sysname}_${{j}}.mdp
sed -i "s/;*.TOPin.*=.*/; TOPin = {sysname}_${{l}}.top/" {sysname}_${{j}}.mdp
sed -i "s/tinit.*=.*/tinit = $(({effective_steps}*$k))/" {sysname}_${{j}}.mdp
""".format(
    rundir=sparams.cur_dir,
    mdp=args.mdp,
    sysname=sparams.sysname,
    effective_steps=int(sparams.effective_time),
)


copy_files = """
# Copy important files for local directory; first gro is called
mkdir -p /tmp/{user}_CpHMD$$

# line needed to mount /programs using autofs
ls /programs/CpH-MD/CpH-MD_v2.04/scripts/ >/dev/null 2>&1; sleep 5

cp {rundir}/{sysname}_${{j}}.mdp /tmp/{user}_CpHMD$$
cp {rundir}/{sysname}_${{l}}.{{top,gro}} /tmp/{user}_CpHMD$$
cp {rundir}/{dat} /tmp/{user}_CpHMD$$
cp {rundir}/{sysname}_${{j}}.blockinfo /tmp/{user}_CpHMD$$
""".format(
    user=sparams.user,
    sysname=sparams.sysname,
    dat=sparams.DATin,
    rundir=sparams.cur_dir,
)

if sparams.NDXin:
    copy_files += f"cp {sparams.cur_dir}/{sparams.NDXin} /tmp/{sparams.user}_CpHMD$$\n"


run_cphmd = """
# Enter local directory
cd /tmp/{user}_CpHMD$$

# Run Constant-pH MD segment:
nice -n 19 python3.8 /home/pedror/MMS@FCUL/cphmd/pypkamd {sysname}_${{j}}.mdp
""".format(
    user=sparams.user, sysname=sparams.sysname
)

copy_output = """
# Copy files and return to shared directory
cp {sysname}_${{j}}* {rundir}
if (for f in {sysname}_${{j}}*; do diff $f {rundir}/$f; done); then
cd {rundir}
rm -rf /tmp/{user}_CpHMD$$
gzip -9 {sysname}_${{j}}.log
else
echo "Error in file copy... please check local files" >> {rundir}/{sysname}_${{j}}.blockinfo
exit 1
fi
""".format(
    user=sparams.user, sysname=sparams.sysname, rundir=sparams.cur_dir
)

finish_blockinfo = """
echo "" >> {rundir}/{sysname}_${{j}}.blockinfo
echo -e "Job finished on: " >> {rundir}/{sysname}_${{j}}.blockinfo
date >> {rundir}/{sysname}_${{j}}.blockinfo
""".format(
    sysname=sparams.sysname, rundir=sparams.cur_dir
)

relaunch = """
# Launch next job before exiting:
if [ ${{i}} -lt {segments} ] # usually 1 segment == 1 ns
then
    cd {rundir}; {runfile} {mdp}
fi
""".format(
    rundir=sparams.cur_dir,
    runfile=sparams.runfile,
    segments=sparams.segments,
    mdp=args.mdp,
)

sections = (
    find_current_segment,
    create_blockinfo,
    edit_mdp,
    copy_files,
    run_cphmd,
    copy_output,
    finish_blockinfo,
    relaunch,
)
with open(sparams.sysname + ".slurm", "w") as f_new:
    content = ""
    for section in sections:
        content += section
    f_new.write(content)

sbatch_cmd = "sbatch"
if sparams.requeue == 1:
    sbatch_cmd += " --requeue"
sbatch_cmd += f" -p {sparams.partition} -N 1 -n {sparams.nCPUs} -o {sparams.sysname}.sout -e {sparams.sysname}.serr {sparams.sysname}.slurm"

os.system(sbatch_cmd)

print(
    f"Job submitted to Partition(s): {sparams.partition} with {sparams.nCPUs} Processors"
)
