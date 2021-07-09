import argparse
import logging
import time

from configs import Config
from cphmd import create_cphdm_directory, final_cleanup, run_cphmd
from misc import get_curtime
from topology import Topology


def CLI():
    init_time = time.time()
    parser = argparse.ArgumentParser(description="PypKa-MD")
    parser.add_argument("mdp", type=str, help="MDP file")

    args = parser.parse_args()

    logging.info("### Starting PypKa-MD ###")

    Config.storeParams(args.mdp)

    info = "Simulation run by {} @ {}\n" "Initial time: {}".format(
        Config.md_configs.USER, Config.md_configs.HOST, get_curtime()
    )
    logging.info(info)

    top = Topology(Config.md_configs)

    create_cphdm_directory()
    top.read_index(Config.md_configs.NDX)

    run_cphmd(top)

    final_cleanup(
        Config.md_configs.sysname,
        Config.md_configs.tmpDIR,
        Config.md_configs.effective_name,
    )

    end_time = time.time()
    duration = end_time - init_time
    duration_formated = time.strftime("%H h %M m %S s", time.gmtime(duration))
    logging.info("\nCpHMD block duration: {}".format(duration_formated))
    logging.info("### PypKa-MD finished succesfully ###")
