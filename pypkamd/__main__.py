import logging
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))

from pypkamd.cli import CLI

logging.basicConfig(format="%(message)s", level=logging.INFO)


def main():
    CLI()


if __name__ == "__main__":
    main()
