import logging

from cli import CLI

logging.basicConfig(format="%(message)s", level=logging.INFO)


def main():
    CLI()


if __name__ == "__main__":
    main()
