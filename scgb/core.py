#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################################################################
# -*- coding:utf8 -*-
# version: 1.0
# script: core.py
# author: sundayzz
# mail: sundayzz@sina.cn
# main modules 
############################################################################


"""Core CLI entrypoint for the simplified SCGB toolkit.

The public workflow is intentionally reduced to three modules:

1. predict
   Predict assembly contigs as plasmid or chromosome contigs.
2. binning
   Bin chromosome contigs with COMEBin, SemiBin2 and MetaDecoder.
3. link
   Use single-cell reads to link plasmid contigs to binned genomes.

Usage
-----
# List available subcommands
python core.py --help

# Run a specific module (example: binning)
python core.py binning [MODULE-SPECIFIC-ARGS...]

# Enable file logging and set console log level
python core.py --logfile run.log --loglevel DEBUG binning [ARGS...]

Subcommands
-----------
- predict
- binning
- link

Notes
-----
Each subcommand is expected to expose a `main` Click command object so it can
be registered here with `add_command(...)`.
"""

from __future__ import annotations
import sys
from pathlib import Path

import click
from loguru import logger

# Import the three public workflow modules.
from modules import Binning, Link, PredictPlasmid


@click.group(
    help=(
        "SCGB simplified CLI. Choose one module: "
        "predict, binning, link"
    ),
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--logfile",
    type=click.Path(dir_okay=False, writable=True),
    default=None,
    help="Optional path to also log output to a file (appends).",
)
@click.option(
    "--loglevel",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    default="INFO",
    show_default=True,
    help="Console log verbosity.",
)
def main(logfile: str | None, loglevel: str) -> None:
    """
    Initialize global logging and act as a container for subcommands.

    Parameters
    ----------
    logfile : str | None
        If provided, logs are additionally written to this file (append mode).
    loglevel : str
        Log level for console output (DEBUG/INFO/WARNING/ERROR).
    """
    # Reset Loguru’s default handler and configure console logging.
    logger.remove()
    logger.add(
        sys.stdout,
        level=loglevel,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}",
    )

    # Optional file logging (append mode).
    if logfile:
        logger.add(
            Path(logfile).expanduser(),
            level=loglevel,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
            mode="a",
            enqueue=True,  # safer in multi-process contexts
        )


# Register the simplified public subcommands.
# Each module exposes a Click command object named `main`.
main.add_command(PredictPlasmid.main, name="predict")
main.add_command(Binning.main,        name="binning")
main.add_command(Link.main,           name="link")


if __name__ == "__main__":
    main()
