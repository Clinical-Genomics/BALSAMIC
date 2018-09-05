#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json
from yapf.yapflib.yapf_api import FormatFile

@click.command(
    "report",
    short_help="Create a sample config file from input sample data")
@click.pass_context

def report(context):
    """
    Prepares a config file for balsamic run_analysis. For now it is just treating json as dictionary and merging them as
it is. So this is just a placeholder for future.

    """
    pass
