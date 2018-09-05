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

from .sample import sample as sample_command
from .report import report as report_command
 
@click.group()
@click.pass_context

def config(context):
    "create config files required for running the pipeline and reporting it"
    pass

config.add_command(sample_command)
config.add_command(report_command)
