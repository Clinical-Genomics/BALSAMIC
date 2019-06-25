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

from BALSAMIC.commands.config.sample import sample as sample_command
 
@click.group()
@click.pass_context

def config(context):
    "create config files required for running the pipeline."
    pass

config.add_command(sample_command)
