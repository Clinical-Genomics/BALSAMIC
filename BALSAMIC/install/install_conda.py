#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json


def conda_env_check(env_prefix):
  """
  Check if conda env exists.
  Input: Conda env prefix extracted from conda yaml file
  Output: True/False
  """
  try:
    p = json.loads(subprocess.check_output(["conda", "env", "list", "--json"], stderr=subprocess.STDOUT)) 
  except subprocess.CalledProcessError as e:
    print(e.output.decode())
    if verbose:
      raise e.output.decode()
   
  return env_prefix in p["envs"]

def conda_default_prefix():
  try:
    p = json.loads(subprocess.check_output(["conda", "info", "--json"], stderr=subprocess.STDOUT))
  except subprocess.CalledProcessError as e:
    print(e.output.decode())
    if verbose:
      raise e.output.decode()

  return p["conda_prefix"]
 
def get_prefix(yaml_handle):
  """
  Return conda env prefix from yaml file. Raise error if it doesn't exist.
  Input: yaml file handle
  Output: conda env prefix string
  """
  try:
    env_prefix = yaml.load(yaml_handle)
  except yaml.YAMLError as e:
    print("Error while reading yaml file")

  if "prefix" in env_prefix:
    env_prefix = env_prefix["prefix"] 
  else:
    env_prefix = False

  return env_prefix

def conda_remove(env_prefix):
  """
  Remove conda environment 
  """
   
  shellcmd = ["conda", "env", "remove", "--prefix", env_prefix, "--yes"]

  try:
    subprocess.check_output(shellcmd, stderr=subprocess.STDOUT)
  except subprocess.CalledProcessError as e:
    print(e.output.decode())
#    if verbose:
#      raise e.output.decode()

def conda_install(conda_yaml, env_prefix):
  """
  Install conda environment 
  """
  
  shellcmd = ["conda", "env", "create", "--file", conda_yaml, "--prefix", env_prefix]

  try:
    subprocess.check_output(shellcmd, stderr=subprocess.STDOUT)
  except subprocess.CalledProcessError as e:
    print(e.output.decode())
#    if verbose:
#      raise e.output.decode()

@click.command("install_conda", short_help="Installs required conda environments")
@click.option('-i','--input-conda-yaml',
              required=True,
              multiple=True,
              type=click.Path(),
              help='Input conda yaml file.')
@click.option('-s','--env-name-suffix',
              required=True,
              type=str,
              help='Mandatory alphanumeric suffix for environment name.')
@click.option('-o','--overwrite-env',
              default=False,
              is_flag=True,
              help='Overwite conda enviroment if it exists. Default = no. WARNING: The environment with matching name will be deleted')
@click.option('-d','--env-dir-prefix',
              type=str,
              help='Conda enviroment directory. It will be ignored if its provided within yaml file. Format: /path/env/envname.')
@click.pass_context

def install_conda(context, input_conda_yaml, env_dir_prefix, overwrite_env, env_name_suffix):
  """
  Installs conda environments from a conda yaml file.

  By default it doesn't overwrite if the environment by the same name exists. If _overwrite_ flag is provided, it tries to remove the
  enviroment first, and then install it in the path provided.
  """
  for fname in input_conda_yaml:
    
    if os.path.exists(fname):

      with open(fname) as yaml_handle:
        env_prefix = get_prefix(yaml_handle)

      env_name = os.path.basename(os.path.splitext(fname)[0]) + env_name_suffix 
      click.echo(click.style("Setting env name as %s" % env_name, fg='yellow' ))
      
      if not env_dir_prefix and not env_prefix:
        env_prefix = os.path.join(*[conda_default_prefix(), "envs", env_name])
        click.echo(click.style("Prefix is missing and no prefix was found in file. Setting prefix to %s" % env_prefix, fg='yellow' ))

      if conda_env_check(env_prefix) and not overwrite_env :
        click.echo(click.style("Conda env %s exists. Try using -o to overwrite." % env_prefix, fg='yellow' ))
        click.echo(click.style("Nothing to do.", fg='yellow'))

      elif not conda_env_check(env_prefix) :
        click.echo(click.style("Installing conda environment from %s" % fname, fg='green' ))
        conda_install(fname, env_prefix)
        
        click.echo("Conda environment %s was installed." % env_prefix)

      elif conda_env_check(env_prefix) and overwrite_env :
        click.echo(click.style("Removing old conda environment: %s" % env_prefix, fg='yellow'))
        conda_remove(env_prefix)
        click.echo(click.style("Installing conda environment from %s" % fname, fg='green' ))
        conda_install(fname, env_prefix)
        click.echo(click.style("Conda environment %s was installed." % env_prefix, fg='green' ))
