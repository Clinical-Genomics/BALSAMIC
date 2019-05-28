#!/usr/bin/env python
import os
import subprocess
import hashlib
import yaml
import click
import logging
import sys
import json


def get_packages(yaml_handle):
    """
    Retrieve dependencies from conda yaml installationfile
    Input: yaml file name
    Output: list of installed packages
    """
    try:
        yaml_in = yaml.load(yaml_handle)
    except yaml.YAMLError as e:
        print("Error while reading yaml file")

    installed_packages = [s.split("=")[0] for s in yaml_in["dependencies"]]

    return installed_packages


def conda_env_check(env_prefix):
    """
    Check if conda env exists.
    Input: Conda env prefix extracted from conda yaml file
    Output: True/False
    """
    try:
        p = json.loads(
            subprocess.check_output(["conda", "env", "list", "--json"],
                                    stderr=subprocess.STDOUT))
    except subprocess.CalledProcessError as e:
        print(e.output.decode())
        if verbose:
            raise e.output.decode()

    return env_prefix in p["envs"]


def conda_default_prefix():
    try:
        p = json.loads(
            subprocess.check_output(["conda", "info", "--json"],
                                    stderr=subprocess.STDOUT))
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

    shellcmd = [
        "conda", "env", "create", "--file", conda_yaml,
        "--prefix", env_prefix
    ]

    try:
        subprocess.check_output(shellcmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(e.output.decode())


#    if verbose:
#      raise e.output.decode()


@click.command("install", short_help="Installs required conda environments")
@click.option('-i',
              '--input-conda-yaml',
              required=True,
              multiple=True,
              type=click.Path(),
              help='Input conda yaml file.')
@click.option('-s',
              '--env-name-suffix',
              required=True,
              help='Mandatory alphanumeric suffix for environment name.')
@click.option(
    '-o',
    '--overwrite-env',
    show_default=True,
    is_flag=True,
    help=
    'Overwite conda enviroment if it exists. Default = no. WARNING: The environment with matching name will be deleted'
)
@click.option(
    '-d',
    '--env-dir-prefix',
    help=
    'Conda enviroment directory. It will be ignored if its provided within yaml file. Format: /path/env/envname.'
)
@click.option(
    '-p',
    '--packages-output-yaml',
    required=True,
    type=click.Path(),
    help=
    'Output a yaml file containing packages installed in each input yaml file.'
)
@click.option(
    '-t',
    '--env-type',
    type=click.Choice(['D', 'P', 'S']),
    default='D',
    help=
    'Environment type. P: Production, D: Development, S: Stage. It will be added to filename: "[D|P|S]_"+filename+env-name-suffix'
)
@click.pass_context
def install(context, input_conda_yaml, env_dir_prefix, overwrite_env,
            env_name_suffix, packages_output_yaml, env_type):
    """
    Installs conda environments from a conda yaml file.
    
    By default it doesn't overwrite if the environment by the same name exists. If _overwrite_ flag is provided, it tries to remove the
    enviroment first, and then install it in the path provided.
    """

    conda_packages = dict()

    for fname in input_conda_yaml:

        if os.path.exists(fname):

            with open(fname) as yaml_handle:
                env_prefix = get_prefix(yaml_handle)

            env_base_name = os.path.basename(os.path.splitext(fname)[0])

            env_name = env_type + "_" + env_base_name + env_name_suffix

            click.echo(
                click.style("Setting env name as %s" % env_name, fg='yellow'))

            if not env_dir_prefix and not env_prefix:
                env_prefix = os.path.join(
                    *[conda_default_prefix(), "envs", env_name])
                click.echo(
                    click.style(
                        "Prefix is missing and no prefix was found in file. Setting prefix to %s"
                        % env_prefix,
                        fg='yellow'))
            elif env_dir_prefix and not env_prefix:
                env_prefix = os.path.join(*[env_dir_prefix, env_name])
                click.echo(
                    click.style(
                        "A conda prefix is provided. Setting prefix to %s" %
                        env_prefix,
                        fg='yellow'))


            if conda_env_check(env_prefix) and not overwrite_env:
                click.echo(
                    click.style(
                        "Conda env %s exists. Try using -o to overwrite." %
                        env_prefix,
                        fg='yellow'))
                click.echo(click.style("Nothing to do.", fg='yellow'))

            elif not conda_env_check(env_prefix):
                click.echo(
                    click.style("Installing conda environment from %s" % fname,
                                fg='green'))

                conda_install(fname, env_prefix)

                with open(fname) as yaml_handle:
                    conda_packages[env_prefix] = get_packages(yaml_handle)

                click.echo("Conda environment %s was installed." % env_prefix)

            elif conda_env_check(env_prefix) and overwrite_env:
                click.echo(
                    click.style("Removing old conda environment: %s" %
                                env_prefix,
                                fg='yellow'))
                conda_remove(env_prefix)
                click.echo(
                    click.style("Installing conda environment from %s" % fname,
                                fg='green'))
                conda_install(fname, env_prefix)

                with open(fname) as yaml_handle:
                    conda_packages[env_prefix] = get_packages(yaml_handle)

                click.echo(
                    click.style("Conda environment %s was installed." %
                                env_prefix,
                                fg='green'))

    yaml.dump(conda_packages, open(packages_output_yaml, 'w'))
