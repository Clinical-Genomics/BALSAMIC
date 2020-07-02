import os
import re
import json
import pytest
from unittest import mock
import click


# Given standard run fixtures, execute command and verify dag graph exists
# Given invalid, existing fastq,  -> invoke cli > error
# Given an empty reference config -> invoke command > error
# Given analysis dir without write permissions -> invoke command > error


