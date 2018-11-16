#script test ccmpred

#! /usr/bin/env python3

# Third-party modules
import os
from multiprocessing import Pool, cpu_count
from functools import partial
from datetime import datetime
from tqdm import tqdm
from docopt import docopt
from schema import Schema, And, Use, SchemaError

# Local modules
import src.parsing as parsing
from src.alignment import process
from src.score import Score
