import os
import argparse
import logging
from pathlib import Path

from openmmtools import cache, utils
from openmmtools.multistate import MultiStateReporter

from perses.dispersed.utils import configure_platform
from perses.samplers.multistate import HybridRepexSampler

# Set up logger
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)

# Configure platform
platform = configure_platform(utils.get_fastest_platform().getName())

# Load arguments
parser = argparse.ArgumentParser(description='resume repex')
parser.add_argument('dir', type=str, help='path to input/output dir')
parser.add_argument('phase', type=str, help='phase of the simulation used in storage filename')
parser.add_argument('total_iterations', type=int, help='total number of iterations desired')
args = parser.parse_args()

# Load repex simulation
directory_number = Path(args.dir).parts[-2]
reporter_file = os.path.join(args.dir, f"{directory_number}_{args.phase}.nc")
reporter = MultiStateReporter(reporter_file, checkpoint_interval=100)
sampler = HybridRepexSampler.from_storage(reporter)
sampler.energy_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)
sampler.sampler_context_cache = cache.ContextCache(capacity=None, time_to_live=None, platform=platform)

# Determine how many more iterations are needed
total_iterations = args.total_iterations
iterations =  total_iterations - sampler.iteration

# Resume simulation
sampler.extend(iterations)
del sampler
