__all__ = ['model_seed_database']

# Add Directories to PYTHONPATH
import sys, os
dir = os.path.dirname(__file__)

for subpackage in __all__:
   sys.path.append(os.path.join(dir,subpackage))

