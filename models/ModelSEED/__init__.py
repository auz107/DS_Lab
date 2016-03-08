__all__ = ['ModelSeed_compounds_website','ModelSeed_reactions_website']

# Add Directories to PYTHONPATH
import sys, os
dir = os.path.dirname(__file__)

for subpackage in __all__:
   sys.path.append(os.path.join(dir,subpackage))

