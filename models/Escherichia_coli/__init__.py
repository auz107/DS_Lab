__all__ = ["iAF1260","iJO1366"]

# Add Directories to PYTHONPATH
import sys, os
dir = os.path.dirname(__file__)

for subpackage in __all__:
   sys.path.append(os.path.join(dir,subpackage))

