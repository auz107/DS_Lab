__all__ = ["core","fba","gap_filling","COMETSpy","GAMETES","ancillary","customError","customWarning","pyomoSolverCreator","importData","strain_design", "globalVariables"]

# Add Directories to PYTHONPATH
import sys, os
dir = os.path.dirname(__file__)

for subpackage in __all__:
   sys.path.append(os.path.join(dir,subpackage))


