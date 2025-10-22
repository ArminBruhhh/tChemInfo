from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import mols2grid
from itertools import product



rxn_smarts = "[#6:6]-[#7H2:5].[#6:4]-[#6:1](-[#8H1])=[O:2]>>[#6:6]-[#7:5]-[#6:1](-[#6:4])=[O:2]"