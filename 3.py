from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem import Draw #RDKit drawing
from rdkit.Chem.Draw import MolsToGridImage
from py2opsin import py2opsin

IPythonConsole.ipython_useSVG = True
IPythonConsole.molSize = 300, 300
rdDepictor.SetPreferCoordGen(True)


# rest of the code
## its better to run this on J-notebook but in case use for python you can:
import os , subprocess


ethanol = Chem.MolFromSmiles("CCO")

Draw.MolToFile(ethanol, "2.png", size=(300, 300))

# subprocess.run(["xdg-open", os.path.abspath("2.png")])








### Exercise

# Write the SMILES and display a table of chemical structures for the following. If you don't know the structures, google the names. 

# 1. 2-menthylpentene
# 2. isopropanal
# 3. 2-pentyn-1-ol
# 4. 1,2,2,3-tetrafluorobutane
# 5. propanoic acid
# 6. 2-t-butyl-3-hydroxy-propane



# 1. 2-menthylpentene
MP12 = Chem.MolFromSmiles("C=C(C)CCC") 
Draw.MolToFile(MP12, "3.png", size=(300, 300))


# define a text buffer for our examples
examples = """C(C)(C)O isopropanol
C(Cl)(Cl)(Cl)Cl carbon tetrachloride
CC(=O)O acetic acid"""
# not that we use the second argumen to split to only return two tokens
smiles_list = [x.split(" ",1) for x in examples.split("\n")]
# print(smiles_list)


## IUPAC to smiles 
from py2opsin import py2opsin

smiles_string = py2opsin(
    chemical_name=["ethane", "isopropanal", "2-pentyn-1-ol", "1,2,2,3-tetrafluorobutane", "propanoic acid", "2-t-butyl-3-hydroxy-propane",
    "1,2-dimethylcyclopropane"
    ,"1-methyl-3,3-dimethylcyclohexane"
    ,"piperazine"
    ,"octahydro-1H-indene"
    ,"norbornane"
    ,"cyclopentanol"
],
    output_format="SMILES"
)
# print(smiles_string)


mol_list = [Chem.MolFromSmiles(x) for x in smiles_string]
print (mol_list) 

for i, mol in enumerate(mol_list, start=4):
    Draw.MolToFile(mol, f"{i}.png", size=(300, 300))

