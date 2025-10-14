from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem import Draw #RDKit drawing

# for JNotebook
# IPythonConsole.ipython_useSVG = True
# IPythonConsole.molSize = 300, 300
# rdDepictor.SetPreferCoordGen(True)

query_mol = Chem.MolFromSmarts("c1ccccc1C")


ethyl_benzene = Chem.MolFromSmiles("c1ccccc1CC")




print(ethyl_benzene.HasSubstructMatch(query_mol))