from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import Draw #RDKit drawing
import os , subprocess

m1 = Chem.MolFromSmiles("OC(=O)c1cccc(C=O)c1")


recursive_acid_smarts = "[$(C(=O)[OH])]"
recursive_acid_query = Chem.MolFromSmarts(recursive_acid_smarts)


match = m1.GetSubstructMatch(recursive_acid_query)


img  = Draw.MolToImage(m1, size=(300, 300), highlightAtoms=match)
img.save( "pics/20-H.png")
subprocess.run(["xdg-open", os.path.abspath("pics/20-H.png")])
