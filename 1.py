import os , subprocess
from rdkit import Chem #RDKit Chemistry
from rdkit.Chem.Draw import IPythonConsole #RDKit drawing
from rdkit.Chem import Draw #RDKit drawing
# A few settings to improve the quality of structures 
from rdkit.Chem import rdDepictor
from rdkit.Chem import PandasTools #Add the ability to add a molecule to a dataframe
import mols2grid #The mols2grid library provides a convenient way of displaying molecules in a grid
import requests

IPythonConsole.ipython_useSVG = True
rdDepictor.SetPreferCoordGen(True)
mol = Chem.MolFromSmiles("c1ccccc1") #from smiles data extractions

# print([mol])
# Draw.MolToFile(mol, "1.png", size=(300, 300))
# import subprocess, os
# subprocess.run(["xdg-open", os.path.abspath("pics/1.png")])



#data extractions from sdf format
url = "https://raw.githubusercontent.com/PatWalters/practical_cheminformatics_tutorials/main/data/example_compounds.sdf"
r = requests.get(url)
bytes_written = open('example_compounds.sdf', 'w').write(r.text)



#printing the data
mols = [ x for x in Chem.SDMolSupplier("example_compounds.sdf") ]

         # #drawing the grid
        # Draw.MolsToGridImage(mols,molsPerRow=4,useSVG=True)
        # mols2grid.display(mols)

mols2grid.display(mols)
