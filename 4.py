from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem import Draw #RDKit drawing
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import io


import os , subprocess



# for JNotebook
# IPythonConsole.ipython_useSVG = True
# IPythonConsole.molSize = 300, 300
# rdDepictor.SetPreferCoordGen(True)

query_mol = Chem.MolFromSmarts("c1ccccc1C")


ethyl_benzene = Chem.MolFromSmiles("c1ccccc1CC")


match = ethyl_benzene.GetSubstructMatch(query_mol)

# print(ethyl_benzene.HasSubstructMatch(query_mol))
# print(ethyl_benzene.GetSubstructMatch(query_mol))



# img  = Draw.MolToImage(ethyl_benzene, size=(300, 300), highlightAtoms=match)
# img.save( "pics/16-H.png")
# subprocess.run(["xdg-open", os.path.abspath("pics/16-H.png")])



# An list of SMILES with acetic acid and acetate anion
smiles_list = ["C(=O)O","C(=O)[O-]"]
# Convert SMILES to RDKit molecules
mol_list = [Chem.MolFromSmiles(x) for x in smiles_list]
# Create a query molecule from SMARTS
anion_query = Chem.MolFromSmarts("[-1]")
# Match the query to the molecules, convert the output to string so that MolsToGrid can display as a legend
match_list = [str(x.HasSubstructMatch(anion_query)) for x in mol_list]
# Draw the structures

img = Draw.MolsToGridImage(mol_list, legends=match_list)


## depend on your rdkit version you can use this or the one below to get the type of image
# print("Returned type:", type(img)) 


# extract PNG bytes from IPython image
img_bytes = img.data  # this is raw PNG data
pil_img = Image.open(io.BytesIO(img_bytes))

# # show and save
# pil_img.show()
# pil_img.save("pics/19.png")