from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdDepictor
import pandas as pd
from rdkit.Chem import EnumerateStereoisomers
from rdkit.Chem.MolStandardize import rdMolStandardize
import subprocess , os





import mols2grid


rdDepictor.SetPreferCoordGen(True)

smi_1 = "F[C@](Cl)(Br)c1ccc(cc1)C(F)(Cl)Br"
mol_1 = Chem.MolFromSmiles(smi_1)
smi_2 = "F[C@@](Cl)(Br)c1ccc(cc1)C(F)(Cl)Br"
mol_2 = Chem.MolFromSmiles(smi_2)



mol_list = [mol_1, mol_2]



data = {
    'Molecule': [mol_1, mol_2],
    'Name': ['Stereoisomer 1', 'Stereoisomer 2'],
    'SMILES': [smi_1, smi_2]
}
df = pd.DataFrame(data)

# grid = mols2grid.MolGrid(df, subset=["img", "Name", "SMILES"])
# grid.save("two_grid.html")

# print("Grid saved as 'two_grid.html'")
# # Open in default browser
# subprocess.run(["xdg-open", "two_grid.html"])
# # still some left



mol_1_isomer_list = [x for x in EnumerateStereoisomers.EnumerateStereoisomers(mol_1)]
print(mol_1_isomer_list)




img = MolsToGridImage(mol_1_isomer_list)
img.save("pics/24-stereoisomers.png")