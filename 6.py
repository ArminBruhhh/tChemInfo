from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import mols2grid
from itertools import product
import subprocess, os
from rdkit.Chem import Draw #RDKit drawing

rxn_smarts = "[#6:6]-[#7H2:5].[#6:4]-[#6:1](-[#8H1])=[O:2]>>[#6:6]-[#7:5]-[#6:1](-[#6:4])=[O:2]"


rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)


img = Draw.ReactionToImage(rxn_mol, subImgSize=(800, 400))
# img.save("pics/21.png")
# subprocess.run(["xdg-open", os.path.abspath("pics/21.png")])



df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/practical_cheminformatics_tutorials/main/reaction/data/amide_reagents.csv")

df['mol'] = df.SMILES.apply(Chem.MolFromSmiles)

acid_df = df.query("Type == 'carboxylic_acid'").copy()
amine_df = df.query("Type == 'primary_amine'").copy()





# grid1 = mols2grid.MolGrid(acid_df, subset=["img", "Name"])
# grid1.save("acid_grid.html")

# print("Grid saved as 'acid_grid.html'")
# # Open in default browser
# subprocess.run(["xdg-open", "acid_grid.html"])




# grid2 = mols2grid.MolGrid(amine_df, subset=["img", "Name"])
# grid2.save("amine_grid.html")

# print("Grid saved as 'amine_grid.html'")
# # Open in default browser
# subprocess.run(["xdg-open", "amine_grid.html"])


acid_list = acid_df[['mol','Name']].values
amine_list = amine_df[['mol','Name']].values

def enumerate_library(rxn_mol, reagent_lol):
    prod_list = []
    # itertools.product generates all combinations of reactants
    for reagents in product(*reagent_lol):
        mol_list = [x[0] for x in reagents]
        name_list = [str(x[1]) for x in reagents]
        name = "_".join(name_list)
        prod = rxn_mol.RunReactants(mol_list)
        if prod is not None and len(prod):
            product_mol = prod[0][0]
            Chem.SanitizeMol(product_mol)
            prod_list.append([Chem.MolToSmiles(product_mol), name])
    return prod_list


prod_list = enumerate_library(rxn_mol,[amine_list, acid_list])


prod_df = pd.DataFrame(prod_list,columns=["SMILES","Name"])




grid3 = mols2grid.MolGrid(prod_df, subset=["img", "Name"])
grid3.save("prod_grid.html")

print("Grid saved as 'prod_grid.html'")
# Open in default browser
subprocess.run(["xdg-open", "prod_grid.html"])
