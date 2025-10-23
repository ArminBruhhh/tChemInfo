from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import mols2grid
from itertools import product
from rdkit.Chem import Draw
import subprocess , os

rxn_smarts = "[#6:1]-[#7H2:2].[#6:4]-[#6:3](-[#8H1])=O.[#7H2:5]-[c:6]:[c:7]-[#6:8](-[#8H1])=[O:9]>>[#6:4]-[c:3]1[n:5][c:6][c:7][c:8](=[O:9])[n:2]1-[#6:1]"

rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)
# img = Draw.ReactionToImage(rxn_mol, subImgSize=(800, 400))

# img.save("pics/22.png")
# subprocess.run(["xdg-open", os.path.abspath("pics/22.png")])



reagent_df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/practical_cheminformatics_tutorials/main/reaction/data/quinazoline_reagents.csv")





reagent_df['mol'] = reagent_df.SMILES.apply(Chem.MolFromSmiles)
# print(reagent_df)


acid_list = reagent_df.query("Type == 'carboxylic_acid'")[["mol","Name"]].values
aminobenzoic_list = reagent_df.query("Type == 'aminobenzoic_acid'")[["mol","Name"]].values
amine_list = reagent_df.query("Type == 'primary_amine'")[["mol","Name"]].values



# get the first reagent from each of the 3 lists
test_list = list([x[0][0] for x in [amine_list, acid_list, aminobenzoic_list]])

# run the reaction
product_mol = rxn_mol.RunReactants(test_list)[0][0]

# the reaction product needs to be cleaned up using SanitizeMol
Chem.SanitizeMol(product_mol)

print(product_mol)


# img  = Draw.MolToImage(product_mol, size=(300, 300))
# img.save( "pics/23.png")
# subprocess.run(["xdg-open", os.path.abspath("pics/23.png")])


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


prod_list = enumerate_library(rxn_mol, [amine_list, acid_list, aminobenzoic_list])
prod_df = pd.DataFrame(prod_list, columns=["SMILES","Name"])

print(prod_df)
