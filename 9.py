
import pandas as pd
import datamol as dm
from tqdm.auto import tqdm
from rdkit.Chem import rdDepictor
import mols2grid



rdDepictor.SetPreferCoordGen(True)
pd.options.display.float_format = '{:,.2f}'.format

url = "https://raw.githubusercontent.com/PatWalters/yamc/main/data/HERG.smi"
df = pd.read_csv(url,sep=" ",names=["SMILES","Name","pIC50"])

# print(df.head())


df['IC50'] = dm.molar.log_to_molar(df.pIC50,'uM')



df['mol'] = dm.from_df(df, smiles_column="SMILES")


def max_ring_size(mol):
    """Get the size of the largest ring in a molecule

    :param mol: input_molecule
    :return: size of the largest ring or 0 for an acyclic molecule
    """
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    if len(atom_rings) == 0:
        return 0
    else:
        return max([len(x) for x in ri.AtomRings()])
    
    
my_prop_dict = {
    "mw" : dm.descriptors.mw,
    "logp" : dm.descriptors.clogp,
    "hbd" : dm.descriptors.n_lipinski_hbd,
    "hba" : dm.descriptors.n_lipinski_hba,
    "max_ring_size" : max_ring_size
}



prop_df = dm.descriptors.batch_compute_many_descriptors(df.mol,properties_fn=my_prop_dict,add_properties=False, progress=True)
# print(prop_df)

df = pd.concat([df,prop_df],axis=1)

print(df)



df_ro5_ok = df.query("mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10")
# print(len(df_ro5_ok))


# still some left 