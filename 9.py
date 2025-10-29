
import pandas as pd
import datamol as dm
from tqdm.auto import tqdm
from rdkit.Chem import rdDepictor
import mols2grid

import os , subprocess

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

df['ro5_ok'] = (df.mw <= 500) & (df.logp <=5) & (df.hbd <= 5) & (df.hba <= 10)





# df.to_csv("caluDF.csv", index=False)
# print("✅ Saved as 'caluDF.csv'")




# cluster the molecules
cluster_list = dm.cluster_mols(df.mol)
# create a list to hold the cluster ids, which we will add to the dataframe
cluster_idx = [-1] * len(df)
for i,cluster in enumerate(tqdm(cluster_list[0])):
    # align the structures for each cluster using Bemis Murcko frameworks
    dm.align.auto_align_many([df.mol.values[x] for x in cluster],copy=False,partition_method='scaffold')
    # add the cluster id to cluster_idx
    for c in cluster:
        cluster_idx[c] = i
# add a column with cluster ids to the dataframe
df['cluster'] = cluster_idx


cluster_sample_df = df.sort_values("cluster").drop_duplicates("cluster").copy()
# print(cluster_sample_df.head())




cluster_count_df = df.cluster.value_counts().to_frame().reset_index()
cluster_count_df.columns = ["cluster","count"]
cluster_sample_df = cluster_sample_df.merge(cluster_count_df,on="cluster")
print(cluster_sample_df.head())



# df.to_csv("clusDF.csv", index=False)
# print("✅ Saved as 'clusDF.csv'")




grid1 = mols2grid.MolGrid(cluster_sample_df,mol_col="mol", subset=["img","cluster","count"])
grid1.save("clus_grid.html")

print("Grid saved as 'clus_grid.html'")
# Open in default browser
subprocess.run(["xdg-open", "clus_grid.html"])