import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib_inline



df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/practical_cheminformatics_tutorials/main/data/ChEMBL_hERG.csv")
# print(df.shape)
# print(df.head())

df.molregno = df.molregno.apply(str)

# print(df.dtypes)

df = df.query("assay_type == 'B'")
# print(df.shape)



gb = df.groupby("molregno")
# print(gb)

# print(gb.size().value_counts())


row_list = []
for k,v in gb:
    row_list.append([k,v.standard_value.mean()])
row_df = pd.DataFrame(row_list,columns=["name","standard_value"])

# print(row_df.shape)



## make it look better
sns.set(rc={'figure.figsize': (15, 12)})
sns.set(font_scale=1.5)
sns.set_style('whitegrid')

# ax = sns.displot(row_df.standard_value,kind="kde")
# ax.set(xlabel='pIC50', ylabel='Density')
# # print(ax)

# plt.show()


row_df["pIC50"] = -np.log10(row_df.standard_value * 1e-9)
# ax = sns.displot(row_df.pIC50,kind="kde")
# plt.show()


# print(row_df.dropna().shape)
# print(row_df.head())


# print(row_df.sort_values("pIC50",ascending=True).head())


score_9 = df.query("confidence_score == 9").copy()
print(score_9.head)