import pandas as pd
import numpy as np
import json
import sys
from GPfates import GPfates

temp_folder = sys.argv[1]

# Load params
p = json.load(open(temp_folder + "/params.json", "r"))

# load data
etpm = pd.read_table(temp_folder + "/expression.tsv", index_col=0)
etpm = etpm[(etpm > p["log_expression_cutoff"]).sum(1) > p["min_cells_expression_cutoff"]]
logexp = np.log10(etpm + 1)

cellinfo = pd.read_table(temp_folder + "/cellinfo.tsv", index_col=0)

m = GPfates.GPfates(cellinfo, logexp)

print("Dimensionality reduction--------------------------------------")
m.dimensionality_reduction()

print("Story DR------------------------------------------------------")
m.store_dr(dims=range(p["ndims"])) # store the dr in the sample table (m.s), so it can be used in the gplvm

print("Infer pseudotime----------------------------------------------")
m.infer_pseudotime(s_columns=["bgplvm_" + str(i) for i in range(ndims)]) # use the first two components to infer pseudotime

print("Model cell fates----------------------------------------------")
m.model_fates(C=p["nfates"])

print("Saving--------------------------------------------------------")
m.s.pseudotime.to_csv(temp_folder + "/pseudotimes.csv")
pd.DataFrame(m.fate_model.phi, index=m.s.pseudotime.index).to_csv(temp_folder + "/phi.csv")

dr = m.dr_models["bgplvm"]
dr2 = dr.X.mean[:, :]
pd.DataFrame(dr2.tolist(), index=m.s.pseudotime.index).to_csv(temp_folder + "/dr.csv")

print("Finished------------------------------------------------------")
