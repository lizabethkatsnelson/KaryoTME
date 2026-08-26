#%%
############### IMPORTING NEEDED LIBRARIES AND MODULES FOR THE SCRIPT ###############
from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statistics
import sys

sys.path.append("E:/Immune_Evasion_Project/Mario/Scripts_Python")
# %%
############### IMPORTING FUNCTIONS LIBRARY ###############

# We will import the 'Functions' library from our directory, make sure the file is present in the defined Directory above.
import Functions_Mario
from Functions_Mario import *
from Functions_Mario import merge_dataframes
from Functions_Mario import drop_columns
import Functions
from Functions import *

# This is for reloading the Library above. 
# import importlib
# importlib.reload(Functions_Mario)

# %%
############### IMPORTING REQUIRED DATA ###############

# You can change the value of cancer_type to the desired cancer type, and the code will automatically update the file paths accordingly. For example, if you want to analyze the data for lung cancer (LUAD), you can change cancer_type = "LUAD".
cancer_type = "SKCM"  # specify the cancer type here


immune_score = pd.read_csv(f"E:/Immune_Evasion_Project/Mario/Cancer_Types/{cancer_type}/Data/{cancer_type}_immune_score.txt", sep="\t", header=0)
SCNA_absolute = pd.read_csv(f"E:/Immune_Evasion_Project/Mario/Cancer_Types/{cancer_type}/Data/{cancer_type}_purity_&_ploidy.txt", sep="\t", header=0) # Gene x sample
#SCNA_absolute =  pd.read_csv(f"E:/Immune_Evasion_Project/Dania/Data/{cancer_type}_all_data_by_genes_purity_rescale.txt", sep = "\t", header =0) # This is actually Pathology data set, not absolute.
RNA = pd.read_csv(f"E:/Immune_Evasion_Project/Mario/Cancer_Types/{cancer_type}/Data/{cancer_type}_rnaseq_RSEM.txt", sep="\t", header=0) # 520x20k
all_genes = pd.read_csv("E:/Immune_Evasion_Project/Dania/original_candidate_lists/tabcat2_HNSC.SKCM.BLCA.ESCA.LUAD.LUSC.PAAD_filtereddel.csv", header=0) # gene x sample
# %%
############### DEFINING ALL FUNCTIONS USED IN THIS SCRIPT ###############

### CORRELATION FUNCTIONS ###
def corr_pvals_SCNA_RNA(df1, df2):
    ''' Function to get correlation between 2 dataframes - mainly SCNA v RNA, and SCNA v protein in both TCGA and CCLE
    df1 = SCNA with IS
    df2 = RNA with IS '''

    data_unfiltered = {'corr_val':None, 'p_val':None} 
    data_filtered = {'corr_val':None, 'p_val':None}
    # dataframes that will be output
    corr_unfiltered = pd.DataFrame(data_unfiltered, index = df1.columns[df1.columns != "Ranked Sum"])
    corr_filtered = pd.DataFrame(data_filtered, index = df1.columns[df1.columns != "Ranked Sum"])

    for i in df1.columns[df1.columns != "Ranked Sum"]: # don't take the last column (ranked sum col)
        # unfiltered
        corr_result_unfiltered = stats.spearmanr(df1[i], df2[i], axis=0)
        corr_unfiltered.corr_val[i] = corr_result_unfiltered[0]
        corr_unfiltered.p_val[i] = corr_result_unfiltered[1]

        # filtered
        #filter out samples w gains ie only take samples that have losses or neutral
        col_filtered = df1[i][(df1[i]<0.2)] 
        col_filtered_indx = col_filtered.index
        df2_filtered = df2[i][df2.index.isin(col_filtered_indx)] # filter df2 only in the indx list

        corr_result_filtered = stats.spearmanr(col_filtered, df2_filtered, axis=0) 
        corr_filtered.corr_val[i] = corr_result_filtered[0]
        corr_filtered.p_val[i] = corr_result_filtered[1]

    corr_filtered = corr_filtered.reset_index().rename(columns = {"index": "Gene-Symbol"})
    corr_unfiltered = corr_unfiltered.reset_index().rename(columns = {"index": "Gene-Symbol"})

    return corr_unfiltered, corr_filtered, col_filtered_indx, df2_filtered # call each df using index [0] or [1]

def corr_pvals_SCNA_IS(df1):
    '''Function to get correlation values and pvals (2 cols)
    Correlation involving IS will always have only 1 argument ie 1 df, since the IS is always in the last column
    - Filtering out samples that have gains for that gene
    - Output: 2 dfs with unfiltered and corr vals after filtered out gains'''
    data_unfiltered = {'corr_val':None, 'p_val':None} 
    data_filtered = {'corr_val':None, 'p_val':None}

    # dataframes that will be output
    corr_unfiltered = pd.DataFrame(data_unfiltered, index = df1.columns[df1.columns != "Ranked Sum"])
    corr_filtered = pd.DataFrame(data_filtered, index = df1.columns[df1.columns != "Ranked Sum"])

    # create dict to store sample names that were filtered for less than 0.2
    keyList = []
    for i in df1.columns[df1.columns != "Ranked Sum"]:
        keyList.append(i)
    loss_samples = {key: None for key in keyList}

    for i in df1.columns[df1.columns != "Ranked Sum"]:
        # unfiltered
        corr_result_unfiltered = stats.spearmanr(df1[i], df1["Ranked Sum"], axis=0)
        corr_unfiltered.corr_val[i] = corr_result_unfiltered[0]
        corr_unfiltered.p_val[i] = corr_result_unfiltered[1]


        # filtered
        #filter out samples w gains ie only take samples that have losses or neutral
        col_filtered = df1[i][(df1[i]<0.2)] 
        col_filtered_indx = col_filtered.index
        IS_filtered = df1["Ranked Sum"][df1["Ranked Sum"].index.isin(col_filtered_indx)] # filter IS only in the indx list

        corr_result_filtered = stats.spearmanr(col_filtered, IS_filtered, axis=0) 
        corr_filtered.corr_val[i] = corr_result_filtered[0]
        corr_filtered.p_val[i] = corr_result_filtered[1]

        loss_samples[i] = col_filtered_indx

    corr_filtered = corr_filtered.reset_index().rename(columns = {"index": "Gene-Symbol"})
    corr_unfiltered = corr_unfiltered.reset_index().rename(columns = {"index": "Gene-Symbol"})

    return corr_unfiltered, corr_filtered, loss_samples, col_filtered, col_filtered_indx, IS_filtered # call each df using index [0] or [1]

def corr_pvals_RNA_IS(df1,df2):
    '''Function to get correlation values and pvals (2 cols)
    Correlation involving IS will always have only 1 argument ie 1 df, since the IS is always in the last column
    - Filtering out samples that have gains for that gene
    - Output: 2 dfs with unfiltered and corr vals after filtered out gains'''

    # dataframes that will be output
    data_unfiltered = {'corr_val': None, 'p_val': None}
    data_filtered = {'corr_val': None, 'p_val': None}
    corr_unfiltered = pd.DataFrame(data_unfiltered, index=df1.columns[df1.columns != "Ranked Sum"])
    corr_filtered = pd.DataFrame(data_filtered, index=df1.columns[df1.columns != "Ranked Sum"])

    # create dict to store sample names that were filtered for less than 0.2
    keyList = [i for i in df1.columns[df1.columns != "Ranked Sum"]]
    loss_samples = {key: None for key in keyList}

    for i in keyList:
        # unfiltered
        corr_result_unfiltered = stats.spearmanr(df1[i], df1["Ranked Sum"], axis=0)
        corr_unfiltered.loc[i, 'corr_val'] = corr_result_unfiltered[0]
        corr_unfiltered.loc[i, 'p_val'] = corr_result_unfiltered[1]

        # filtered
        # filter out samples w gains ie only take samples that have losses or neutral
        col_filtered = df1[i][(df1[i] < 0.2)]
        col_filtered_indx = col_filtered.index
        loss_samples[i] = col_filtered_indx

        # filter RNAval_IS based on filtered samples from df1
        RNAval_IS_filtered = df2.loc[col_filtered_indx, :]

        # filter IS only in the indx list
        IS_filtered = df1["Ranked Sum"][df1["Ranked Sum"].index.isin(col_filtered_indx)]

        # calculate correlation between RNAval_IS_filtered and its corresponding "Ranked Sum" column
        corr_result_RNA_IS_filtered = stats.spearmanr(RNAval_IS_filtered[i], IS_filtered, axis=0)
        corr_filtered.loc[i, 'corr_val'] = corr_result_RNA_IS_filtered[0]
        corr_filtered.loc[i, 'p_val'] = corr_result_RNA_IS_filtered[1]

    corr_filtered = corr_filtered.reset_index().rename(columns={"index": "Gene-Symbol"})
    corr_unfiltered = corr_unfiltered.reset_index().rename(columns={"index": "Gene-Symbol"})


    return corr_unfiltered, corr_filtered, loss_samples, col_filtered, IS_filtered, col_filtered_indx

### MERGE FUNCTION ###
def get_merge_survival(df1, df2):
    '''Function to merge survival data geens with our genes'''
    merged = pd.merge(df1, df2, how ='inner', on =['Gene-Symbol'])
    return(merged)

### SURVIVAL FUNCTION ###
def get_survival_loss(mergedf):
    cox_loss = mergedf.loc[mergedf["items"] == "loss"]
    cox_loss = cox_loss[["Gene-Symbol", "z"]]
    return cox_loss

def get_survival_gain(mergedf):
    cox_gain = mergedf.loc[mergedf["items"] == "gain"]
    cox_gain = cox_gain[["Gene-Symbol", "z"]]
    return cox_gain

def get_survival_exp(mergedf):
    cox_exp = mergedf.loc[mergedf["items"] == "exp.continuous"]
    cox_exp = cox_exp[["Gene-Symbol", "z"]]
    return cox_exp
# %%
############### DATA MANIPULATION: KEEPING GENE SYMBOLS ###############

# We can start by filtering out all genes in the 'Gene-Symbol' column. We need these to create the final full data set. 
# Let's start by using creating a "key" of the columns we want to keep. In our case we only want to keep the 'Gene-Symbol' column.
to_keep = ['Gene-Symbol'] # Key
genes = all_genes[to_keep] # This will output a single column with all the 16K gene symbols that we need to build our metadata from.
# %%
############### DATA MANIPULATION: MERGING DATA FRAMES ###############

##### MERGING 'all_genes' gene names and 'SCNA_absolute' based on 'Gene-Symbol' column #####
# Now lets merge those data frames with 'Gene-Symbol' columns.
fulldata = merge_dataframes(genes,SCNA_absolute,'Gene-Symbol')
# %%
############### DATA MANIPULATION: REMOVING DATA ###############

# We don't need the 'Loucs-ID' and 'Cytoband' columns from our data, so lets remove these columns.
columns_to_drop = ['Locus-ID', 'Cytoband'] # Key to define the name of the columns that we want to drop. 
fulldata = drop_columns(fulldata, columns_to_drop)
# %%
############### DATA MANIPULATION: TRANSPOSING DATA ###############
# make list of colnames into strings
df_t = RNA.T
new_genenames = df_t.index.str.split("|").str[0] # split at "|" and only take the first part of the split
# change col names
RNA.columns = new_genenames

### CLEAN UP INDEX (SAMPLE ID) NAMES ####
#new_samplenames = df_t.index.str.split("|").str[0]

# Transpose and make RNA df into genes x Sample ID. Change index name to "Gene-Symbol"
RNA_T = RNA.T
RNA_T.index.name = "Gene-Symbol"
RNA_T = RNA_T.reset_index() # gene x sample ID
# %%
############### DATA MANIPULATION: MERGING DATA FRAMES BASED ON GENE-SYMBOL ###############
RNA_T = merge_dataframes(genes, RNA_T, 'Gene-Symbol')
# %%
############### DATA MANIPULATION: RESETING 'immune_score' INDEX ###############

# Manipulating 'immune_score' data frame so that we can merge it with 'RNA_T' data frame.
immune_score = immune_score.T
immune_score.index.name = 'Gene-Symbol'
immune_score.reset_index() # 8 rows x 521 columns
# %%
############### DATA MANIPULATION: RESETING 'RNA_T' INDEX ###############
RNA_T = RNA_T.set_index('Gene-Symbol')
RNA_T # 16466 rows x 250 columns
# %%
############### DATA MANIPULATION: EXTRACTING 'RANKED_SUM' DATA ###############
immune_score = immune_score.T
Ranked_Sum = immune_score[["Ranked Sum"]] # 520 rows x 1 columns
Ranked_Sum
# %%
############### DATA MANIPULATION: MERGING DATA FRAMES AND DROPING NA VALUES ###############
SCNA_IS_df0 = pd.concat([RNA_T.T, Ranked_Sum], axis = 1) # sample ID x genes

SCNA_IS_df = SCNA_IS_df0[SCNA_IS_df0.index.isin(immune_score.index)] 
# sort according to index (sample name)
SCNA_IS_df = SCNA_IS_df.sort_index(axis = 0)

## GETTING RID OF NA ROWS
# first, drop any rows with genes/columns having NA values >8000
SCNA_IS_df = SCNA_IS_df.dropna(thresh = 8000) # 520 rows
# next, drop rows where immune score is unavailable
SCNA_IS_df = SCNA_IS_df.dropna(subset = ['Ranked Sum']) # 512 rows

SCNA_IS_df # 520 rows x 16467 columns
# %%
############### DATA MANIPULATION: SELECTING GENES THAT ARE PRESENT IN BOTH DATA FRAMES ###############

SCNA_IS_df[SCNA_IS_df.index.isin(immune_score.index)] 
# sort according to index (sample name)
RNAval_IS = SCNA_IS_df.sort_index(axis = 0)

SCNA_IS_df # 520 rows x 16467 columns
RNAval_IS # 520 rows x 16467 columns
# %%
############### CODE TO CHECK IF THERE ARE ANY MISSING VALUES LEFT ###############

# check if column 'Ranked Sum' has any missing values or values equal to zero
if SCNA_IS_df['Ranked Sum'].isna().any() or (SCNA_IS_df['Ranked Sum'] == 0).any():
    print("Column Ranked Sum has missing values or values equal to zero.")
else:
    print("Column Ranked Sum does not have missing values or values equal to zero.")
# %%
############### SCNA_IS ###############

# do same for SCNA data
SCNA_only = merge_dataframes(SCNA_absolute,all_genes,'Gene-Symbol')
# Resetting Index for 'SCNA_only' data frame 
SCNA_only = SCNA_only.set_index('Gene-Symbol')
# We don't need the 'Loucs-ID' and 'Cytoband' columns from our data, so lets remove these columns.
columns_to_drop = ['Locus-ID', 'Cytoband'] # Key to define the name of the columns that we want to drop. 
SCNA_only = drop_columns(SCNA_only , columns_to_drop)

# merge with ranked sum 
SCNA_IS_df0 = pd.concat([SCNA_only.T, Ranked_Sum], axis = 1) # sample ID x genes

SCNA_IS_df = SCNA_IS_df0[SCNA_IS_df0.index.isin(immune_score.index)] 
# sort according to index (sample name)
SCNA_IS_df = SCNA_IS_df.sort_index(axis = 0)

## GETTING RID OF NA ROWS
# first, drop any rows with genes/columns having NA values >8000
SCNA_IS_df = SCNA_IS_df.dropna(thresh = 8000) # 520 rows
# next, drop rows where immune score is unavailable
SCNA_IS = SCNA_IS_df.dropna(subset = ['Ranked Sum']) # 512 rows
# sort according to index (sample name)
SCNA_IS = SCNA_IS.sort_index(axis = 0) # axis=0 is according to row

SCNA_IS # 512 rows x 16468
# %%
############### CODE TO CHECK IF THERE ARE ANY MISSING VALUES LEFT ###############

# check if column 'Ranked Sum' has any missing values or values equal to zero
if SCNA_IS['Ranked Sum'].isna().any() or (SCNA_IS['Ranked Sum'] == 0).any():
    print("Column Ranked Sum has missing values or values equal to zero.")
else:
    print("Column Ranked Sum does not have missing values or values equal to zero.")

# %% 
############### CODE TO CHECK IF THERE ARE ANY LOSSES (>0.2) ###############
if (SCNA_IS > 0.2).any().any():
    print("SCNA_IS does have values above 0.2")
else:
    print("SCNA_IS data frame does not have values above 0.2")

# %%
############### SCNA_IS VS RNAval_IS ###############

# subset SCNA only samples and genes that are present in both
# do first for SCNA table
SCNA_IS = SCNA_IS[SCNA_IS.index.isin(RNAval_IS.index)] # for samples
SCNA_IS = SCNA_IS[SCNA_IS.columns[SCNA_IS.columns.isin(RNAval_IS.columns)]] # for genes

# then do for RNA table
RNAval_IS = RNAval_IS[RNAval_IS.index.isin(SCNA_IS.index)] 
RNAval_IS = RNAval_IS[RNAval_IS.columns[RNAval_IS.columns.isin(SCNA_IS.columns)]]

SCNA_IS # 512 rows x 16467 columns
RNAval_IS #512 rows x 16467 columns
# %%
############### DATA MANIPULATION: SORTING BY ALPHABETICAL ORDER ###############

SCNA_IS = SCNA_IS.sort_index(axis=1)
RNAval_IS = RNAval_IS.sort_index(axis=1)

SCNA_IS # 512 rows x 16467 columns
RNAval_IS #512 rows x 16467 columns

# %%
############### CALCULATING CORRELATIONS: SORTING BY ALPHABETICAL ORDER ###############

# correlations
SCNA_IScorr = corr_pvals_SCNA_IS(SCNA_IS)
SCNA_RNAcorr = corr_pvals_SCNA_RNA(SCNA_IS, RNAval_IS) # this function filters both DNA and RNA similarly to get same samples
RNA_IScorr = corr_pvals_RNA_IS(SCNA_IS,RNAval_IS)
#%%
############### IMPORTING REQUIRED DATA FOR SURVIVAL ANALYSIS ###############
cox = pd.read_csv(f"E:/Immune_Evasion_Project/Mario/Cancer_Types/{cancer_type}/Data/{cancer_type}.COX.HR.p.txt", sep="\t")
all_genes = pd.read_csv("E:/Immune_Evasion_Project/Data/all_genes/DEL_HNSC.SKCM.BLCA.ESCA.LUAD.LUSC.PAAD_filtereddel.csv") # gene x sample
# %%
############### DATA MANIPULATION: RENAMING COLUMN ###############
cox = cox.rename(columns = {"gene": "Gene-Symbol"})
# %%
############### DATA MANIPULATION: MERGING DATA ###############
cox_merge = get_merge_survival(cox, all_genes)

# %%
############### CALCULATING Z-SCORES ###############

cox_loss = get_survival_loss(cox_merge)
cox_gain = get_survival_gain(cox_merge)
cox_exp = get_survival_exp(cox_merge)
# %%
#####  FINAL FULL DATA SET #####
model = merge_dataframes(SCNA_IScorr[1], RNA_IScorr[1], 'Gene-Symbol')
model = merge_dataframes(model,SCNA_RNAcorr[1], 'Gene-Symbol')
model = merge_dataframes(model,cox_loss, 'Gene-Symbol')
model = merge_dataframes(model,cox_exp, 'Gene-Symbol')
model = model.iloc[:, [0, 1, 3,5,7,8]]
model.columns = ['Gene-Symbol','Corr_SCNA_IS','Corr_RNA_IS', 'Corr_SCNA_RNA', 'Z-score_Loss', 'Z-score_Exp']
model
# %%
# Converting DataFrame into 'CSV' file
model.to_csv(f'E:/Immune_Evasion_Project/Mario/Output/{cancer_type}_model.csv', index = False)