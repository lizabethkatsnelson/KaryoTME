## FUNCTIONS ##
import pandas as pd
import numpy as np
import glob
from operator import index
import matplotlib.pyplot as plt
from scipy import stats
import statistics

# Function to get merged df on Gene-Symbol - only keep genes in gene list, then transpose
def get_merged_df(df1, df2):
    ''' Function to only keep genes from TCGA dataset that are in our gene list, 
    then transpose. Merged df will be (sample ID x genes). Merge on "Gene-Symbol".
    df1 = TCGA dataset (genes x patient)
    df2 = Smaller dataset (genes x patient) 
    Output  = sample x genes'''

    cols_to_throw = ["Locus-ID", "Cytoband"]
    if sum(df1.columns.isin(cols_to_throw)) > 0: # isin values
        df1 = df1.drop( cols_to_throw, axis = 1) # make if not exist, skip this
    merge = pd.merge(df1, df2, how ='inner', on =['Gene-Symbol'])
    merge = merge.iloc[:,0:-17]
    merged_df = merge.set_index("Gene-Symbol").T

    # find n of genes that did not merge
    not_merged = len(df1[~df1["Gene-Symbol"].isin(merged_df.columns)]) 
    print("n of genes that did not merge:", not_merged)

    # find genes that were duplicated
    dup_index = merged_df.columns.duplicated()
    dup_index = np.where(dup_index == True) # get index of genes /columns that are duplicated
    print("n of genes that have duplicates:", len(dup_index[0][:])) # length of genes that are duplicated
    
    col_list = merged_df.columns
    genes_duplicated = []
    # list of genes that were duplicated
    for i in dup_index:
        genes_duplicated.append(col_list[i])
    print(genes_duplicated)

    # drop rows that the genes are duplicated (but we wont know which were dropped)
    merged_df_T = merged_df.T.reset_index().drop_duplicates(subset=['Gene-Symbol'])
    merged_df = merged_df_T.set_index('Gene-Symbol').T
    return merged_df

def get_HPVneg(HNSC_df):
    '''Get only HPV- samples. Where it's comparing it from has to be chanegd manually'''
    HNSC_neg = pd.read_csv("E:/Immune_Evasion_Project/Mario/Cancer_Types/HNSC/Data/HNSC_neg_list.csv")
    # subset only HPV- columns/samples
    HNSC_neg = HNSC_df[HNSC_df.columns[HNSC_df.columns.isin(HNSC_neg)]] # CHANGE HERE
    HNSC_hpvneg = pd.concat([HNSC_df, HNSC_neg], axis = 1) #.iloc[:,0:3]
    return HNSC_hpvneg

def merge_IS_to_df(parameter_df, IS_df):
    ''' Function to add Immune score column to the parameter df. Any samples with >8000 
    NA for genes and any samples with no Immune score (NA) will be dropped.
    
    parameter_df: gene x sample
    IS_df: sample x 8 cols '''

    df_t = parameter_df.T.reset_index()
    df_t = df_t.set_index("Gene-Symbol") # gene x sample

    Ranked_Sum = IS_df[["Ranked Sum"]] # 520 samples x 8
    
    # merge with ranked sum 
    SCNA_IS_df0 = pd.concat([df_t.T, Ranked_Sum], axis = 1) # sample ID x genes

    SCNA_IS_df = SCNA_IS_df0[SCNA_IS_df0.index.isin(IS_df.index)] 
    # sort according to index (sample name)
    SCNA_IS_df = SCNA_IS_df.sort_index(axis = 0)

    ## GETTING RID OF NA ROWS
    # first, drop any rows with genes/columns having NA values >8000
    SCNA_IS_df = SCNA_IS_df.dropna(thresh = 8000) # 520 rows
    # next, drop rows where immune score is unavailable
    SCNA_IS_df = SCNA_IS_df.dropna(subset = ['Ranked Sum']) # 512 rows
    return SCNA_IS_df

def get_freq_del(SCNA_table, cutoff_val):
    '''Function to get freq of deletion for a df. 
    rows are patients, columns are genes'''
    n_del_pergene = {} # list of n of patients with deletion for EACH GENE
    pct_del_pergene = {}
    for i in SCNA_table.columns:
        deletions = SCNA_table[ (SCNA_table[i]) < cutoff_val]
        n_pct = {i:len(deletions)/len(SCNA_table)}
        n_del = {i:len(deletions)}

        n_del_pergene.update(n_del)
        pct_del_pergene.update(n_pct)
        
    freq_del = pd.DataFrame(columns = ("Gene-Symbol", "freq_del"))
    freq_del["Gene-Symbol"] = pct_del_pergene.keys()
    freq_del["freq_del"] = pct_del_pergene.values()
    freq_del = freq_del.set_index("Gene-Symbol")
        
    return freq_del

def get_corr_pvals_IS(df1):
    '''Function to get correlation values and pvals (2 cols)
    Correlation involving IS will always have only 1 argument ie 1 df, since the IS is always in the last column
    - Input: sample x genes (+ Ranked Sum col)
    - Filtering out samples that have gains for that gene
    - Output: 2 dfs with unfiltered and corr vals after filtered out gains'''

    data_unfiltered = {'corr_val':None, 'p_val':None} 
    data_filtered = {'corr_val':None, 'p_val':None}

    # dataframes that will be output
    corr_unfiltered = pd.DataFrame(data_unfiltered, index = df1.columns[:-1])
    corr_filtered = pd.DataFrame(data_filtered, index = df1.columns[:-1])

    for i in df1.columns[:-1]: # don't take the last column (ranked sum col)
        # unfiltered
        corr_result_unfiltered = stats.spearmanr(df1[i], df1["Ranked Sum"], axis=0) #If axis=0 (default), then each column represents a variable
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

    corr_filtered = corr_filtered.reset_index().rename(columns = {"index": "Gene-Symbol"})
    corr_unfiltered = corr_unfiltered.reset_index().rename(columns = {"index": "Gene-Symbol"})

    return corr_unfiltered, corr_filtered # call each df using index [0] or [1]

def get_corr_pvals_IS(df1):
    '''Function to get correlation values and pvals (2 cols)
    Correlation involving IS will always have only 1 argument ie 1 df, since the IS is always in the last column
    - Filtering out samples that have gains for that gene
    - Output: 2 dfs with unfiltered and corr vals after filtered out gains'''
    data_unfiltered = {'corr_val':None, 'p_val':None} 
    data_filtered = {'corr_val':None, 'p_val':None}

    # dataframes that will be output
    corr_unfiltered = pd.DataFrame(data_unfiltered, index = df1.columns[:-1])
    corr_filtered = pd.DataFrame(data_filtered, index = df1.columns[:-1])

    # create dict to store sample names that were filtered for <-0.2
    keyList = []
    for i in df1.columns[:-1]:
        keyList.append(i)
    loss_samples = {key: None for key in keyList}

    for i in df1.columns[:-1]: # don't take the last column (ranked sum col)
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

    return corr_unfiltered, corr_filtered,loss_samples # call each df using index [0] or [1]

def cleanup_RNA(RNA_df):
    ### CLEAN UP GENE NAMES ####
    # make list of colnames into strings
    df_t = RNA_df.T
    new_genenames = df_t.index.str.split("|").str[0] # split at "|" and only take the first part of the split
    # change col names
    RNA_df.columns = new_genenames

    ### CLEAN UP INDEX (SAMPLE ID) NAMES ####
    #new_samplenames = df_t.index.str.split("|").str[0]

    # Transpose and make RNA df into genes x Sample ID. Change index name to "Gene-Symbol"
    RNA_df_T = RNA_df.T
    RNA_df_T.index.name = "Gene-Symbol"
    RNA_df_T = RNA_df_T.reset_index() # gene x sample ID

    return RNA_df_T

def match_sample_genes(df1, df2):
    '''Function to ensure samples and genes that are in 2 dfs match and have the sample number. Here usage is for SCNA_IS df
    and SCNA_IS df. Do before correlation between the 2 dfs.
    Output: sample x genes'''
    # subset SCNA only samples and genes that are present in both
    # do first for SCNA table
    df1 = df1[df1.index.isin(df2.index)] # for samples
    df1_new = df1[df1.columns[df1.columns.isin(df2.columns)]].drop_duplicates().sort_index() # for genes. Sort alphabetically
    df1_new = df1_new.sort_index(axis =1)
    # then do for RNA table
    df2 = df2[df2.index.isin(df1.index)] 
    df2_new = df2[df2.columns[df2.columns.isin(df1.columns)]].drop_duplicates().sort_index()
    df2_new = df2_new.sort_index(axis =1)


    return df1_new, df2_new # correlation test will index the output of this as [0] or [1]

def get_corr_2dfs(df1, df2):
    ''' Function to get correlation between 2 dataframes - mainly SCNA v RNA in TCGA
    df1 = SCNA with IS
    df2 = RNA with IS '''

    data_unfiltered = {'corr_val':None, 'p_val':None} 
    data_filtered = {'corr_val':None, 'p_val':None}
    # dataframes that will be output
    corr_unfiltered = pd.DataFrame(data_unfiltered, index = df1.columns[:-1])
    corr_filtered = pd.DataFrame(data_filtered, index = df1.columns[:-1])

    for i in df1.columns[:-1]:
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
    
    return corr_unfiltered, corr_filtered, df2_filtered # call each df using index [0] or [1]

### CCLE FUNCTIONS ###

def CCLE_cleanup(CCLE_df):
    '''Function to cleanup gene names in CCLE data. Input should be RNA exp df'''
    new_genenames = CCLE_df.columns.str.split(" ").str[0] # split at " " and only take the first part of the split
    CCLE_df.columns = new_genenames
    return CCLE_df

def CCLE_filter(CCLE_df, sample_info, disease, disease_subtype):

    # Get IDs first
    xtumor_only = sample_info[sample_info.disease.str.contains(f'^{disease}') & sample_info.disease_subtype.str.contains(f"{disease_subtype}", regex = True, na=False)]
    if disease == "Head and Neck Cancer":
        CCLE_HPVneglist = pd.read_excel("/Users/daniaannuar/Documents/Thesis/HNSC/HNSC_CELL_HPVneg.xlsx")
        xtumor_only = sample_info[sample_info.stripped_cell_line_name.isin(CCLE_HPVneglist["Cell name"])]
    else:
        xtumor_only = sample_info[sample_info.disease.str.contains(f'^{disease}') & sample_info.disease_subtype.str.contains(f"{disease_subtype}", regex = True, na=False)]

    CCLE_IDs = xtumor_only.DepMap_ID.drop_duplicates()
    # Get only x type, make gene column name as Gene-Symbol
    CCLE_filter = CCLE_df[CCLE_df.iloc[:,0].isin(CCLE_IDs)] # 65 samples
    CCLE_filter_T = CCLE_filter.T
    CCLE_filter_T.columns = CCLE_filter.T.iloc[0,:] # change column names 
    CCLE_filter_T = CCLE_filter_T.iloc[1:,:] # get rid of 1st row w column name
    CCLE_filter_T.index.name = "Gene-Symbol"
    CCLE_filter_T = CCLE_filter_T.reset_index()

    return CCLE_filter_T

def CCLE_cleanup_prot(CCLE_prot):
    CCLE_prot_T = CCLE_prot.T
    CCLE_prot_T.columns = CCLE_prot.T.iloc[0,:] # change column names 
    CCLE_prot_T = CCLE_prot_T.iloc[1:,:].reset_index()
    
    #CCLE_prot = CCLE_df.rename(columns = {"Gene_Symbol": "Gene-Symbol"})
    #CCLE_prot1 = CCLE_prot.reset_index().T
    return CCLE_prot_T

## SURVIVAL FUNCTIONS ##

def get_merge_survival(df1, df2):
    '''Function to merge survival data geens with our genes'''
    merged = pd.merge(df1, df2, how ='inner', on =['Gene-Symbol'])
    return(merged)

def get_survival_loss(mergedf):
    HNSC_cox_loss = mergedf.loc[mergedf["items"] == "loss"]
    HNSC_cox_loss = HNSC_cox_loss[["Gene-Symbol", "z"]]
    return HNSC_cox_loss

def get_survival_gain(mergedf):
    HNSC_cox_gain = mergedf.loc[mergedf["items"] == "gain"]
    HNSC_cox_gain = HNSC_cox_gain[["Gene-Symbol", "z"]]
    return HNSC_cox_gain

def get_survival_exp(mergedf):
    HNSC_cox_exp = mergedf.loc[mergedf["items"] == "exp.continuous"]
    HNSC_cox_exp = HNSC_cox_exp[["Gene-Symbol", "z"]]
    return HNSC_cox_exp



## DISTRIBUTIONS ##

import matplotlib.pyplot as plt
import seaborn as sns
import statistics

def get_hist_corr(parameter_corr_table):
    for i in [parameter_corr_table]:
        plt.xlabel("Correlation values", size=10)
        plt.ylabel("n of genes", size=10)
        plt.hist(i.corr_val) # corr values
        sig_corr = i[i.p_val < 0.05]
        plt.hist(sig_corr.corr_val)
        plt.legend(["corr values", "corr values (p<0.05)"], loc='upper right')

def get_hist_freqdel(freq_del_col, label_name):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import statistics

    freq_del_cols = freq_del_col.iloc[:,1]
    #plt.figure(figsize=(8,6))
    plt.legend(loc='upper right')
    plt.xlabel("Freq deletion", size=14)
    plt.ylabel("n of genes", size=14)
    #plt.hist(freq_del_cols, alpha=0.5, label=freq_del_col)
    return plt.hist(freq_del_cols, alpha=0.5, label=label_name)

## Model building functions ##

# Concat params
def concat_params_new(params_list, params_names):
    params_list2 = []
    for n,i in zip(params_list, params_names):
        # for i in params_names:
        df = n.rename(columns = {n.columns[2]: i}) # rename the cols according to its parameter name
        df1 = df.iloc[:,1:3].set_index("Gene-Symbol") # set index, and only get the gene-symbol and param value columns
        params_list2.append(df1)
    print(params_list2)

    dfs_concat = pd.concat(params_list2, axis = 1)
    return(dfs_concat)

# Get neg/pos only df from all_params df
def get_control_df(pos_control_df, neg_control_df, all_params_df, *subset):
    '''Function to add class column to control genes and output the pos and neg control list
    It will also concatenate the controls into a "mix" df containing all your controls 
    so you can use for the model.
    Outputs:
    [0] neg control genes only df
    [1] pos control genes only df
    [2] df of all control genes and the predictors (x of model)
    [3] df of all control genes and their classes (0 or 1) (y of model) '''

    # get only neg control genes from the all_params df
    neg_only = neg_control_df.drop_duplicates(subset = "Gene-Symbol")
    neg_only = all_params_df[all_params_df.index.isin(neg_control_df["Gene-Symbol"])].fillna(0)
    neg_only["class"] = 0
    # get only pos control genes from the all_params df
    pos_only = pos_control_df.drop_duplicates(subset = "Gene-Symbol")
    pos_only = all_params_df[all_params_df.index.isin(pos_control_df["Gene-Symbol"])].fillna(0)
    pos_only["class"] = 1

    if subset == 79:
    # if you want to subset
        import random
        pos_subset_index = random.choices(pos_only.index, k = subset) # get 83 random indices
        pos_subset = pos_only.loc[pos_subset_index,:] # get the rows usinf index list
        pos_subset["class"] = 1

        neg_subset_index = random.choices(neg_only.index, k = subset) # get 83 random indices
        neg_subset = neg_only.loc[neg_subset_index,:] # get the rows usinf index list
        neg_subset["class"] = 0

        # if subseting, change "pos_only" to "pos_subset"
        mix = pd.concat([neg_subset, pos_subset])
        mix['class'] = pd.Categorical(mix['class'])
        mix_x = mix.iloc[:,:-1] # dont take class column
        mix_y = mix.iloc[:,-1] # only take class column
    else:
        mix = pd.concat([neg_only, pos_only])
        mix['class'] = pd.Categorical(mix['class'])
        mix_x = mix.iloc[:,:-1] # dont take class column
        mix_y = mix.iloc[:,-1] # only take class column

    return(neg_only, pos_only, mix_x, mix_y)

# Get log reg model
def get_model(x, y):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix
    from sklearn.model_selection import train_test_split
    from sklearn.model_selection import cross_val_score
    from sklearn.model_selection import cross_validate
    from sklearn.metrics import recall_score
    from sklearn.metrics import accuracy_score
    from sklearn.metrics import precision_score
    from sklearn.metrics import roc_auc_score
    import sklearn.metrics as metrics
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV
    
    # FOR CROSS VAL
    # getting train, test split for all data
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42, stratify = y)

    # fit the model with ALL variables and calculate auc score
    clf_logistic = LogisticRegression(penalty='l1', solver='liblinear').fit(x_train,y_train)
    #preds = clf_logistic.predict_proba(x_test)[:,1]
    #roc_auc_score(y_test, preds)

    lassocv = LassoCV(alphas = None, cv = 10, max_iter = 100000, normalize = True)
    lassocv.fit(X_train, y_train)

    lasso.set_params(alpha=lassocv.alpha_)
    lasso.fit(X_train, y_train)

    pred = clf_logistic.predict(x_test)
    print("proportion of predicted as class 0:", round(sum(pred == 0)/len(pred), 2))
    print("proportion of predicted as class 1:", round(sum(pred == 1)/len(pred), 2))
    print("accuracy score:", round(accuracy_score(pred, y_test), 2))
    # So our model has x% overall accuracy, but is it because it's predicting only 1 class?
    print("unique values in predicted classes are:", np.unique(pred)) # list unique values in the predicted values

    # using cross validation
    k_fold_res = cross_validate(clf_logistic, x_train, y_train, scoring=['precision', 'recall','accuracy'], cv=5)
    print("Precision: %0.2f \u00B1 %0.2f" % (k_fold_res["test_precision"].mean(), k_fold_res["test_precision"].std()))

    print("Recall %0.2f \u00B1 %0.2f"% (k_fold_res["test_recall"].mean(), k_fold_res["test_recall"].std()))

    print("Accuracy %0.2f \u00B1 %0.2f" % (k_fold_res["test_accuracy"].mean(), k_fold_res["test_accuracy"].std()))
    return pred, clf_logistic

