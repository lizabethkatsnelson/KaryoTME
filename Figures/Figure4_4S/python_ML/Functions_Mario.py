# %%
############### IMPORTING NEEDED MODULES FOR THE SCRIPT ###############
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, spearmanr, pearsonr
# %%
############### DEFINING '-HPV' FUNCTION: ###############
def HPVneg(HNSC_df):
    '''Get only HPV- samples. Where it's comparing it from has to be chanegd manually'''
    HNSC_neg = pd.read_csv("E:/Immune_Evasion_Project/Mario/Cancer_Types/HNSC/Data/HNSC_neg_list.csv")
    # subset only HPV- columns/samples
    HNSC_neg = HNSC_df[HNSC_df.columns[HNSC_df.columns.isin(HNSC_neg)]] # CHANGE HERE
    HNSC_hpvneg = pd.concat([HNSC_df, HNSC_neg], axis = 1) #.iloc[:,0:3]
    return HNSC_hpvneg
# %%
############### DEFINING MERGING FUNCTION 1: ###############
def merged_df(df1, df2):
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
# %%
############### DEFINING MERGING FUNCTION 2: ###############
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
# %%
############### DEFINING 'CORRELATION' FUNCTION: ###############

### The code below is to get an example and get an 'output' for the 'calculate_correlation' function defined below. ###

# data = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [2, 3, 4, 5]})
# correlation_type = 'spearman'
# correlation, p_value = calculate_correlation(data, correlation_type)
# print(f'Correlation: {correlation}, P-value: {p_value}')

def calculate_correlation(data, correlation_type):
    '''In this script, the calculate_correlation function takes two arguments: a DataFrame of the data, and a string indicating the type of correlation to calculate ('spearman' or 'pearson'). The function uses an if-elif block to determine which correlation to calculate and returns the correlation coefficient and the p-value. You can call the function and pass the data and the correlation type to get the correlation and p-value.'''
    if correlation_type == 'spearman':
        correlation, p_value = spearmanr(data)
    elif correlation_type == 'pearson':
        correlation, p_value = pearsonr(data)
    else:
        raise ValueError("Invalid correlation type")
    return correlation, p_value
# %%
############### DEFINING 'MERGE' FUNCTION: ###############

# Merge function of two data frames based on a name of a column
def merge_dataframes(df1, df2, on):
    return pd.merge(df1, df2, on=on)
# %%
############### DEFINING 'DROP COLUMNS' FUNCTION: ###############

# Droping columns based on their name from the data frame. 
def drop_columns(df, columns):
    return df.drop(columns, axis=1)
# %%
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
# %%
############### DISTRIBUTION HISTOGRAMS FUNCTION ###############
def plot_controls(data, file_name):
    # Create separate DataFrames for positive and negative controls
    positive_controls = data[data['Label'] == 1].iloc[:, 1:]
    negative_controls = data[data['Label'] == 0].iloc[:, 1:]

    # Create a figure with subplots for each feature
    num_features = len(positive_controls.columns)
    fig, axs = plt.subplots(nrows=1, ncols=num_features, figsize=(15, 3))

    # Loop over each feature and plot the distribution for positive and negative controls
    for i, feature in enumerate(positive_controls.columns):
        # Plot the histograms for positive and negative controls
        axs[i].hist(positive_controls[feature], bins=20, alpha=0.5, density=True, label='Positive controls')
        axs[i].hist(negative_controls[feature], bins=20, alpha=0.5, density=True, label='Negative controls')

        # Plot the normal distribution curve for positive and negative controls
        x_min = np.min(data[feature])
        x_max = np.max(data[feature])
        x_range = np.linspace(x_min, x_max, 100)
        mu_positive = positive_controls[feature].mean()
        sigma_positive = positive_controls[feature].std()
        mu_negative = negative_controls[feature].mean()
        sigma_negative = negative_controls[feature].std()
        axs[i].plot(x_range, norm.pdf(x_range, mu_positive, sigma_positive), 'r-', label='Positive control curve')
        axs[i].plot(x_range, norm.pdf(x_range, mu_negative, sigma_negative), 'b-', label='Negative control curve')

        axs[i].set_xlabel(feature)
        axs[i].legend(loc='center', bbox_to_anchor=(0.5, -0.4))

    # Add the file name to the plot title in bold
    fig.suptitle(f"\n\nDistribution of positive and negative controls for '{file_name}'", fontweight="bold", y=2.2)

    # Show the plot
    plt.tight_layout()
    plt.subplots_adjust(top=1.85, bottom=1.)
    plt.show()

# In this script, the os.listdir() function is used to get a list of CSV file names in the current directory. 
# The script then loops over each file, loads it into a pandas DataFrame using the pd.read_csv() function, and plots the distribution for positive and negative controls using the plot_controls() function. 
# You can modify the plot_controls() function and the CSV files in the directory to suit your needs.
# %%
############### DEFINING CSV FILES TO USE ###############
if __name__ == "__main__":
    ## Get a list of CSV file names in the current directory
    # csv_files = [f for f in os.listdir(".") if f.endswith(".csv")] #If you have a directory full of 'csv' files that you want to run, you can use this code. Else, use the code below.
    csv_files  = [
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "Danias_allmetadata.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Young.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Dubrot_50.csv"}
    ] #Use this code if you have a list of specific csv files that you want to use. You specify the path and the name of each csv file.

    # Loop over each CSV file and plot the distribution for positive and negative controls
    for csv_file in csv_files:
        # Load the data into a pandas DataFrame
        data = pd.read_csv(csv_file["path"] + "/" + csv_file["name"])

        # Plot the distribution for positive and negative controls
        plot_controls(data, csv_file["name"])

# The if __name__ == "__main__": statement is not part of the plot_controls() function. It's a separate part of the script that is executed when the script is run as the main program.
# The if __name__ == "__main__": statement is a common Python idiom that allows the script to be used as both a standalone program and a module that can be imported into another script. 
# When the script is run as the main program, the code within the if __name__ == "__main__": block is executed. 
# If the script is imported into another script, the code within the if __name__ == "__main__": block is not executed.
# In this script, the if __name__ == "__main__": block is used to get a list of CSV file names in the current directory, load each file into a pandas DataFrame, and plot the distribution for positive and negative controls for each file. 
# When you run the script, the code within the if __name__ == "__main__": block is executed, which calls the plot_controls() function for each CSV file in the directory.
# %%
