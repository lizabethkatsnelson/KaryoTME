#%%
############### IMPORTING NEEDED MODULES FOR THE SCRIPT ###############
from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statistics
import sys


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, spearmanr, pearsonr

sys.path.append("E:/Immune_Evasion_Project/Mario/Scripts_Python/")
# %%
# We will import the 'Functions' library from our directory, make sure the file is present in the defined Directory above.
import Functions_Mario
from Functions_Mario import *

# This is for reloading the Library above. 
# import importlib
# importlib.reload(Functions_Mario)
# %%
from Functions_Mario import plot_controls

# %%
############### LOADING THE DATA ###############
if __name__ == "__main__":
    ## Get a list of CSV file names in the current directory
    # csv_files = [f for f in os.listdir(".") if f.endswith(".csv")] #If you have a directory full of 'csv' files that you want to run, you can use this code. Else, use the code below.
    csv_files  = [
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "Danias_allmetadata.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Young.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Dubrot_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Dubrot.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_Young_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_SKCM_Dubrot.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_Young.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_KearneyMC38.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_SKCM_Dubrot.csv"},

        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_KearneyMC38_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_PanOTI_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_Pan_Pmenl1_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_Dubrot_50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_COADREAD_Dubrot_Top50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_SKCM_Dubrot_Top50.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_SKCM_AP2_MC38.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_COADREAD_AP2_MC38.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_AP2_MC38.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_PAAD_AP2_MC38_Law.csv"},
        
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_BRCA_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_COADREAD_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_KIRC_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_KIRP_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_LUAD_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_LUSC_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_PAAD_AP2_MEL.csv"},
        {"path": "E:/Immune_Evasion_Project/Output/", "name": "DataControls_SKCM_AP2_MEL.csv"}
    ] #Use this code if you have a list of specific csv files that you want to use. You need to specify the path and the name of each csv file. Adjust accordingly.

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
