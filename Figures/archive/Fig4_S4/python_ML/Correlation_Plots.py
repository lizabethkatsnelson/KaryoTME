# %%
import matplotlib.pyplot as plt
import numpy as np

### IMPORTING ALL NECESSARY MODULES

from operator import index
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import statistics
import sys
sys.path.append("E:/Immune_Evasion_Project/Scripts")

# LogReg modules
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import roc_auc_score, f1_score
import sklearn.metrics as metrics
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV
from sklearn.preprocessing import scale 
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing
# %%

## Definining Control List
HNSC_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_SKCM_Dubrot_50.csv")
corr = HNSC_metadata['Corr_RNA_IS']
labels = HNSC_metadata['Label']
# %%
# Assume that the correlation measurement has been calculated and stored in a variable called 'corr'
# Assume that the positive and negative controls are labeled and stored in a variable called 'labels'

# %%
# Create a figure and axes object
fig, ax = plt.subplots()

# Plot the correlation measurement as a scatter plot
ax.scatter(range(len(corr)), corr, c=labels)

# Add labels and title
ax.set_xlabel('Data Points')
ax.set_ylabel('Correlation')
ax.set_title('BRCA_Dubrot Top 50: Plot of Correlation Data from SCNA vs RNA with Positive/Negative controls')

# Add a colorbar
cbar = plt.colorbar()
cbar.set_label('Control')

# Show the plot
plt.show()
# %%
