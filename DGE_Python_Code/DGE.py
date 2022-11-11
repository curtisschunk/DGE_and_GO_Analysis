# Import Packages
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy import stats
import argparse

# Inputs
parser = argparse.ArgumentParser()
parser.add_argument('-c','--control',required=True,type=str,help="File for Control Gene Expression")
parser.add_argument('-t','--test',required=True,type=str,help="File for Test Gene Expression")
parser.add_argument('-g','--gene_list',required=True,type=str,help="File for list of Genes")
parser.add_argument('-n','--n_dif_genes',required=False,default=200,type=int,help="Number of Differential Genes to Export")
parser.add_argument('-s','--save_path',required=False,default='Differential_Genes',type=str,help="Save Path and Save File Name")
args = parser.parse_args()

## Functions ##

def test(data1, data2):
    stat, p = stats.mannwhitneyu(data1, data2)
    return p

#for specified gene, finds p-val and whether it is over or under expressed

def p_val_finder(test_array,ctrl_array,gene_test_num):
    
    test_vals=test_array[gene_test_num,:]
    ctrl_vals=ctrl_array[gene_test_num,:]

    p_val=test(test_vals,ctrl_vals)
    log_p_val=-np.log10(p_val)

    if np.mean(test_vals)>np.mean(ctrl_vals): #overexpressed gene in test group
        cur=True
    else: #underexpressed gene in test group
        cur=False

    return p_val,log_p_val,cur


# cur -> if = 1, expression is higher in the group being compared to the other groups
#        if = 0, expression is lower in the group being comared to the other groups


def dif_gene_analysis(test_array,ctrl_array,gene_list): 
    gene_len=len(gene_list)
    dif_gene_array=np.zeros([gene_len,3])
    dif_gene_list=list()
    for i in range(gene_len):

        dif_gene_list.append(gene_list[i])
        dif_gene_array[i,0],dif_gene_array[i,1],dif_gene_array[i,2]=p_val_finder(test_array,ctrl_array,i)
        
    tot_df=pd.DataFrame([dif_gene_list,dif_gene_array[:,0],dif_gene_array[:,1],dif_gene_array[:,2]]).transpose()
    tot_df=tot_df.set_axis(['Gene','p-val','-log10(p-val)','cur'], axis=1, inplace=False)
    
    sorted_df=tot_df.sort_values(by=['cur', '-log10(p-val)'],ascending=False)
    
    top200_genes=sorted_df['Gene'][0:args.n_dif_genes]
    
    #return sorted_df
    return top200_genes

# Import Files
ctrl_array=np.array(pd.read_csv(args.control))
test_array=np.array(pd.read_csv(args.test))
genes_array=np.array(pd.read_csv(args.gene_list))

# Run Differential Gene Expression Analysis
dif_genes=dif_gene_analysis(test_array,ctrl_array,genes_array)

# Fix Labeling
dif_genes=np.array(dif_genes)

fixed_label_dif_genes=list()
for i in range(len(dif_genes)):
    current_gene_working_with=str(dif_genes[i])
    current_gene_working_with=current_gene_working_with.replace("['","")
    current_gene_working_with=current_gene_working_with.replace("']","")
    fixed_label_dif_genes.append(current_gene_working_with)

# Export List of Differential Genes
pd.DataFrame(fixed_label_dif_genes).to_csv(args.save_path+'.csv',index=None)

