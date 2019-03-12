from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pickle
#The goal of this Rscript is to compare variance within vs variance between

conditions = pd.read_csv("/ebio/abt6_projects9/tnseq/data/fitness_datasets/all_conditions.txt", sep="\t")

#this is a dictionary with the set and gene cluster as the key and the fitness measurement as a value
full_tag_dict = pickle.load(open('/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.cpk', 'rb'))

'''full_tag_pd = pd.DataFrame(columns=set([line[0] for line in list(full_tag_dict.keys())]), index=set([line[1] for line in list(full_tag_dict.keys())]))
for col in full_tag_pd.columns:
    for index in full_tag_pd.index:
        try:
            full_tag_pd[col].ix[index] = full_tag_dict[(col,index)]
            print("SUCCESS")
        except KeyError:
            full_tag_pd[col].ix[index] = "NaN"
full_tag_pd.to_csv('/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.csv', sep="\t")
'''
full_tag_pd = pd.read_csv('/ebio/abt6_projects9/tnseq/data/fitness_datasets/fitness_tables/full_tag_dict_file.csv', sep="\t", index_col=0)

#
def give_shared(conditions):
    '''Find experiments in which the same organism was measured in the same condition more than once'''
    shared = {}
    for strain in set(conditions['Organism']):
        for condition in set(conditions['Condition']):
            subset=[]
            subset=conditions.loc[(conditions['Organism'] == strain) & (conditions['Condition'] == condition)]
            if len(subset)>1:
                shared["_".join((strain+"__"+str(condition)).split(" "))] = list(subset.Name)
            else:
                pass
    return shared

def give_shared_between(conditions):
    '''Find experiments in which the same condition was measured more than once'''
    shared_between = {}
    for condition in set(conditions['Condition']):
            subset = conditions.loc[conditions['Condition'] == condition]
            if len(subset)>1:
                shared_between[condition] = list(subset.Name)
    return shared_between


def calc_corr_within(shared, full_tag_pd)
    '''Calculate correlation between fitness of same gene in 'same condition'''
    #full correlation matrix
    corr_pd = full_tag_pd.astype(float).corr()
    corr_dict = {}
    group_dict={}
    #first find relevant subsets
    for key, group in shared.items():
        temp=[]
        for set1 in group:
            for set2 in group:
                if set1 != set2:
                    try:
                        temp.append(corr_pd[set1].loc[set2])
                    except KeyError:
                        pass
        corr_dict[key] = np.mean(temp)
        for val in group:
            group_dict[val]=np.mean(temp)
    group_dict=dict([rec for rec in group_dict.items() if not math.isnan(rec[1])])
    return group_dict

def correct_large_matrix(group_dict, full_tag_pd):
'''this function takes the correlation matrix and corrects'''
    corr_pd = full_tag_pd.astype(float).corr()
    corr_normalized=pd.DataFrame(columns=[rec for rec in group_dict.keys() if rec in corr_pd.index], index=[rec for rec in group_dict.keys() if rec in corr_pd.index])
    for col in corr_normalized.columns:
        col_cor="NaN"
        ind_cor="NaN"
        yes="No"
        for ind in corr_normalized.index:
            if col in list(group_dict.keys()):
                col_cor = group_dict[col]
                yes="yes"
            else:
                col_cor = "NaN"
            if ind in list(group_dict.keys()):
                ind_cor = group_dict[ind]
                yes="yes"
            else:
                ind_cor = "NaN"
            if yes=="yes":
                denom = np.sqrt(ind_cor^2*col_cor^2)
                corr_normalized[col].loc[ind]=corr_pd[col].loc[ind]/denom
            else:
                corr_normalized[col].loc[ind]="NaN"
    return corr_normalized





#1424 conditions have some repetition
shared = give_shared(conditions)
group_dict = calc_corr_within(shared, full_tag_pd)
corr_normalized = correct_large_matrix(group_dict, full_tag_pd)


#calculate average correlation within treatment and strain
calc_corr_witin(shared, full_tag_pd)


#build a dictionary for every set the average correlation


#now we want to look at the relationship across organisms
