from pickle import FALSE
from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
import math 
import numpy as np
import tkinter as tk 
import random
np.seterr(divide = 'ignore') 

def df_add_stats(df_table):
    n_t = df_table.shape[0]
    
    col_sum = df_table.sum(axis = 0, numeric_only = True) #count up the number of times each mut appears and store it in a frame
    df_ret = col_sum.to_frame('n(t_L)')
    df_ret = df_ret.reset_index()
    df_ret.rename(columns={'index':'Mutation List'}, inplace=True)
    df_ret['n(t_R)'] =  n_t - df_ret['n(t_L)']
    
    df_c = df_table[df_table["Unnamed: 0"].str.contains(pat = "NC") == False].reset_index() #create a dataframe of just c
    df_c = df_c.drop("index", axis = 1)
    df_c_sum = df_c.sum(axis = 0).to_frame('n(t_L, C)') #count up the number of times each mut appears for c and store it in a frame
    df_c_sum = df_c_sum.reset_index()
    df_c_sum.rename(columns={'index':'Mutation List'}, inplace=True)
    
    df_ret = df_ret.merge(df_c_sum, on = "Mutation List")
    df_ret['n(t_R, C)'] = df_c.shape[0] - df_ret['n(t_L, C)']
    df_nc = df_table[df_table["Unnamed: 0"].str.contains(pat = "NC")].reset_index() #create a dataframe of just nc
    df_nc = df_nc.drop("index", axis = 1)
    
    df_nc_sum = df_nc.sum(axis = 0).to_frame('n(t_L, NC)') #count up the number of times each mut appears for nc and store it in a frame
    df_nc_sum = df_nc_sum.reset_index()
    df_nc_sum.rename(columns={'index':'Mutation List'}, inplace=True)
   
    df_ret = df_ret.merge(df_nc_sum, on = "Mutation List") # merge all of the tables
    df_ret['n(t_R, NC)'] = df_nc.shape[0] - df_ret['n(t_L, NC)']
    df_ret['P_L'] = df_ret['n(t_L)'] / n_t
    df_ret['P_R'] = df_ret['n(t_R)'] / n_t
    
    #P(C | tL) = n(tL, C) / n(tL)
    df_ret['C | t_L'] = df_ret['n(t_L, C)'] / df_ret['n(t_L)']
    #P(NC | tL) = n(tL, NC) / n(tL)
    df_ret['NC | t_L'] = df_ret['n(t_L, NC)'] / df_ret['n(t_L)']
    #P(C | tR) = n(tR, C) / n(tR)
    df_ret['C | t_R'] = df_ret['n(t_R, C)'] / df_ret['n(t_R)']
    #P(NC | tR) = n(tR, NC) / n(tR)
    df_ret['NC | t_R'] = df_ret['n(t_R, NC)'] / df_ret['n(t_R)']
    #Q(s|t) = |P(C | tL) - P(C | tR)| + |P(NC | tL) - P(NC | tR)|
    df_ret['Q(s|t)'] = abs(df_ret['C | t_L'] - df_ret['C | t_R']) + abs(df_ret['NC | t_L'] - df_ret['NC | t_R'])
    #F(s,t) = 2PLPR * Q(s|t)
    df_ret['F(s, t)'] = (2 * df_ret['P_L'] * df_ret['P_R']) * df_ret['Q(s|t)']


    #pC,t = n(t, C) / n(t) main
    df_ret['p_(C,t)'] = (df_ret['n(t_L, C)'] + df_ret['n(t_R, C)']) / n_t
    #pNC,t = n(t, NC) / n(t)  
    df_ret['p_(NC,t)'] = (df_ret['n(t_L, NC)'] + df_ret['n(t_R, NC)']) / n_t

    #pC,t = n(t, C) / n(t) left
    df_ret['p_(C,t_L)'] = (df_ret['n(t_L, C)']) / n_t
    #pNC,t = n(t, NC) / n(t)  
    df_ret['p_(NC,t_L)'] = (df_ret['n(t_L, NC)']) / n_t

    #pC,t = n(t, C) / n(t)  right
    df_ret['p_(C,t_R)'] = (df_ret['n(t_R, C)']) / n_t
    #pNC,t = n(t, NC) / n(t)  
    df_ret['p_(NC,t_R)'] = (df_ret['n(t_R, NC)']) / n_t

    #H(t) = -[pC,t log2(pC,t) + pNC,t log2(pNC,t)]
    df_ret['H(t)'] = -((df_ret['p_(C,t)'] * np.log2(df_ret['p_(C,t)'].astype(float))) + (df_ret['p_(NC,t)'] * np.log2(df_ret['p_(NC,t)'].astype(float))))

    #H(t) = -[pC,t log2(pC,t) + pNC,t log2(pNC,t)] left
    df_ret['H(t_L)'] = -((df_ret['p_(C,t_L)'] * np.log2(df_ret['p_(C,t_L)'].astype(float))) + (df_ret['p_(NC,t_L)'] * np.log2(df_ret['p_(NC,t_L)'].astype(float))))
    
    #H(t) = -[pC,t log2(pC,t) + pNC,t log2(pNC,t)] right
    df_ret['H(t_R)'] = -((df_ret['p_(C,t_R)'] * np.log2(df_ret['p_(C,t_R)'].astype(float))) + (df_ret['p_(NC,t_R)'] * np.log2(df_ret['p_(NC,t_R)'].astype(float))))

    #H(s,t) = PLH(tL) + PRH(tR)
    df_ret['H(s,t)'] = ((df_ret['P_L'] * df_ret['H(t_L)']) + (df_ret['P_R'] * df_ret['H(t_R)']))

    #gain(s) = H(t) - H(s,t)
    df_ret['gain(s)'] = df_ret['H(t)'] - df_ret['H(s,t)']

    return df_ret

def get_left_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 1].reset_index(drop=True)
    return df_ret

def get_right_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 0].reset_index(drop=True)
    return df_ret

def make_tree(head, df_training):
    df_left = get_left_df(head, df_training)
    #print(df_left)

    df_right = get_right_df(head, df_training)
    #print(df_right)

    df_left_with_stats = df_add_stats(df_left)
    #print(df_left_with_stats)
    df_right_with_stats = df_add_stats(df_right)
    #print(df_right_with_stats)

    #df_left_with_stats.drop(index = 0, inplace = True)
    #df_left_with_stats.reset_index(inplace = True, drop = True)
    
    df_left_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_left_with_stats.reset_index(inplace = True, drop = True)
    #print(df_left_with_stats)
    left_head = df_left_with_stats['Mutation List'][0]
    #print(left_head)
    
    df_right_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_right_with_stats.reset_index(inplace = True, drop = True)
    #print(df_right_with_stats)
    right_head = df_right_with_stats['Mutation List'][0]
    #print(right_head)

    df_left_left = get_left_df(left_head, df_left)
    #print(df_left_left)

    df_left_right = get_right_df(left_head, df_left)
    #print(df_left_right)

    df_right_left = get_left_df(right_head, df_right)
    #print(df_right_left)

    df_right_right = get_right_df(right_head, df_right)
    #print(df_right_right)

    df_left_left_with_stats = df_add_stats(df_left_left)
    df_left_left_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_left_left_with_stats.reset_index(inplace = True, drop = True)
    #print(df_left_left_with_stats.head(10))

    df_left_right_with_stats = df_add_stats(df_left_right)
    df_left_right_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_left_right_with_stats.reset_index(inplace = True, drop = True)
    #print(df_left_right_with_stats.head(10))

    df_right_left_with_stats = df_add_stats(df_right_left)
    df_right_left_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_right_left_with_stats.reset_index(inplace = True, drop = True)
    #print(df_right_left_with_stats.head(10))


    df_right_right_with_stats = df_add_stats(df_right_right)
    df_right_right_with_stats.sort_values('gain(s)', ascending = False, inplace = True)
    df_right_right_with_stats.reset_index(inplace = True, drop = True)
    #print(df_right_right_with_stats.head(10))

    h_list = [head, left_head, right_head, df_left_left_with_stats, df_left_right_with_stats, df_right_left_with_stats, df_right_right_with_stats, df_left_left, df_left_right, df_right_left, df_right_right]
    return h_list

def tree_output(tree, df_out_of_bag):
    temp_df1 = tree[3]
    temp_df2 = tree[4]
    temp_df3 = tree[5]
    temp_df4 = tree[6]

    print("Head = ", tree[1])
    # print("Out of bag set size = ", df_out_of_bag.shape[0])
    # print("Out of bag samples = ")
    # print(df_out_of_bag['Unnamed: 0'])
    # print("Left of the root node = ", tree[1])
    # print("Right of the root node = ", tree[2])
    # print("Left of the left child node = ", temp_df1['Mutation List'][0])
    # print("Right of the left child node = ", temp_df2['Mutation List'][1])
    # print("Left of the right child node = ", temp_df3['Mutation List'][1])
    # print("Right of the right child node = ", temp_df4['Mutation List'][0])


def make_single_tree(df_initial):
    df_bootstrap = df_initial.sample(n = df_initial.shape[0], replace = True)

    df_out_of_bag = df_initial.drop(df_bootstrap.index)

    unnamed0_series = df_bootstrap.iloc[:,0]

    df_bootstrap.drop('Unnamed: 0', inplace = True, axis = 1)

    df_bootstrap = df_bootstrap.sample(axis = 1, n = int(math.sqrt(df_initial.shape[1])))

    df_bootstrap.insert(0, "Unnamed: 0", unnamed0_series)

    df_bootstrap.reset_index(drop = True, inplace = True)

    #print(df_bootstrap)

    df_bootstrap_sorted = df_add_stats(df_bootstrap)

    df_bootstrap_sorted.sort_values('gain(s)', ascending = False, inplace = True)

    df_bootstrap_sorted.reset_index(drop = True, inplace = True)

    #print(df_bootstrap_sorted)

    head = df_bootstrap_sorted['Mutation List'][0]

    #print(head)

    tree = make_tree(head, df_bootstrap)

    return tree

def make_forest(df_initial, size):
    tree_list = list()
    i = 0

    while i < size:
        tree_temp = make_single_tree(df_initial)
        tree_list.append(tree_temp)
        i += 1
    
    return tree_list

def get_node_status(tree):
    node_status = list()
    i = 7

    while i < 11:
        df_node = tree[i]
        df_c = df_node[df_node["Unnamed: 0"].str.contains(pat = "NC") == False].reset_index(drop = True) 
        df_nc = df_node[df_node["Unnamed: 0"].str.contains(pat = "NC")].reset_index(drop = True) 
        if df_c.shape[0] > df_nc.shape[0]:
            node_status.append('cancer')
        else:
            node_status.append('not cancer')
        i += 1

    return node_status

def classifier(forest, size):
    node_status_list = list()
    i = 0

    while i < size:
        working_tree = forest[i]
        working_tree_status = get_node_status(working_tree)
        node_status_list.append(working_tree_status)
        i += 1

    return node_status_list


def forest_output(forest, size):
    important_mutations_list = list()
    i = 0

    while i < size:
        working_tree = forest[i]
        important_mutations_list.append(working_tree[0])
        important_mutations_list.append(working_tree[1])
        important_mutations_list.append(working_tree[2])
        i += 1
    
    return important_mutations_list

def process_sample(head, left_head, right_head, sample, df_test, node_list):
    if((df_test.loc[df_test['Unnamed: 0'] == sample, head].iloc[0] == 1) & (df_test.loc[df_test['Unnamed: 0'] == sample, left_head].iloc[0] == 1)):
        return node_list[0]
    elif((df_test.loc[df_test['Unnamed: 0'] == sample, head].iloc[0] == 1) & (df_test.loc[df_test['Unnamed: 0'] == sample, left_head].iloc[0] == 0)):
        return node_list[1]
    elif((df_test.loc[df_test['Unnamed: 0'] == sample, head].iloc[0] == 0) & (df_test.loc[df_test['Unnamed: 0'] == sample, right_head].iloc[0] == 1)):
        return node_list[2]
    else:
        return node_list[3]

def eval_sample(classifier_list, sample, forest, size, df_test):
    i = 0
    results_list = list()

    while i < size:
        working_tree = forest[i]
        working_node_list = classifier_list[i]
        result = process_sample(working_tree[0], working_tree[1], working_tree[2], sample, df_test, working_node_list)
        print("Sample ", sample, " according to tree ", i, " = ", result)
        results_list.append(result)
        i += 1

    return results_list

df_initial = pd.read_csv("/home/kilkenney1/datamining2/mutations(1).csv")

size = 10000
sample1 = 'C1'
sample2 = 'C10'
sample3 = 'C50'
sample4 = 'NC5'
sample5 = 'NC15'

forest = make_forest(df_initial, size)

classifier_list = classifier(forest, size)
#print(classifier_list[1])

#final_list = eval_sample(classifier_list, sample1, forest, size, df_initial)

# i = 0
# ccount = 0
# nccount = 0
# while i < size:
#     if final_list[i] == 'cancer':
#         ccount += 1
#     else:
#         nccount += 1
#     i += 1

# print("Total amount of cancer trees = ", ccount)
# print("Total amount of non cancer trees = ", nccount)

# if ccount > nccount:
#     print('Majority says cancer')
# else: 
#     print('Majority says non cancer')


temp_list = forest_output(forest, size)

important_muts_df = pd.DataFrame(temp_list)
    
print(important_muts_df.value_counts())