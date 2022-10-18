from cProfile import label
from pickle import FALSE
from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk 
import random

df_initial = pd.read_csv("/home/kilkenney1/datamining2/mutations(1).csv")

#feature table that shows top 10 mutations ranked by accuracy and show their values
#forumla -> accuracy = tp/(tp+tn)*100


def spiltdf(df_table):
    shuffled = df_table.sample(frac=1)
    df_ret = np.array_split(shuffled, 3)  
    return df_ret

df_split = spiltdf(df_initial)


df_split_sub1 = df_split[0]
df_split_sub2 = df_split[1]
df_split_sub3 = df_split[2]

#print(df_split_sub3)

df_training1 = pd.concat([df_split_sub1, df_split_sub2], axis=0, ignore_index=True)
df_training2 = pd.concat([df_split_sub1, df_split_sub3], axis=0, ignore_index=True)
df_training3 = pd.concat([df_split_sub2, df_split_sub3], axis=0, ignore_index=True)

df_test1 = df_split_sub1.reset_index(drop = True)
df_test2 = df_split_sub2.reset_index(drop = True)
df_test3 = df_split_sub3.reset_index(drop = True)

#print(df_test1)

#print(df_training)
    

def df_add_stats(df_table):
    col_sum = df_table.sum(axis = 0, numeric_only = True) #count up the number of times each mut appears and store it in a frame
    df_ret = col_sum.to_frame('Total')
    df_ret = df_ret.reset_index()
    df_ret.rename(columns={'index':'Mutation List'}, inplace=True)

    df_nc = df_table[df_table["Unnamed: 0"].str.contains(pat = "NC")].reset_index() #create a dataframe of just nc
    df_nc = df_nc.drop("index", axis = 1)

    df_c = df_table[df_table["Unnamed: 0"].str.contains(pat = "NC") == False].reset_index() #create a dataframe of just c
    df_c = df_c.drop("index", axis = 1)

    df_nc_sum = df_nc.sum(axis = 0).to_frame('False Positive') #count up the number of times each mut appears for nc and store it in a frame
    df_nc_sum = df_nc_sum.reset_index()
    df_nc_sum.rename(columns={'index':'Mutation List'}, inplace=True)

    df_c_sum = df_c.sum(axis = 0).to_frame('True Positive') #count up the number of times each mut appears for c and store it in a frame
    df_c_sum = df_c_sum.reset_index()
    df_c_sum.rename(columns={'index':'Mutation List'}, inplace=True)

    df_ret = df_ret.merge(df_nc_sum, on = "Mutation List") # merge all of the tables
    df_ret = df_ret.merge(df_c_sum, on = "Mutation List")

    df_ret['True Negative'] = df_nc.shape[0] - df_ret['False Positive'] #find tn and fn and then return
    df_ret['False Negative'] = df_c.shape[0] - df_ret['True Positive'] 
    df_ret['TP-FP'] = df_ret['True Positive'] - df_ret['False Positive']

    df_ret['Accuracy'] = (df_ret['True Positive'] / (df_ret['True Positive'] + df_ret['True Negative'] + df_ret['False Positive'] + df_ret['False Negative']))
    df_ret = df_ret.sort_values('TP-FP', ascending = False).reset_index(drop=True) 

    return df_ret

#df_with_all_stats = df_add_stats(df_initial)
df_with_all_stats1 = df_add_stats(df_training1)
df_with_all_stats2 = df_add_stats(df_training2)
df_with_all_stats3 = df_add_stats(df_training3)
#print(df_with_all_stats[['Mutation List','Accuracy']])


#print(df_with_all_stats.head(10))

head1 = df_with_all_stats1['Mutation List'][0]
head2 = df_with_all_stats2['Mutation List'][0]
head3 = df_with_all_stats3['Mutation List'][0]

#print(head1)


def get_left_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 1].reset_index(drop=True)
    return df_ret

def get_right_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 0].reset_index(drop=True)
    return df_ret

#make this return a list of right and left head, replace df_initial with the training set
def make_tree(head, df_training):
    df_left = get_left_df(head, df_training)
    #print(df_left)

    df_right = get_right_df(head, df_training)
    #print(df_right)

    df_left_with_stats = df_add_stats(df_left)
    #print(df_left_with_stats)
    df_right_with_stats = df_add_stats(df_right)
    #print(df_right_with_stats)

    df_left_with_stats.drop(index = 0, inplace = True)
    df_left_with_stats.reset_index(inplace = True, drop = True)
    #print(df_left_with_stats)
    left_head = df_left_with_stats['Mutation List'][0]
    #print(left_head)

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
    #print(df_left_left_with_stats.head(10))
    df_left_right_with_stats = df_add_stats(df_left_right)
    #print(df_left_right_with_stats.head(10))

    df_right_left_with_stats = df_add_stats(df_right_left)
    #print(df_right_left_with_stats.head(10))
    df_right_right_with_stats = df_add_stats(df_right_right)
    #print(df_right_right_with_stats.head(10))

    h_list = [head, left_head, right_head]
    return h_list

tree1 = make_tree(head1, df_training1)
tree2 = make_tree(head2, df_training2)
tree3 = make_tree(head3, df_training3)

#def eval_sample(sample, head1, left_head, right_head, test_sample):
    #if((sample.loc[sample['Unnamed: 0'] == test_sample, head1].iloc[0] == 1) & (sample.loc[sample['Unnamed: 0'] == test_sample, left_head].iloc[0] == 1)):
        #return 'cancer'
    #elif((sample.loc[sample['Unnamed: 0'] == test_sample, head1].iloc[0] == 1) & (sample[sample['Unnamed: 0'] == test_sample, left_head].iloc[0] == 0)):
        #return 'not cancer'
   #elif((sample.loc[sample['Unnamed: 0'] == test_sample, head1].iloc[0] == 0) & (sample.loc[sample['Unnamed: 0'] == test_sample, right_head].iloc[0] == 1)):
        #return 'cancer'
    #else:
        #return 'not cancer'

def eval_sample(head, left_head, right_head, sample):
    if((df_initial.loc[df_initial['Unnamed: 0'] == sample, head].iloc[0] == 1) & (df_initial.loc[df_initial['Unnamed: 0'] == sample, left_head].iloc[0] == 1)):
        return 'cancer'
    elif((df_initial.loc[df_initial['Unnamed: 0'] == sample, head].iloc[0] == 1) & (df_initial.loc[df_initial['Unnamed: 0'] == sample, left_head].iloc[0] == 0)):
        return 'not cancer'
    elif((df_initial.loc[df_initial['Unnamed: 0'] == sample, head].iloc[0] == 0) & (df_initial.loc[df_initial['Unnamed: 0'] == sample, right_head].iloc[0] == 1)):
        return 'cancer'
    else:
        return 'not cancer'

def eval_out(tree, df_training, df_test):
    counter = 0
    true_pos = 0
    true_neg = 0
    false_pos = 0
    false_neg = 0

    while counter < df_test.shape[0]:
        test_sample = df_test['Unnamed: 0'][counter]
        sample_eval = eval_sample(tree[0], tree[1], tree[2], test_sample)
        if((sample_eval == 'cancer') & ('NC' not in test_sample)):
            true_pos += 1
        elif((sample_eval == 'cancer') & ('NC' in test_sample)):
            false_pos += 1
        elif((sample_eval == 'not cancer') & ('NC' not in test_sample)):
            false_neg += 1
        else:
            true_neg += 1
        counter += 1

    print("Running the file through the decision tree shows that")
    print("True Positives = ", true_pos)
    print("False Positives = ", false_pos)
    print("False Negatives = ", false_neg)
    print("True Negatives = ", true_neg)
    df_con = pd.DataFrame(columns = ['POSITIVE', 'NEGATIVE'], index = ['POSITIVE', 'NEGATIVE'])
    df_con['POSITIVE']['POSITIVE'] = true_pos
    df_con['POSITIVE']['NEGATIVE'] = false_pos
    df_con['NEGATIVE']['POSITIVE'] = false_neg
    df_con['NEGATIVE']['NEGATIVE'] = true_neg
    print(df_con)
    print("From these numbers, we can determine these statistics")
    print("Accuracy = ", ((true_pos + true_neg) / (true_pos + true_neg + false_pos + false_neg)) * 100)
    acc = (true_pos + true_neg) / (true_pos + true_neg + false_pos + false_neg) * 100
    print("Sensitivity = ", (true_pos / (true_pos + false_neg)) * 100)
    print("Specificity = ", (true_neg / (true_neg + false_pos)) * 100)
    print("Precision = ", (true_pos / (true_pos + false_pos)) * 100)
    print("Miss rate = ", (false_neg / (false_neg + true_pos)) * 100)
    print("False discovery rate = ", (false_pos / (false_pos + true_pos)) * 100)
    print("False omission rate = ", (false_neg / (false_neg + true_neg)) * 100)
    return acc

acc1 = eval_out(tree1, df_training1, df_test3)
acc2 = eval_out(tree2, df_training2, df_test2)
acc3 = eval_out(tree3, df_training3, df_test1)

avgacc = (acc1 + acc2 + acc3)/3
print("AVG ACCURACY = ", avgacc)

