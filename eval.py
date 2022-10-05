from cProfile import label
from pickle import FALSE
from turtle import color
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk 

df_initial = pd.read_csv("/mnt/c/Users/Ryan/Desktop/mutations(1).csv")

#feature table that shows top 10 mutations ranked by accuracy and show their values
#forumla -> accuracy = tp/(tp+tn)*100

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

df_with_all_stats = df_add_stats(df_initial)
#print(df_with_all_stats[['Mutation List','Accuracy']])


#print(df_with_all_stats.head(10))

head = df_with_all_stats['Mutation List'][0]

#print(head)


def get_left_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 1].reset_index(drop=True)
    return df_ret

def get_right_df(mut, df_table):
    df_ret = df_table.loc[df_table[mut] == 0].reset_index(drop=True)
    return df_ret

df_left = get_left_df(head, df_initial)
#print(df_left)

df_right = get_right_df(head, df_initial)
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

def eval_sample(sample, head, left_head, right_head):
    if((df_initial.loc[df_initial['Unnamed: 0'] == test_sample, head].iloc[0] == 1) & (df_initial.loc[df_initial['Unnamed: 0'] == test_sample, left_head].iloc[0] == 1)):
        return 'cancer'
    elif((df_initial.loc[df_initial['Unnamed: 0'] == test_sample, head].iloc[0] == 1) & (df_initial.loc[df_initial['Unnamed: 0'] == test_sample, left_head].iloc[0] == 0)):
        return 'not cancer'
    elif((df_initial.loc[df_initial['Unnamed: 0'] == test_sample, head].iloc[0] == 0) & (df_initial.loc[df_initial['Unnamed: 0'] == test_sample, right_head].iloc[0] == 1)):
        return 'cancer'
    else:
        return 'not cancer'

counter = 0
true_pos = 0
true_neg = 0
false_pos = 0
false_neg = 0
while counter < 230:
    test_sample = df_initial['Unnamed: 0'][counter]
    sample_eval = eval_sample(test_sample, head, left_head, right_head)
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
print("Sensitivity = ", (true_pos / (true_pos + false_neg)) * 100)
print("Specificity = ", (true_neg / (true_neg + false_pos)) * 100)
print("Precision = ", (true_pos / (true_pos + false_pos)) * 100)
print("Miss rate = ", (false_neg / (false_neg + true_pos)) * 100)
print("False discovery rate = ", (false_pos / (false_pos + true_pos)) * 100)
print("False omission rate = ", (false_neg / (false_neg + true_neg)) * 100)
#just gotta get this on to my laptop
