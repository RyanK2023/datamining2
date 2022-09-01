from pickle import FALSE
import pandas as pd
import matplotlib as mpl
import numpy as np

df = pd.read_csv("/mnt/c/Users/Ryan/Desktop/mutations(1).csv")
cols = df.shape[1] - 1 #number of unique mutations = 3816
rows = df.shape[0] - 1 #number of individual samples = 229

df_c1 = df.iloc[0].value_counts() #shows that c1 has ten mutations 
df_nc1 = df.iloc[3].value_counts() #shows that c1 has seven mutations 

df_drop = df.drop("Unnamed: 0", axis = 1)

df_mean = df_drop.values.sum()/rows #says that a person has an average of 43.46 mutations

#print(cols)
#print(rows)

df.loc[:,'row count'] = df.sum(numeric_only=True, axis=1)

#print(df[df['row count']==df['row count'].max()]) #says that c8 has a max of 415 mutations 

#rint(df[df['row count']==df['row count'].min()]) #says that nc110 has a min of 0 mutations 

df = df.drop("row count", axis = 1)

col_sum = df_drop.sum(axis = 0)

#print(col_sum.max()) #says that one mutation has a max of 42 samples

#print(col_sum.min()) #says many different mutations have a min of 2 samples 

df_mean_cols = df_drop.values.sum()/cols

#print(df_mean_cols)  #says that a specific mutation on average is present in 2.608 of the samples 


df_sub = df.columns[df.columns.str.contains(pat = "BRAF")] #says that 34 people have the mutation 

#print(df[df_sub].sum())

df_sub2 = df.columns[df.columns.str.contains(pat = "KRAS")]
#print(df[df_sub2].drop_duplicates().sum()) #says that 13 people have some variation of KRAS

#df['NumMutPerSample'] = df.sum(numeric_only=True, axis=1)

#df.plot.scatter(x = 'Unnamed: 0', y = 'NumMutPerSample')
#add a line for second scatter plot, will figure something out later 
#print(col_sum)
df_table = col_sum.to_frame('Total count')
df_table = df_table.reset_index()
df_table.rename(columns={'index':'Mutation List'}, inplace=True)


df_nc = df[df["Unnamed: 0"].str.contains(pat = "NC")].reset_index()
df_nc = df_nc.drop("index", axis = 1)
#print(df_nc)

df_c = df[df["Unnamed: 0"].str.contains(pat = "NC") == False].reset_index()
df_c = df_c.drop("index", axis = 1)
#print(df_c)

df_nc_sum = df_nc.sum(axis = 0).to_frame('Total nc')
df_nc_sum = df_nc_sum.reset_index()
df_nc_sum.rename(columns={'index':'Mutation List'}, inplace=True)

df_c_sum = df_c.sum(axis = 0).to_frame('Total c')
df_c_sum = df_c_sum.reset_index()
df_c_sum.rename(columns={'index':'Mutation List'}, inplace=True)


#this feels wrong but it works 
df_table = df_table.merge(df_nc_sum, on = "Mutation List")
df_table = df_table.merge(df_c_sum, on = "Mutation List")

df_table['%C'] = df_table['Total c']/110
df_table['%NC'] = df_table['Total nc']/120
df_table['%C-%NC'] = df_table['%C'] - df_table['%NC']
df_table['%C/%NC'] = df_table['%C'] / df_table['%NC']

print(df_table.head(10).sort_values('Total count', ascending = False))
print(df_table.head(10).sort_values('Total c', ascending = False))
print(df_table.head(10).sort_values('Total nc', ascending = False))
print(df_table.head(10).sort_values('%C', ascending = False))
print(df_table.head(10).sort_values('%NC', ascending = False))
print(df_table.head(10).sort_values('%C-%NC', ascending = False))
print(df_table.head(10).sort_values('%C/%NC', ascending = False))
