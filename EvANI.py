#!/usr/bin/env python3


"""


How to run


python3 EvANI.py mytool  mutation 

The last argument could be one of the following mutation duplication lgt
mytool includes tsv files of the ani values. The format should be this
SE001.fa SE002.fa 98

folder structure:
Check the github folder precomputed_distances including mash and dashing results

"""



import scipy
from ete3 import Tree
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math





def read_matrix_dashing(filename):
    sample_names = []
    dic_distance = {}
    file_handle =  open(filename,'r')
    line_counter = -2
    for line in file_handle:
        line_counter += 1
        if line.startswith("#"): # \t  # skip sample name line
            sample_names_raw = line.strip().split('\t')
            sample_names = [i.split("/")[-1] for i in sample_names_raw] # split(".")[0] 
            #print("number of samples ",len(sample_names))
            continue        
        elements = line.strip().split('\t')
        for i in range(line_counter+1,len(sample_names)):
            s_min=min(sample_names[line_counter],sample_names[i])
            s_max=max(sample_names[line_counter],sample_names[i])
            dic_distance[(s_min,s_max)]=float(elements[i+1])
            #dic_distance[(s_max,s_min)]=float(elements[i+1])
    #print(elements)
    return dic_distance


def read_matrix_mash(filename):
    sample_names = []
    dic_distance = {}
    file_handle =  open(filename,'r')
    line_counter = -2
    for line in file_handle:
        line_counter += 1
        if line.startswith("#"): # \t  # skip sample name line
            sample_names_raw = line.strip().split('\t')
            sample_names = [i.split("/")[-1] for i in sample_names_raw] # .split(".")[0]
            sample_names= sample_names[1:]
            #print("number of samples ",len(sample_names))
            continue        
        elements = line.strip().split('\t')
        for i in range(line_counter+1,len(sample_names)):
            s_min=min(sample_names[line_counter],sample_names[i])
            s_max=max(sample_names[line_counter],sample_names[i])
            dic_distance[(s_min,s_max)]=float(elements[i+1])
    #print(elements)
    return dic_distance



def dist_tree_calc(filename,extension=".fa"): # , quoted_node_names=False)
    tree1= Tree(filename,format=1) # , quoted_node_names=quoted_node_names
    #print(filename, len(tree1))    
    dist_tree={}
    dist_tree_top={}
    samples = [i.name for i in tree1.get_leaves()]
    num_sample = len(samples)
    for i in range(num_sample):
        tax_i= samples[i]
        for j in range(i+1,num_sample):
            tax_j= samples[j]
            dist_toplg = tree1.get_distance(tax_i,tax_j, topology_only=True) 
            dist_ = tree1.get_distance(tax_i,tax_j, topology_only=False) 
            dist_tree[(min(tax_i,tax_j)+extension,max(tax_i,tax_j)+extension)]=dist_
            dist_tree_top[(min(tax_i,tax_j)+extension,max(tax_i,tax_j)+extension)]=dist_toplg
    #print(num_sample,num_sample*(num_sample-1)/2, len(dist_tree_top))
    return dist_tree


def read_fastani(filename):
    dic_ani = {}
    file_handle =  open(filename,'r')
    for line in file_handle:            
        elements = line.strip().split('\t')
        sp1= elements[0].split("/")[-1]#.split(".")[0]
        sp2= elements[1].split("/")[-1]#.split(".")[0]
        ani= float(elements[2]) /100
        if sp1!=sp2:
            s_min=min(sp1,sp2)
            s_max=max(sp1,sp2)            
            dic_ani[(s_min,s_max)]=ani
    #print(elements)
    return dic_ani







if __name__ == "__main__":

    folder_sample=sys.argv[1] #"precomputed_distances/sample"
    study=sys.argv[2] # study, either mutation duplication lgt

    rates_dic={'mutation':['5','10','25','50','100'],'lgt':["0.0001", "0.0005", "0.0010", "0.0020"],
             'duplication':["0.0000","0.0005", "0.0010", "0.0020"]}
    rates=rates_dic[study]
    replicates=list(range(1,6))



    ###### read  folder  of the tool under test including tsv files
    dic_sample={}
    for rate in rates:
        for replicate in replicates:
            file_name=folder_sample+"/"+study+"_"+str(rate)+"_"+str(replicate)+".tsv"
            dic_rate_replicate_sample = read_fastani(file_name)
            dic_sample[str(rate)+"_"+str(replicate)] = dic_rate_replicate_sample
    print("For input tool we parse  "+str(len(dic_sample))+" cases. For case "+str(rate)+"_"+str(replicate)+", we parsed "+str(len(dic_sample[str(rate)+"_"+str(replicate)]))+" species pairs"   )


    ###### read  folder  of the dashing and mash
    dic_dashing={}
    dic_mash={}
    for rate in rates:
        for replicate in replicates:
            tool="dashing-k21"
            file_name="precomputed_distances/"+tool+"/"+tool+"_"+study+"_"+str(rate)+"_"+str(replicate)+".tsv"
            dic_rate_replicate = read_matrix_dashing(file_name)
            dic_dashing[str(rate)+"_"+str(replicate)]= dic_rate_replicate

            tool="mash-k21-s10k"
            file_name="precomputed_distances/"+tool+"/"+tool+"_"+study+"_"+str(rate)+"_"+str(replicate)+".tsv"
            dic_rate_replicate = read_matrix_mash(file_name)
            dic_mash[str(rate)+"_"+str(replicate)]= dic_rate_replicate
    print("Mash and dashing distances are paresed")


    ###### read trees using ete3
    trees={}
    if study=="mutation":
        for rate in rates:
            filename="trees/simulated_tree_mutation_"+str(rate)+".nwk"
            tree_distances=dist_tree_calc(filename,extension=".fa")
            trees[study+"_"+str(rate)]=tree_distances

    if study=="duplication" or study=="lgt" :
        filename="trees/simulated_tree_duplication_lgt.nwk"
        tree_distances=dist_tree_calc(filename,extension=".fa")
        trees[study]=tree_distances
    print("Trees were parsed")

    
    ###### create a dataframe of pvalue of spearman correlation test
    df_stat = pd.DataFrame({'study':[],'rate':[],'tool':[], 'logpval': [], 'statistics': [] })
    if study=="duplication" or study=="lgt" :
        tree_distances = trees[study]

    for rate in rates:
        if study=="mutation":
            tree_distances = trees[study+"_"+str(rate)]

        species_pair_list = list(tree_distances.keys())
        tree_list =[tree_distances[pair] for pair in species_pair_list]

        for replicate in replicates:

            tool="dashing-k21"
            dic_rate_replicate=dic_dashing[str(rate)+"_"+str(replicate)]
            ani_list =[dic_rate_replicate[pair] for pair in species_pair_list]
            result_spearman=scipy.stats.spearmanr(ani_list,tree_list)
            logpval = math.log10(result_spearman.pvalue)
            statistics = result_spearman.statistic
            df_stat = pd.concat([pd.DataFrame([[study,float(rate),tool,logpval,statistics]], columns=df_stat.columns), df_stat], ignore_index=True)

            tool="mash-k21-s10k"
            dic_rate_replicate=dic_mash[str(rate)+"_"+str(replicate)]
            ani_list =[1- dic_rate_replicate[pair] for pair in species_pair_list]
            result_spearman=scipy.stats.spearmanr(ani_list,tree_list)
            logpval = math.log10(result_spearman.pvalue)
            statistics = result_spearman.statistic
            df_stat = pd.concat([pd.DataFrame([[study,float(rate),tool,logpval,statistics]], columns=df_stat.columns), df_stat], ignore_index=True)

            tool=folder_sample
            dic_rate_replicate=dic_sample[str(rate)+"_"+str(replicate)]
            ani_list=[]
            for pair in species_pair_list:
                if pair in dic_rate_replicate:
                    ani_list.append(dic_rate_replicate[pair])
                else:
                    ani_list.append(0)

            result_spearman=scipy.stats.spearmanr(ani_list,tree_list)
            logpval = math.log10(result_spearman.pvalue)
            statistics = result_spearman.statistic
            df_stat = pd.concat([pd.DataFrame([[study,float(rate),tool,logpval,statistics]], columns=df_stat.columns), df_stat], ignore_index=True)

    ###### creat plots and save them as pdf file

    for statistics_logpval in ['logpval','statistics']:
        plt.figure()
        sns.pointplot(data=df_stat, x='rate', y=statistics_logpval, hue='tool', native_scale=True ,errorbar=('pi', 50))
        plt.savefig("EvANI_output_"+folder_sample.split("/")[-1]+"_"+study+"_"+statistics_logpval+".pdf")

        print("Figure is stored in " + "EvANI_output_"+folder_sample.split("/")[-1]+"_"+study+"_"+statistics_logpval+".pdf" )



