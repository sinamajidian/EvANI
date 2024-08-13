


# use the functions in ani_parser.py


root_folder="realdata/"
dataset_="o__Bacillales_A2" 

fasta_folder = root_folder+dataset_+"/cds/" # genomes_ncbi fastoma_100
extension= ""
files=os.listdir(fasta_folder)
fastas= [i for i in files if i.endswith(extension)]




tree2_f=root_folder+dataset_+"/GTDB_249o__Bacillales_A.nwk"
dist_tree= dist_tree_calc(tree2_f,extension="")

species_pair_list = list(dist_tree.keys())
print(species_pair_list[:2], len(species_pair_list))

species_pair_list_used = species_pair_list_
dic_dashing ={}
ext=".fna"
df_pval = pd.DataFrame({'k': [] , 'logpval': [],'pval':[], 'sta':[] ,"dataset":[]}) 
extension= ".fna"
dist_tree_list = [dist_tree[pair] for pair in species_pair_list_used]

dic_dashing ={}
for k in range(5,33): 
    #print(k)

    filename1 = root_folder+dataset_+"/dashingfull_cds100/"+"/distance_k"+str(k)+".txt" # all_withplasmid
    dic_dashing["k"+str(k)] = read_matrix_dashing(filename1) 
    ani_list =[dic_dashing["k"+str(k)][pair] for pair in species_pair_list_used]     # (pair[0]+ext,pair[1]+ext)
    
    res=scipy.stats.spearmanr(dist_tree_list, ani_list)
    
    logpval = math.log10(res.pvalue)  # res.pvalue  # math.log10(res.pvalue) #math.log10(res.pvalue) # res.statistic#res.statistic#math.log10(res.pvalue) #res.statistic # sta#math.log10(res.pvalue)#math.log10(res.pvalue) #res.statistic #res.pvalue #math.log10(res.pvalue)
    df_pval = pd.concat([pd.DataFrame([[k,logpval,res.pvalue, res.statistic ,"100 random cds genes"]], columns=df_pval.columns), df_pval], ignore_index=True)
    print("k"+str(k), round(res.statistic,3), round(logpval,1),res.pvalue)




fig,ax = plt.subplots(figsize=(9,6))


g=sns.lineplot(data=df_pval, x='k', y='logpval', marker='o' ,hue='dataset') #    logpval  sta , hue='size'sta
plt.ylabel("Statistic", fontsize=18) # Log P-value Statistic | whole genome with Plasmids | no plasmid  only species GC<0.45 (1580 pairs) #| species pairs from all except Nostocales/Synechococcales

ax.set_xticks(range(5,33))
plt.xlabel("k", fontsize=18) 
plt.title("Spearman: Dashing-FullKhash vs GTDB Distance tree  |"+dataset_) # #  SynechococcalesDistance tree | GC < 0.48 38 out of 81 species (CDS >4.5Mbp) 
# genome NCBI (67 out of 81 species with GC >0.38) #  | removing outlier of uniq-11mers (72 out of 81)
plt.show() 
# genomes_ncbi



