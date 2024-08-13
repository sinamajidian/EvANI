


root_folder="/vast/blangme2/smajidi5/lca/real/"

dataset_="o__Bacillales_A2" # f__Neisseriaceae o__Bacillales_A
# tree_file={"Cyanobacteria":"NCBI_84Cyanobacteria.nwk", "Chlamydiia":"NCBI_101Chlamydiia.nwk",
#            "Synechococcales":"NCBI_32Synechococcales.nwk"}

fasta_folder = root_folder+dataset_+"/cds/" # cds_ncbi  genomes_ncbi bac120
extension= "" #".fna"
files=os.listdir(fasta_folder)
fastas= [i for i in files if i.endswith(extension)]
len(fastas)


