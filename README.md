# EvANI


This repository includes the instructions how to run EvANI benchmarking pipeline on the benchmarking datasets. 

## Preprint
```
S. Majidian, , S. Hwang, M. Zakeri, & B. Langmead, (2025)
EvANI benchmarking workflow for evolutionary distance estimation
https://www.biorxiv.org/content/10.1101/2025.02.23.639716
```


### requirements
EvANI framework is a basic python script that performs a Spearman rank correlation test. 

```
conda create -n evani python=3.12
conda activate evani

conda install conda-forge::ete3
conda install conda-forge::seaborn
conda install conda-forge::matplotlib
```
The package ete3 is used to parse phylogenetic trees and calculate tree distances. 



## EvANI benchmarking pipeline


### Step 1: download simulated dataset

First, go to our [zenodo page](https://zenodo.org/records/14579845) and download the simulated dataset


The dataset includes three evolutionary scenarios: 

1. varying mutation rates including 5,10,25,50,100,
2. varying duplication rates including 0.0000,0.0005, 0.0010, 0.0020, 
3. varying rates of lateral gene transfer (LGT), also known as horizontal gene transfer (HGT) including 0.0001, 0.0005, 0.0010, 0.0020. Note that for lgt=0 you could use duplication=0 dataset.

For each stud/rate, there are five replicates of evolution simulation. Each case contain 15 genomes (DNA fasta files).


### Step 2: run your ANI tool 

Run your ANI tool on the simulated datasets and report the results in TSV files. The format of the file name is `study_rate_replicate.tsv`  e.g. `duplication_0.0005_2.tsv`. 
We provided a folder `sample_tool` including the outputs of a sample tool. 

(if you just want to try the benchmarking pipleline on the provided test data, see step 3).



### Step 3: run EvANI benchmarking  (test example)
First, clone this github repo:

```
git clone git@github.com:sinamajidian/EvANI.git

```

This provides you with the python script and the precomputed simulated dataset. Make sure you have installed the requirements.


Make sure you are in the correct folder where you see at least the folowing file and folders
```
cd EvANI
$ ls 
EvANI.py  precomputed_distances sample_tool  trees
```


The python code `EvANI.py` needs two positional arguments: the folder name `sample_tool` where the TSV files of ANI values are stored, and one of the studies `duplication`, `mutation` or `lgt`. Now run it as  


```
python EvANI.py sample_tool mutation

```


For a test of the script, you can run the same command line  `python EvANI.py sample_tool mutation`  using files in this repo without running your tools. 


This will output two figures in PDF, one for log p-values and one for statistics versus the rates.  For the sample_tool, the output will be 


<div align="center">
  <img width="300px" src="./sample_tool/expected_output/EvANI_output_sample_tool_mutation_logpval.jpg" alt="EvANI output figure" />
</div>


Sample tool here is [FastANI](https://github.com/ParBLiSS/FastANI) with Min Fraction of genome shared =0.1 and fragment length =3000.





## Note on EvANI benchmarking datasets


For generating simulated data, we benefited from [ALF simulator]((https://github.com/DessimozLab/ALF)) and ran locally using [the script](https://github.com/sinamajidian/EvANI/blob/main/scripts/run_ALF_locally.sh), based on the [parameter files](https://github.com/sinamajidian/EvANI/blob/main/scripts/ALF_sim-params.drw). Note that small dataset can also be generated online [here](http://alf.cs.ucl.ac.uk/ALF/).   Simulated datasets are freely available on our [zenodo page](https://zenodo.org/records/14579845).


The bash script to download the real genomes are provided in the folder `real_data`. For example for the Caldisericia clade, we have  a bash script and the phylogeny in newick format [here](https://github.com/sinamajidian/EvANI/tree/main/real_data/c__Caldisericia)
The bash script includes command line to download the NCBI genomes using esearch

```
wget  `esearch -db assembly -query ${i} | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`  -O ${i}.fna.gz 
```
A fast alternative is to use `datasets` package after installing with pip install datasets or `conda install conda-forge::ncbi-datasets-cli`.
```
echo "GCA_017999835.1" > acc.txt
datasets download genome accession --inputfile acc.txt --dehydrated
unzip ncbi_dataset.zip
datasets rehydrate --directory .
```

For real data we used the GTDB tree, available [here](https://data.gtdb.ecogenomic.org/releases/release202/202.0/). 

We also used [ete3](https://github.com/etetoolkit/ete/tree/3.0) to download the [NCBI taxanomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) in Python.
```
import ete3
ncbi = ete3.NCBITaxa() 
ncbi_sub_tree = ncbi.get_topology(ncbi_taxon_list)
```







