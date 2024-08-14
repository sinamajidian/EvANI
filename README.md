# EvANI


## EvANI benchmarking pipeline

This repo includes scripts for benchmarking ANI tools. We provided a [bash script] to run each tool (including k-mer-based: dashing, mash, [here](https://github.com/sinamajidian/EvANI/blob/main/scripts/kmer_tools.sh), fastANI, orthoANI,  ANIm [here](https://github.com/sinamajidian/EvANI/blob/main/scripts/ANI_tools.sh) on genomes in Fasta fromat. 

The output of each tool is paresed in python using this [piece](https://github.com/sinamajidian/EvANI/blob/main/scripts/ani_parser.py). We compared the genomic distance found by each tool with the true distance on the phylogenetic tree using [this code](https://github.com/sinamajidian/EvANI/blob/main/scripts/rank_correlation_test.py) and visualized the results using seabron.


For inferring orthologous genes we used [FastOMA](https://github.com/DessimozLab/fastoma).



## EvANI benchmarking dataset


### Simualted data

We beneftited from [ALF simulator]((https://github.com/DessimozLab/ALF)) and ran locally using [the script](https://github.com/sinamajidian/EvANI/blob/main/scripts/run_ALF_locally.sh), based on the [parameter files](https://github.com/sinamajidian/EvANI/blob/main/scripts/ALF_sim-params.drw). Note that small dataset can be also generated online [here](http://alf.cs.ucl.ac.uk/ALF/).  

We varied parameters  mutation rate (`mutRate` from 2 to 20, varying ANI values 50-100), gene duplication rate (`geneDuplRate` from 0.0001 to 0.01), and rate of horizontal gene transfer (`lgtGRate` from 0.00001 to 0.001). A diverse range of datasets are provided on [EvANI Zenodo](https://zenodo.org/records/13308784).


### Real data

We downloaded the NCBI genomes using esearch

```
wget  `esearch -db assembly -query ${i} | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`  -O ${i}.fna.gz 
```

For real data we used the GTDB tree, available [here](https://data.gtdb.ecogenomic.org/releases/release202/202.0/).  For each accession ID, 

We also used [ete3](https://github.com/etetoolkit/ete/tree/3.0) to download the [NCBI taxanomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi) in Python.
```
import ete3
ncbi = ete3.NCBITaxa() 
ncbi_sub_tree = ncbi.get_topology(ncbi_taxon_list)
```



