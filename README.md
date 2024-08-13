# EvANI


## EvANI benchmarking pipeline

This repo include scripts for benchmarking ANI tools.




## EvANI benchmarking dataset

We beneftited from ALF simulator and ran locally using https://github.com/DessimozLab/ALF. Small dataset can be also generated online [here](http://alf.cs.ucl.ac.uk/ALF/).  

A diverse range of datasets are provided on [Zenodo](https://zenodo.org/records/13308784)


For real data we used the GTDB tree
https://data.gtdb.ecogenomic.org/releases/release202/202.0/ 

We downloaded the genomes using esearch

```
wget  `esearch -db assembly -query ${i} | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`  -O ${i}.fna.gz 
```

