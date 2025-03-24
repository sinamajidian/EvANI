
This folder contains information on how to downloaded the genomes of real data from NCBI. 
We provided the how to use the `esearch`. The fast alternative is to use datasets package after installing with `pip install datasets`

```
$ datasets download genome accession --inputfile acc.txt --dehydrated
$ unzip ncbi_dataset.zip
$ datasets rehydrate --directory .
```
