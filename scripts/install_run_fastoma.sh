

conda create -n f39p python=3.9
conda activate f39p
git clone https://github.com/DessimozLab/FastOMA.git
conda install bioconda::fasttree

conda install bioconda::mafft

conda install conda-forge::openjdk


pip install .[report,nextflow] 




# run fastOMA


nextflow run FastOMA.nf --input_folder testdata/in_folder --output_folder out_folder --omamer_db LUCA.h5