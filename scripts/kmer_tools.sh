


#dashing

threads=10
size=14
k=21

dashing dist  -k $k -p $threads -O dashing_out_s14/distance.txt -o dashing_out_s14/sizes.txt -S $size --full-mash-dist --sketch-size $size --full-tsv fasta/*




# kmc

for file in $(ls -1 ../genomes_ncbi/); do 

for k in {5..33}; do
kmc  -k${k} -t39 -fm  ../genomes_ncbi/$file out_${i} out_kmc > ${file}_${k}.stat.kmc

done
done



# mash

mkdir mash_out
threads=10
size=10000
k=21
mash sketch -p ${threads} -o mash_out/ref -s ${size} -k $k fasta/* 
mash dist -p ${threads} mash_out/ref.msh mash_out/ref.msh -t > mash_out/distances.tab

