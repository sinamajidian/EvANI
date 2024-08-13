


### fastANI

extension="fa"
out="fastani/"
mkdir $out; threads=35
ls fasta/*${extension}| wc -l
ls -1a fasta/*${extension} > $out/fasta_list
ls $out
fastANI --ql $out/fasta_list --rl $out/fasta_list -o $out/distances.tab --threads $threads --minFraction 0.1 --fragLen 2000




### NUCmer



out=nucmer
mkdir ${out}

for r in $(seq -w 1 15); do 
for l in 8 12 14 16 18 20 24 ; do 
(
for q in $(seq -w 1 15); do 

ref="SE0"$r.fa
query="SE0"$q.fa

nucmer -l ${l} --mum -p $out/${ref}_${query}_l${l} fasta/${ref} fasta/${query}
delta-filter -1 $out/${ref}_${query}_l${l}.delta > $out/${ref}_${query}_l${l}.filter 
done

) & 
done
done
wait




## orthoANI 

oat="OAT_cmd.jar"
blastbin="blast/ncbi-blast-2.15.0+/bin/"
for r in $(seq -w 1 15); do 
for q in $(seq -w 1 15); do 
ref="SE0"$r.fa
query="SE0"$q.fa
java -jar ${oat} -blastplus_dir ${blastbin}  -num_threads  15 -fasta1 fasta/${query} -fasta2 fasta/${ref} > ortho/${ref}_${query}.out

done
done