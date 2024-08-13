

# python functions for reading output of tools



def get_sgc(folder,extension=".fa",full=False):# -> Dict[str, int]:
    fastafilenames= [i for i in os.listdir(folder) if i.endswith(extension)]
    gc_dic = {} # from pyani package   
    for fname in fastafilenames:
        str1=""
        ff=SeqIO.parse(folder+fname, "fasta")
        for rec in ff:
            str1 += rec.seq
        gc= Bio.SeqUtils.gc_fraction(str1)
        if full:
            gc_dic[fname] = gc
        else:
            gc_dic[fname.split(".")[0]] = gc
        
    return gc_dic
    
def get_sequence_lengths(folder,fastafilenames,full=False):# -> Dict[str, int]:
      
    tot_lengths = {} # from pyani package   
    if full:
        for fname in fastafilenames:
            tot_lengths[fname] = sum([len(s) for s in SeqIO.parse(folder+fname, "fasta")])
    else:
        for fname in fastafilenames:
            tot_lengths[fname.split(".")[0]] = sum([len(s) for s in SeqIO.parse(folder+fname, "fasta")])
    return tot_lengths



def parse_mums2(filename):
    file_handle= open(filename,"r")
    len_mums= []
    for line in file_handle:
        if line.startswith(">"):
            if "Reverse" in line:
                print("error")
        else: 
            # position in reference, position in query ,  length of the match 
            values = line.strip().split()  #  ['G98_SE011,'1102,'1114, 21']
            len_mums.append(int(values[-1]))  #  ['1102,'1114, 21']

    return len_mums

def parse_delta(filename): # -> Tuple[int, int, float, int]: # from pyani package 
    """Return (reference alignment length, query alignment length, average identity, similarity erors)
    :param filename: Path to the input .delta file
    Calculates the aligned lengths for reference and query and average nucleotide identity, and returns the cumulative total for each as a tuple.
    The delta file format contains seven numbers in the lines of interest: see http://mummer.sourceforge.net/manual/ for specification
    - start on query
    - end on query
    - start on target
    - end on target
    - error count (non-identical, plus indels)
    - similarity errors (non-positive match scores) [NOTE: with PROmer this is equal to error count]
    - stop codons (always zero for nucmer)
    We report ANIm identity by finding an average across all alignments using the following formula:
    sum of weighted identical bases / sum of aligned bases from each fragment
    For example: reference.fasta query.fasta NUCMER
    >ref_seq_A ref_seq_B 40 40
    1 10 1 11 5 5 0
    -1
    0
    15 20 25 30 0 0 0

    The delta file tells us there are two alignments. The first alignment runs from base 1 to base 10 in the reference sequence, and from base 1 to 11 in the query sequence with a similarity error of 5. The second alignment runs from base 15 to 20 in the reference, and base 25 to 30 in the query with 0 similarity errors. To calculate the %ID, we can:
    - Find the number of all aligned bases from each sequence:
    aligned reference bases region 1 = 10 - 1 + 1 = 10
    aligned query bases region 1 = 11 - 1 + 1 = 11
    aligned reference bases region 2 = 20 - 15 + 1 = 6
    aligned query bases region 2 = 30 - 25 + 1 = 6
    - Find weighted identical bases
    alignment 1 identity weighted = (10 + 11) - (2 * 5) = 11
    alignment 2 identity weighted = (6 + 6) - (2 * 0) = 12
    - Calculate %ID
    (11 + 12) / (10 + 11 + 6 + 6) = 0.696969696969697
    To calculate alignment lengths, we extract the regions of each alignment (either for query or reference) provided in the .delta file and merge the overlapping regions with IntervalTree. Then, we calculate the total sum of all aligned regions.
    """

    current_ref, current_qry, raln_length, qaln_length, sim_error, avrg_ID = ( None, None,0,0, 0, 0.0, )
    regions_ref = defaultdict(list)  # Hold a dictionary for query regions
    regions_qry = defaultdict(list)  # Hold a dictionary for query regions
    aligned_bases = []  # Hold a list for aligned bases for each sequence
    weighted_identical_bases = []  # Hold a list for weighted identical bases
    file_handle= open(filename,"r")
    file_lines_split= [_.strip().split() for _ in file_handle.readlines()]
    for line in file_lines_split:
        if line[0] == "NUCMER":  # Skip headers
            continue # Lines starting with ">" indicate which sequences are aligned
        if line[0].startswith(">"):
            current_ref = line[0].strip(">")
            current_qry = line[1]
        
        if len(line) == 7: # Lines with seven columns are alignment region headers: # Obtaining aligned regions needed to check for overlaps
            regions_ref[current_ref].append( tuple(sorted(list([int(line[0]), int(line[1])]))) )  # aligned regions reference
            regions_qry[current_qry].append( tuple(sorted(list([int(line[2]), int(line[3])]))) )  # aligned regions qry
            # Calculate aligned bases for each sequence
            ref_aln_lengths = abs(int(line[1]) - int(line[0])) + 1
            qry_aln_lengths = abs(int(line[3]) - int(line[2])) + 1
            aligned_bases.append(ref_aln_lengths)
            aligned_bases.append(qry_aln_lengths)
            # Calculate weighted identical bases
            sim_error += int(line[4])
            weighted_identical_bases.append( (ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4])) )
    # Calculate average %ID
    if aligned_bases: # when no aligmnet in the nucmer file 
        avrg_ID = sum(weighted_identical_bases) / sum(aligned_bases)    
    # Calculate total aligned bases (no overlaps)
    for seq_id in regions_qry:
        qry_tree = intervaltree.IntervalTree.from_tuples(regions_qry[seq_id])
        qry_tree.merge_overlaps(strict=False)
        for interval in qry_tree:
            qaln_length += interval.end - interval.begin + 1
    for seq_id in regions_ref:
        ref_tree = intervaltree.IntervalTree.from_tuples(regions_ref[seq_id])
        ref_tree.merge_overlaps(strict=False)
        for interval in ref_tree:
            raln_length += interval.end - interval.begin + 1
    return (raln_length, qaln_length, avrg_ID, sim_error)

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



def read_matrix(filename, extenstion=".fa"):
    sample_names = []
    dic_distance = {}
    file_handle =  open(filename,'r')
    line_counter = -2
    for line in file_handle:
        line_counter += 1
        if line.startswith("\t"):  # skip sample name line
            sample_names_raw = line.strip().split('\t')
            sample_names = [i.split("/")[-1]+extenstion for i in sample_names_raw] # .split(".")[0]
            continue        
        elements = line.strip().split('\t')       
        for i in range(line_counter+1,len(sample_names)):
            s_min=min(sample_names[line_counter],sample_names[i])
            s_max=max(sample_names[line_counter],sample_names[i])
            dic_distance[(s_min,s_max)]=float(elements[i+1])    
    #print(elements)
    return dic_distance
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
    #print(elements)
    return dic_distance
def read_fastani(filename):
    dic_distance = {}
    file_handle =  open(filename,'r')
    for line in file_handle:            
        elements = line.strip().split('\t')
        sp1= elements[0].split("/")[-1]#.split(".")[0]
        sp2= elements[1].split("/")[-1]#.split(".")[0]
        ani= float(elements[2])/100
        if sp1!=sp2:
            s_min=min(sp1,sp2)
            s_max=max(sp1,sp2)            
            dic_distance[(s_min,s_max)]=ani
    #print(elements)
    return dic_distance
    
def read_ortho(filename):
    file_handle =  open(filename,'r')
    found=False
    for line in file_handle:   
        if "OrthoANI" in line: # OrthoANI : 75.1627 (%)
            ani=float(line.split(":")[1].split("(")[0].strip())
            found=True
    if not found:
        ani=0
        print("no line with OrthoANI in it for ", filename)
    return ani

def dist_tree_calc_obj(tree1,extension=".fa"): # , quoted_node_names=False)
    #tree1= Tree(filename1,format=1) # , quoted_node_names=quoted_node_names
    #print(filename1, len(tree1))    
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
    print(num_sample,num_sample*(num_sample-1)/2, len(dist_tree_top))
    return dist_tree


def dist_tree_calc(filename1,extension=".fa"): # , quoted_node_names=False)
    tree1= Tree(filename1,format=1) # , quoted_node_names=quoted_node_names
    print(filename1, len(tree1))    
    dist_tree= dist_tree_calc_obj(tree1,extension)
    return dist_tree


def nucmer_calc(folder, l_list, fastas, tot_lengths):
    dic_ani_dic={}
    dic_cov_dic={}
    perc_ids_dic={}
    dic_anicov_dic = {}
    for l in  l_list: # [8,10,"12","13"  ,"14", "16",  "20"]
        print("l",l)
        l=str(l)
        num_sp= len(fastas)
        perc_ids= []
        dic_ani={}
        dic_cov={}
        dic_anicov= {}
        for i in range(num_sp): # 
            spi=fastas[i]
            if i%30==0: print(i)
            for j in range(i): # num_sp i 
                try:
                    spj=fastas[j]
                    #try:
                    deltafile= folder+spj+"_"+spi+"_l"+l+".filter"
                    (query_tot_length, subject_tot_length, weighted_identity, tot_sim_error) = parse_delta(deltafile)
                    # except:
                    #     deltafile= folder+spi+"_"+spj+"_l"+l+".filter"
                    #     (query_tot_length, subject_tot_length, weighted_identity, tot_sim_error) = parse_delta(deltafile)        
                    query_cover = float(query_tot_length) / tot_lengths[spj] # spj=a.aa for sp=a.aa.fa
                    sbjct_cover = float(subject_tot_length) / tot_lengths[spi]
                    perc_id = weighted_identity
                    perc_ids.append(perc_id)
                    #assert abs(query_cover - sbjct_cover)< 0.2, spi+spj                 #print(spi,spj, perc_id, query_cover , sbjct_cover)
                    s_min=min(spi,spj)
                    s_max=max(spi,spj)            
                    dic_ani[(s_min,s_max)]= perc_id 
                    #print(perc_id)
                    dic_cov[(s_min,s_max)]= query_cover
                    dic_anicov[(s_min,s_max)]= perc_id* query_cover
                except :
                    print("not found",spi,spj)
        perc_ids_dic["l"+str(l)]=perc_ids
        dic_ani_dic["l"+str(l)]=dic_ani
        dic_cov_dic["l"+str(l)]=dic_cov
        dic_anicov_dic["l"+str(l)]=dic_anicov
    print( len(dic_ani_dic),len(dic_ani_dic["l"+str(l)]) )          
    return dic_ani_dic, dic_cov_dic, dic_anicov_dic


def mum_calc(folder, fastas,m,l, tot_lengths):    
    dic_mum_sum_dbygenome = {} 
    dic_mum_sum_dbymum={}
    dic_mum_len={}
    num_sp= len(fastas)
    print(num_sp,fastas[0])
    for i in range(num_sp): 
        spi=fastas[i] # +".fna"
        if i%10==0: print(i)
        for j in range(i): 
            len_mums=[]
            spj=fastas[j] # +".fna"
            mumfile= folder+spi+"_"+spj+"_"+m+"l"+str(l)+".mums"
            mumfile_= folder+spj+"_"+spi+"_"+m+"l"+str(l)+".mums"
            if os.path.exists(mumfile):
                len_mums = parse_mums2(mumfile)
            elif os.path.exists(mumfile_):
                len_mums = parse_mums2(mumfile_)      
            #sum_len_mum= sum(len_mums)
            
            avg_len= (tot_lengths[spj] +tot_lengths[spi] )/2
            #print(avg_len)
            s_min=min(spi,spj)
            s_max=max(spi,spj)            
            dic_mum_sum_dbygenome[(s_min,s_max)]= sum(len_mums)/ avg_len
            if len_mums:
                dic_mum_sum_dbymum[(s_min,s_max)]= sum(len_mums)/len(len_mums)
            else:
                dic_mum_sum_dbymum[(s_min,s_max)]=0
            dic_mum_len[(s_min,s_max)]= len(len_mums)
    print(len(dic_mum_sum_dbygenome))  
    return dic_mum_sum_dbygenome,dic_mum_sum_dbymum,dic_mum_len


def ortho_calc(folder,fastas):
    dic_ortho = {} 
    num_sp= len(fastas)
    #print(num_sp,fastas[0])
    
    for i in range(num_sp): 
        spi=fastas[i] # +".fna"
        if i%30==0: print(i)
        for j in range(i): 
            spj=fastas[j] # +".fna"
            ani=0
            try:
                file= folder+spi+"_"+spj+".out"
                ani = read_ortho(file)
            except:
                try:
                    file= folder+spj+"_"+spi+".out"
                    ani = read_ortho(file)
                except:
                    print("not found", spi, spj)
        
            s_min=min(spi,spj)
            s_max=max(spi,spj)    
            dic_ortho[(s_min,s_max)  ]= ani/100
    print(len(dic_ortho))      
    return dic_ortho