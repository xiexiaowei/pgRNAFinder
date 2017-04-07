########################################
#1. get fasta sequence
#command: python step1_get_fasta.py input.txt 01.fa -g hg19,hg38... -f gene_name,genomic_region,target_sequence
########################################
import sys,os
genome=sys.argv[4]
flag=sys.argv[6]
gene_pos='/data8/xiexw/gRNA_tool/genomes/'+genome+'.pos'       #modifiable 
genome_fa='/data8/xiexw/gRNA_tool/genomes/'+genome+'.genome'   #modifiable
fastaFromBed='fastaFromBed'                                    #modifiable



def convert(fragment):   #chr:start-end into chr\tstart\tend
    chr=fragment.split(':')[0]
    start=fragment.split(':')[1].split('-')[0]
    end=fragment.split(':')[1].split('-')[1]
    return chr+'\t'+start+'\t'+end



#get fasta for gene_name
if flag=='gene_name':
    gene_list=','.join(open(sys.argv[1],"r").read().split('\n')).upper().split(',')
 
    infile=open(gene_pos,"r"); outfile=open("gene.bed","w")
    lines=infile.read().split('\n')
    for line in lines:
        if len(line)>0:
            items=line.split('\t')
            if items[0].upper() in gene_list:
                outfile.write(convert(items[1].strip('\r'))+'\t'+items[1].strip('\r')+'#'+items[0].upper()+'\t'+'1'+'\t'+'.'+'\n')
    infile.close(); outfile.close()
    os.system(fastaFromBed+' -fi '+genome_fa+' -bed gene.bed -fo '+sys.argv[2]+' -name')
    os.remove('gene.bed')
    
    

#get fasta for genomic_region     
if flag=='genomic_region':
    infile=open(sys.argv[1],"r"); outfile=open("gene.bed","w")
    lines=infile.read().split('\n')
    for line in lines:
        if len(line)>0:
            outfile.write(convert(line.strip())+'\t'+line.strip()+'\t'+'1'+'\t'+'.'+'\n')
    infile.close(); outfile.close()
    os.system(fastaFromBed+' -fi '+genome_fa+' -bed gene.bed -fo '+sys.argv[2]+' -name')
    os.remove('gene.bed')
 
       

#get fasta for target_sequence       
if flag=='target_sequence':
    infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
    lines=infile.read().split('\n')
    for i in range(len(lines)):
        if len(lines[i])>0:
            if i==0:
                outfile.write(lines[i]+'\n') 
            else:
                if '>' in lines[i]:
                    outfile.write('\n'+lines[i].strip()+'\n')   
                else:
                    outfile.write(lines[i].strip())
    outfile.write('\n')
    infile.close(); outfile.close()



