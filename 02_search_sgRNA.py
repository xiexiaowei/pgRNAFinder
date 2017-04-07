########################################
#1. get gRNA(gRNA_length+PAM+downstream_length; 30bp for SSC,20bp for off-spotter); first for + strand, then for - strand; 
#2. get the upper of gRNA + count sequence length
#command: python step2_get_sgRNA.py 01.fa 02.sgRNA -f flag -pam NGG,NAG,NNGRRT,NNNNACA,self-defined -ul gRNA_length -dl downstream_length
########################################
import sys,os
flag=sys.argv[4]
PAM=sys.argv[6]
gRNA_length=sys.argv[8]
downstream_length=sys.argv[10]



def reverse_complement(seq):
    compl={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    lis=[]
    for base in seq.upper():
        lis.append(compl[base])
    lis.reverse()
    return ''.join(lis)



#PAM is comprised of base: A,G,C,T,N and R. N is A or G or C or T; R is A or G.
#PAM is input-str, pams is middle-set, PAMs is output-list
def get_PAMs(PAM):
    def replace_N(seq):
        for i in range(len(seq)):
            if seq[i]=='N':
                pams.add(''.join(seq[:i])+'A'+''.join(seq[i+1:]))
                pams.add(''.join(seq[:i])+'G'+''.join(seq[i+1:]))
                pams.add(''.join(seq[:i])+'C'+''.join(seq[i+1:]))
                pams.add(''.join(seq[:i])+'T'+''.join(seq[i+1:]))
        return pams
    def replace_R(seq):
        for i in range(len(seq)):
            if seq[i]=='R':
                pams.add(''.join(seq[:i])+'A'+''.join(seq[i+1:]))
                pams.add(''.join(seq[:i])+'G'+''.join(seq[i+1:]))
        return pams    
          
    pams=set(PAM.split(','))
    N=list(PAM).count('N')
    while N!=0:
        for seq in list(pams):
            replace_N(seq)
        N-=1
    R=list(PAM).count('R')
    while R!=0:
        for seq in list(pams):
            replace_R(seq)
        R-=1
    
    PAMs=[]
    for seq in list(pams):
        if 'N' not in seq and 'R' not in seq:
            PAMs.append(seq)    
    return PAMs


    
#get gRNA: end with PAM and start with reverse_complement(PAM)
infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
lines=infile.read().split('\n')
for i in range(len(lines)):
    if len(lines[i])>0:
        if '>' not in lines[i]:
            for j in range(len(lines[i])):
                PAMs=get_PAMs(PAM)
                
                if lines[i][j:j+len(PAM)].upper() in PAMs:  #eg. PAM is NGG, PAMs contains AGG,TGG,CGG,GGG
                    if len(lines[i][int(j)-int(gRNA_length):int(j)+len(PAM)+int(downstream_length)])==int(gRNA_length)+len(PAM)+int(downstream_length):
                        outfile.write(lines[i][int(j)-int(gRNA_length):int(j)+len(PAM)+int(downstream_length)].upper()+'\t')
                        if flag!='target_sequence': 
                            gRNA_start=int(lines[i-1].split(':')[1].split('-')[0])+int(j)-int(gRNA_length)
                            gRNA_end=int(lines[i-1].split(':')[1].split('-')[0])+int(j)+len(PAM)+int(downstream_length)
                            outfile.write(str(gRNA_start)+'\t'+str(gRNA_end)+'\t'+'+'+'\t'+lines[i-1][1:].strip()+'~~'+str(len(lines[i]))+'\n') 
                        else:
                            gRNA_start=int(j)-int(gRNA_length)
                            gRNA_end=int(j)+len(PAM)+int(downstream_length)
                            outfile.write(str(gRNA_start)+'\t'+str(gRNA_end)+'\t'+'+'+'\t'+lines[i-1][1:].strip()+'~~'+str(len(lines[i]))+'\n')                              
                
                if reverse_complement(lines[i][j-len(PAM):j].upper()) in PAMs:
                    if len(lines[i][int(j)-len(PAM)-int(downstream_length):int(j)+int(gRNA_length)])==int(gRNA_length)+len(PAM)+int(downstream_length):	
                        outfile.write(reverse_complement(lines[i][int(j)-len(PAM)-int(downstream_length):int(j)+int(gRNA_length)]).upper()+'\t')
                        if flag!='target_sequence':
                            gRNA_start=int(lines[i-1].split(':')[1].split('-')[0])+int(j)-len(PAM)-int(downstream_length)
                            gRNA_end=int(lines[i-1].split(':')[1].split('-')[0])+int(j)+int(gRNA_length)
                            outfile.write(str(gRNA_start)+'\t'+str(gRNA_end)+'\t'+'-'+'\t'+lines[i-1][1:].strip()+'~~'+str(len(lines[i]))+'\n') 
                        else:
                            gRNA_start=int(j)-len(PAM)-int(downstream_length)
                            gRNA_end=int(j)+int(gRNA_length)                            
                            outfile.write(str(gRNA_start)+'\t'+str(gRNA_end)+'\t'+'-'+'\t'+lines[i-1][1:].strip()+'~~'+str(len(lines[i]))+'\n')                 
infile.close(); outfile.close()
