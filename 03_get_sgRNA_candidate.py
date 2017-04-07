########################################
#1. delete those gRNAs with >=4T, >=5A,C,G and >6 dinucleotide or trinucleotide repeats
#2. delete those gRNAs containing restriction enzyme site, eg. BsaI sequences: GAGACC
#command: python step3_sgRNA_filter.py 02.sgRNA 03.filter -enzys enzymes(can be no or GAGACC or GAGACC,GGTCTC or e1,e2,...en)
########################################
import sys,os
enzymes=sys.argv[4].split(',') 



di=[] #create di-repeats
for base1 in ['A','T','C','G']:
    for base2 in ['A','T','C','G']:
        if base1!=base2:
            di.append(base1+base2)
tri=[] #create tri-repeats
for base1 in ['A','T','C','G']:
    for base2 in ['A','T','C','G']:
        for base3 in ['A','T','C','G']: 
            if base1!=base2 or base2!=base3:
                tri.append(base1+base2+base3)
                
discard=[] #gRNA with low repeats
infile=open(sys.argv[1],"r")
lines=infile.read().split('\n')
for line in lines:
    if len(line)>0:
        gRNA_seq=line.split('\t')[0]
        for n in range(4,21):
            if 'T'*n in gRNA_seq or 'A'*(n+1) in gRNA_seq or 'C'*(n+1) in gRNA_seq or 'G'*(n+1) in gRNA_seq:
                discard.append(gRNA_seq)
                break
        for n in range(6,10):
            for base in di:
                if base*n in gRNA_seq:
                    discard.append(gRNA_seq)
                    break
        for n in range(5,7):
            for base in tri:
                if base*n in gRNA_seq:
                    discard.append(gRNA_seq)
                    break
infile.close()



discard2=[] #gRNA with restriction enzyme site
infile=open(sys.argv[1],"r")
lines=infile.read().split('\n')
for line in lines:
    if len(line)>0:
        gRNA_seq=line.split('\t')[0]
        for enzy in enzymes: 
            if enzy in gRNA_seq:
                discard2.append(gRNA_seq)
                break
infile.close()
            


infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
lines=infile.read().split('\n')
for line in lines:
    if len(line)>0:
        gRNA_seq=line.split('\t')[0]
        if gRNA_seq not in discard: #not contain low repeats
            if gRNA_seq not in discard2: #not contain restriction enzyme site
                outfile.write(line+'\n')
infile.close(); outfile.close()




