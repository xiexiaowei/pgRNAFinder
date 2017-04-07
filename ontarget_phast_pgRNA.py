########################################
#1. judge whether knock-out regions locate in exon region and overlap with phast conserved element
#command: python ontarget_phast_pgRNA.py 06.pre.otnum 06.otnum -g genome -f flag
########################################
import sys,os
genome=sys.argv[4]
flag=sys.argv[6]



if flag!='target_sequence':
    infile=open(sys.argv[1],"r"); outfile=open("06.bed","w")
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
	    pos=[]
            items=line.strip('\n').split('\t'); chr=items[0].split(':')[0]
	    pos.append(int(items[3].split('-')[0])); pos.append(int(items[3].split('-')[1]))
	    pos.append(int(items[4].split('-')[0])); pos.append(int(items[4].split('-')[1]))
	    pos.sort()
	    start=pos[1]; end=pos[2]
	    outfile.write(str(chr)+'\t'+str(start)+'\t'+str(end)+'\t'+'||'.join(items)+'\t'+'1'+'\t'+'.'+'\n')
	else:
	    i=0
    infile.close(); outfile.close()

    #judge whether knock-out regions locate in exon region and overlap with phast conserved element
    os.system('python ontarget_phast.py 06.bed 06.ontarget.phast -g '+genome)
    infile=open("06.ontarget.phast","r"); outfile=open(sys.argv[2],"w")
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
            items=line.strip('\n').split('\t'); ontarget=items[3].split('|~|')[1]
	    infos=items[3].split('|~|')[0].split('||'); pair_dis=int(items[2])-int(items[1])
            if 'ENSG' in ontarget and '(-1,1)' not in items[4]:
		outfile.write('\t'.join(infos[0:11])+'\t'+ontarget+'//'+items[4]+'\t'+'\t'.join(infos[11:19])+'\t'+str(pair_dis)+'\n')
            elif 'ENSG' in ontarget and '(-1,1)' in items[4]:
                outfile.write('\t'.join(infos[0:11])+'\t'+ontarget+'//'+'non-phastCE'+'\t'+'\t'.join(infos[11:19])+'\t'+str(pair_dis)+'\n')
            elif 'ENSG' not in ontarget and '(-1,1)' not in items[4]:
                outfile.write('\t'.join(infos[0:11])+'\t'+'non-exonic'+'//'+items[4]+'\t'+'\t'.join(infos[11:19])+'\t'+str(pair_dis)+'\n')
            else:
                outfile.write('\t'.join(infos[0:11])+'\t'+'non-exonic'+'//'+'non-phastCE'+'\t'+'\t'.join(infos[11:19])+'\t'+str(pair_dis)+'\n')
	else:
	    i=0
    infile.close(); outfile.close()
    os.remove('06.bed'); os.remove('06.ontarget.phast')



#target_sequence has no position information, so the on-target site and phastCE are unknown.
else:   ##for target_sequence
    infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
	    pos=[]
            items=line.strip('\n').split('\t')
	    pos.append(int(items[3].split('-')[0])); pos.append(int(items[3].split('-')[1]))
	    pos.append(int(items[4].split('-')[0])); pos.append(int(items[4].split('-')[1]))
	    pos.sort()
	    pair_dis=int(pos[2])-int(pos[1])
            outfile.write('\t'.join(items[0:11])+'\t'+'unknown'+'//'+'unknown'+'\t'+'\t'.join(items[11:19])+'\t'+str(pair_dis)+'\n')
	else:
	    i=0
    infile.close(); outfile.close()
     




