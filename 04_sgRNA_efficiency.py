########################################
#1. SSC score for 30bp-gRNA
#2. judge whether gRNAs locate in exon region and overlap with phast conserved element
#command: python step4_SSC_ontarget_phast.py 03.filter 04.ssc.ontarget.phast -g genome -f flag -pam PAM -ul gRNA_length -dl downstream_length -s SSC_score
########################################
import sys,os
genome=sys.argv[4]
flag=sys.argv[6]
PAM=sys.argv[8]
gRNA_length=sys.argv[10]
downstream_length=sys.argv[12]
SSC_score=sys.argv[14]
SSC='/data8/xiexw/SSC0.1/bin/SSC'                                           #modifiable
SSC_matrix='/data8/xiexw/SSC0.1/matrix/human_mouse_CRISPR_KO_30bp.matrix'   #modifiable



def BED(in_name,out_name,SSC_score):   ##file into bed format
    infile=open(in_name,"r"); outfile=open(out_name,"w")
    lines=infile.read().split('\n')
    for line in lines:
        if len(line)>0:
            items=line.split('\t')
            if SSC_score!='non-SSC':
                if float(items[5])>float(SSC_score):
                    outfile.write(items[4].split(':')[0]+'\t'+items[1]+'\t'+items[2]+'\t'+'||'.join(items[0:5])+'||'+'%.2f' % float(items[5])+'\t'+'1'+'\t'+'.'+'\n')
            else:
                outfile.write(items[4].split(':')[0]+'\t'+items[1]+'\t'+items[2]+'\t'+'||'.join(items[0:5])+'||'+'non-SSC'+'\t'+'1'+'\t'+'.'+'\n')
    infile.close(); outfile.close() 
    return outfile



#arrange gRNAs locating in exon region and overlaping with phast conserved element
def ARRANGE(in_name,out_name):
    infile=open(in_name,"r"); outfile=open(out_name,"w")
    lines=infile.read().split('\n')
    for line in lines:
        if len(line)>0:
            items=line.split('\t'); ontarget=items[3].split('|~|')[1]
            infos=items[3].split('|~|')[0].split('||')
            if 'ENSG' in ontarget and '(-1,1)' not in items[4]:
                outfile.write('\t'.join(infos[0:6])+'\t'+ontarget+'//'+items[4]+'\n')
            elif 'ENSG' in ontarget and '(-1,1)' in items[4]:
                outfile.write('\t'.join(infos[0:6])+'\t'+ontarget+'//'+'non-phastCE'+'\n')
            elif 'ENSG' not in ontarget and '(-1,1)' not in items[4]:
                outfile.write('\t'.join(infos[0:6])+'\t'+'non-exonic'+'//'+items[4]+'\n')
            else:
                outfile.write('\t'.join(infos[0:6])+'\t'+'non-exonic'+'//'+'non-phastCE'+'\n')   
    infile.close(); outfile.close()
    return outfile



#target_sequence has no position information, so the on-target site and phastCE are unknown.
def ARRANGE_target(in_name,out_name,SSC_score):
    infile=open(in_name,"r"); outfile=open(out_name,"w")
    lines=infile.read().split('\n')
    for line in lines:
        if len(line)>0:
            items=line.split('\t')
            if SSC_score!='non-SSC':
                if float(items[5])>float(SSC_score):   ##SSC_score
                    outfile.write('\t'.join(items[0:5])+'\t'+'%.2f' % float(items[5])+'\t'+'unknown'+'//'+'unknown'+'\n')
            else:
                outfile.write('\t'.join(items[0:5])+'\t'+'non-SSC'+'\t'+'unknown'+'//'+'unknown'+'\n')
    infile.close(); outfile.close() 
    return outfile  
    


#just NGG for hg19,hg38,mm10; score for 30bp gRNA by means of SSC software 
if PAM=='NGG' and genome in ['hg19','hg38','mm10'] and int(gRNA_length)==20 and int(downstream_length)==7:        
    os.system(SSC+' -l 30 -m '+SSC_matrix+' -i '+sys.argv[1]+' -o 03.ssc') 
    
    if flag!='target_sequence':
        BED('03.ssc','03.bed',SSC_score)
        os.system('python ontarget_phast.py 03.bed 03.ontarget.phast -g '+genome)
        ARRANGE('03.ontarget.phast',sys.argv[2])
        os.remove('03.ssc'); os.remove('03.bed'); os.remove('03.ontarget.phast') 
    else: 
        ARRANGE_target('03.ssc',sys.argv[2],SSC_score)
        os.remove('03.ssc')



else:
    if flag!='target_sequence':
        BED(sys.argv[1],'03.bed','non-SSC')
        os.system('python ontarget_phast.py 03.bed 03.ontarget.phast -g '+genome)
        ARRANGE('03.ontarget.phast',sys.argv[2])
        os.remove('03.bed'); os.remove('03.ontarget.phast') 
    else: 
        ARRANGE_target(sys.argv[1],sys.argv[2],'non-SSC')   
    
