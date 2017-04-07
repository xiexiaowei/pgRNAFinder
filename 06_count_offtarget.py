########################################
#1. count off-targets for sgRNAs
#2. combine pgRNAs and count off-targets for pgRNAs
#command: python 06_count_offtarget.py 05.otsite 06.otnum -g genome -f flag -gRNA sg,pg
#when deletion_size is 10kb, deletion frequency is about 15%. So we set 10b as maximum.
########################################
import sys,os
genome=sys.argv[4]
flag=sys.argv[6]
gRNA=sys.argv[8]
deletion_size=10000     #modifiable



def count_OT(offtarget_dic):
    ot_num=0; ot_num_in_exon=0 
    chr_s=set(offtarget_dic.keys())
    for chr in chr_s:
        sites=offtarget_dic[chr].split(',')
        ot_num+=len(sites)
        for site in sites:
            if '_exon' in site:
                ot_num_in_exon+=1
    return str(ot_num)+';'+str(ot_num_in_exon)



#count off-targets for sgRNAs
if gRNA=='sg':
    infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
            items=line.strip('\n').split('\t')
            offtarget_dic='0'; ot_num=0; ot_num_in_exon=0 
            
            if '0'!=items[7]:
                offtarget_dic=eval(items[7])   ##eval: change the str into dict
                ot_num=count_OT(offtarget_dic).split(';')[0]
                ot_num_in_exon=count_OT(offtarget_dic).split(';')[1]
            outfile.write('\t'.join(items[0:6])+'\t'+items[6].strip('\r')+'\t'+str(ot_num)+'\t'+str(ot_num_in_exon)+'\t'+str(offtarget_dic)+'\n') 
        else:
            i=0
    infile.close; outfile.close    



#combine pgRNAs and count off-targets for pgRNAs
if gRNA=='pg':
    #stat gRNA infos for each gene or fragment
    infile=open(sys.argv[1],"r"); outfile=open('06.pre.otnum',"w") 
    nlist=[]; gRNA_dic={}; pos_dic={}; strand_dic={}; SSC_dic={}; ontarget_dic={}; offtarget_dic={}
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
            items=line.strip('\n').split('\t') 
            if items[4] not in nlist:
                nlist.append(items[4]) 
                gRNA_dic[items[4]]=items[0]; pos_dic[items[4]]=items[1]+'-'+items[2];strand_dic[items[4]]=items[3]    
                SSC_dic[items[4]]=items[5]; ontarget_dic[items[4]]=items[6]; offtarget_dic[items[4]]=items[7]
            else:
                gRNA_dic[items[4]]+=';'+items[0]; pos_dic[items[4]]+=';'+items[1]+'-'+items[2]; strand_dic[items[4]]+=';'+items[3]
                SSC_dic[items[4]]+=';'+items[5]; ontarget_dic[items[4]]+=';'+items[6]; offtarget_dic[items[4]]+=';'+items[7]
        else:
            i=0
    
    
    #combine pgRNAs and count off-targets for pgRNAs
    for key in nlist:
        gRNAs=gRNA_dic[key].split(';'); pos=pos_dic[key].split(';'); strands=strand_dic[key].split(';')
        SSCs=SSC_dic[key].split(';'); ontargets=ontarget_dic[key].split(';'); offtargets=offtarget_dic[key].split(';')
        for j in range(len(offtargets)):
            for k in range(len(offtargets)):
                if k>j:
                    offtarget_dic1='0'; offtarget_dic2='0'
                    ot_num1=0; ot_num2=0; ot_num_pair=0   ##1,2,pair
                    ot_num1_in_exon=0; ot_num2_in_exon=0; ot_num_pair_in_exon=0
    
                    if '0'!=offtargets[j] and '0'!=offtargets[k]:
                        offtarget_dic1=eval(offtargets[j]); offtarget_dic2=eval(offtargets[k])    ##eval: change the str into dict
                        ot_num1=count_OT(offtarget_dic1).split(';')[0]; ot_num2=count_OT(offtarget_dic2).split(';')[0]
                        ot_num1_in_exon=count_OT(offtarget_dic1).split(';')[1]; ot_num2_in_exon=count_OT(offtarget_dic2).split(';')[1]
                        chr_s1=set(offtarget_dic1.keys()); chr_s2=set(offtarget_dic2.keys())
                        for chr in chr_s1.intersection(chr_s2):
                            sites1=offtarget_dic1[chr].split(',')
                            sites2=offtarget_dic2[chr].split(',')
                            for site1 in sites1:
                                for site2 in sites2:
                                    if abs(int(site2.split('__')[0])-int(site1.split('__')[0]))<int(deletion_size):  
                                        ot_num_pair+=1
                                        if '_exon' in site2.split('__')[1] or '_exon' in site1.split('__')[1]:
                                            ot_num_pair_in_exon+=1
                      
                    if '0'!=offtargets[j] and '0'==offtargets[k]:
                        offtarget_dic1=eval(offtargets[j]) 
                        ot_num1=count_OT(offtarget_dic1).split(';')[0]
                        ot_num1_in_exon=count_OT(offtarget_dic1).split(';')[1]
                                            
                    if '0'==offtargets[j] and '0'!=offtargets[k]:  
                        offtarget_dic2=eval(offtargets[k])
                        ot_num2=count_OT(offtarget_dic2).split(';')[0]
                        ot_num2_in_exon=count_OT(offtarget_dic2).split(';')[1]
                         
                    outfile.write(str(key)+'\t'+gRNAs[j]+'\t'+gRNAs[k]+'\t'+pos[j]+'\t'+pos[k]+'\t'+strands[j]+'\t'+strands[k]+'\t'+SSCs[j]+'\t'+SSCs[k]+'\t'+ontargets[j]+'\t'+ontargets[k]+'\t')   
                    outfile.write(str(ot_num1)+'\t'+str(ot_num2)+'\t'+str(ot_num_pair)+'\t'+str(ot_num1_in_exon)+'\t'+str(ot_num2_in_exon)+'\t'+str(ot_num_pair_in_exon)+'\t'+str(offtarget_dic1)+'\t'+str(offtarget_dic2)+'\n')
    infile.close(); outfile.close()  
    os.system('python ontarget_phast_pgRNA.py 06.pre.otnum '+sys.argv[2]+' -g '+genome+' -f '+flag)
    os.remove('06.pre.otnum')
    





