########################################
#command: python ontarget_phast.py .bed .ontarget.phast -g hg19...
########################################
import sys,os
file_prefix=sys.argv[1].split('.')[0]
genome=sys.argv[4]
exon_for_ontarget='/data8/xiexw/gRNA_tool/genomes/'+genome+'.on.exon'       #modifiable
phastCE='/data8/xiexw/gRNA_tool/genomes/'+genome+'.phast'                   #modifiable
intersectBed='intersectBed'                                                 #modifiable
groupBy='groupBy'                                                           #modifiable



#get on-target site of gRNA
os.system(intersectBed+' -a '+file_prefix+'.bed -b '+exon_for_ontarget+' -wa -wb -loj > '+file_prefix+'.ontarget')



#whether gRNA overlap with phast conserved element
infile=open(file_prefix+".ontarget","r"); outfile=open(file_prefix+".ontarget.bed","w")
i=1
while i!=0:
    line=infile.readline()
    if len(line)>0:
        items=line.strip('\n').split('\t')
        outfile.write('\t'.join(items[0:4])+'|~|'+items[9]+'\t'+items[4]+'\t'+items[5]+'\n')
    else:
        i=0
infile.close(); outfile.close()

os.system(intersectBed+' -a '+file_prefix+'.ontarget.bed -b '+phastCE+' -wa -wb -loj > '+file_prefix+'.phast')
os.system(groupBy+' -i '+file_prefix+'.phast -g 1,2,3,4 -c 11 -o sum,count > '+file_prefix+'.phast.sc')   ##because one gRNA may overlap with several phastCEs

infile=open(file_prefix+".phast.sc","r"); outfile=open(file_prefix+".ontarget.phast","w")
i=1
while i!=0:
    line=infile.readline()
    if len(line)>0:
        items=line.strip('\n').split('\t')
        outfile.write('\t'.join(items[0:4])+'\t'+'('+items[4]+','+items[5]+')'+'\n')
    else:
        i=0
infile.close(); outfile.close()
os.remove(file_prefix+'.ontarget'); os.remove(file_prefix+'.ontarget.bed'); os.remove(file_prefix+'.phast'); os.remove(file_prefix+'.phast.sc')



