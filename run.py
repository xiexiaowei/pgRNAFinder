########################################           
#1. read config
#2. process
#3. selection
#4. remove intermediate files
#command: python run.py -i infile -o outfile -g genome -f flag -pam PAM -ul gRNA_length -dl downstream_length -enzys enzymes \
#-s score -mis mismatch -gRNA sg,pg -ontarget exon,all \
#-ota ot_num -ote ot_num_in_exon -cg CGmin~CGmax -dis min~max -cover min~max -strand same,diff,all -otsite yes,no
########################################
import sys,os
infile=sys.argv[2]; outfile=sys.argv[4]; genome=sys.argv[6]; flag=sys.argv[8]
PAM=sys.argv[10]; gRNA_length=sys.argv[12]; downstream_length=sys.argv[14]; enzymes=sys.argv[16]
score=sys.argv[18]
mismatch=sys.argv[20]
gRNA=sys.argv[22]
ontarget=sys.argv[24]
ot_num=sys.argv[26]
ot_num_in_exon=sys.argv[28]
CG=sys.argv[30]
if gRNA=='sg':
    otsite=sys.argv[32]
else:
    pair_dis=sys.argv[32] #only for pair
    coverage=sys.argv[34] #only for pair
    strand=sys.argv[36]   #only for pair
    otsite=sys.argv[38]    





print 'Process start<br>' 
#for hg19,hg38,mm10 + bosTau,canFam,danRer,galGal,ratNor,sacCer,susScr;  still can add other genome
os.system('python 01_get_sequence_for_search.py '+infile+' 01.fa -g '+genome+' -f '+flag)
print 'Step1: getting fasta done!<br>'

#for NGG,NAG,NNGRRT,NNNNACA and self-defined
os.system('python 02_search_sgRNA.py 01.fa 02.sgRNA -f '+flag+' -pam '+PAM+' -ul '+gRNA_length+' -dl '+downstream_length) 
print 'Step2: getting sgRNA done!<br>'

#gRNA filter
os.system('python 03_get_sgRNA_candidate.py 02.sgRNA 03.filter -enzys '+enzymes)
print 'Step3: sgRNA filteration done!<br>'

#count SSC score + ontarget in exon + phast conserved element
os.system('python 04_sgRNA_efficiency.py 03.filter 04.ssc.ontarget.phast -g '+genome+' -f '+flag+' -pam '+PAM+' -ul '+gRNA_length+' -dl '+downstream_length+' -s '+score)
print 'Step4: scoring for sgRNA by SSC, judging whether on-target sites locate in exon region and phast conserved element done!<br>'
    
#get offtarget sites; 4 PAMs for all genomes      
if PAM in ['NGG','NAG','NNGRRT','NNNNACA']:     
    os.system('python 05_offtarget_for_sgRNA.py 04.ssc.ontarget.phast 05.otsite -g '+genome+' -pam '+PAM+' -ul '+gRNA_length+' -dl '+downstream_length+' -mis '+mismatch+' -ota '+ot_num)
    print 'Step5: getting off-target sites done!<br>'
else:
    print 'Step5: can not count off-target sites!<br>'
    
#count offtarget sites    
os.system('python 06_count_offtarget.py 05.otsite 06.otnum -g '+genome+' -f '+flag+' -gRNA '+gRNA)
print 'Step6: counting off-target for sgRNA or pgRNA done!<br>'

#selection
if gRNA=='sg':
    os.system('python 07_select_output.py 06.otnum '+outfile+' -gRNA '+gRNA+' -ontarget '+ontarget+' -ota '+ot_num+' -ote '+ot_num_in_exon+' -cg '+CG+' -otsite '+otsite)
else:    
    os.system('python 07_select_output.py 06.otnum '+outfile+' -gRNA '+gRNA+' -ontarget '+ontarget+' -ota '+ot_num+' -ote '+ot_num_in_exon+' -cg '+CG+' -dis '+pair_dis+' -cover '+coverage+' -strand '+strand+' -otsite '+otsite)
print 'Step7: selection for sgRNA or pgRNA done!<br>'
print 'Process done'



os.remove('01.fa'); os.remove('02.sgRNA'); os.remove('03.filter')
os.remove('04.ssc.ontarget.phast'); os.remove('05.otsite'); os.remove('06.otnum')
print 'The end'

