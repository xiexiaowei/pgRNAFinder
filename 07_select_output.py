########################################           
#1. count CG + count coverage
#2. selection
#command: python 07_select_output.py 06.otnum output.txt -gRNA sg -ontarget exon,all -ota number -ote number -cg 0.2~0.8 -otsite yes,no
#command: python 07_select_output.py 06.otnum output.txt -gRNA pg -ontarget exon,all -ota number -ote number -cg 0.2~0.8 -dis 0~1000 -cover 0.0~1.0 -strand same,diff,all -otsite yes,no
########################################
import sys,os
gRNA=sys.argv[4]
ontarget=sys.argv[6]
ot_num=sys.argv[8]
ot_num_in_exon=sys.argv[10]
CG=sys.argv[12]
CGmin=CG.split('~')[0]
CGmax=CG.split('~')[1]
if gRNA=='sg':
    otsite=sys.argv[14]
else:
    pair_dis=sys.argv[14]
    coverage=sys.argv[16]
    strand=sys.argv[18]
    otsite=sys.argv[20]  
    pair_dis_min=pair_dis.split('~')[0]
    pair_dis_max=pair_dis.split('~')[1]
    coverage_min=coverage.split('~')[0]
    coverage_max=coverage.split('~')[1]

    

def new_line_sg(line,otsite):
    items=line.strip('\n').split('\t')
    if otsite=='yes':
        new_line=items[4].split('~~')[0]+'\t'+items[0]+'\t'+items[1]+'\t'+items[3]+'\t'+items[5]+'\t'+items[6].split('//')[0]+'\t'+items[6].split('//')[1]+'\t'+items[7]+'\t'+items[8]+'\t'+'%.2f' % float(CG)+'\t'+items[9]+'\n'
    else:
	new_line=items[4].split('~~')[0]+'\t'+items[0]+'\t'+items[1]+'\t'+items[3]+'\t'+items[5]+'\t'+items[6].split('//')[0]+'\t'+items[6].split('//')[1]+'\t'+items[7]+'\t'+items[8]+'\t'+'%.2f' % float(CG)+'\n'
    return new_line

def new_line_pg(line,otsite):
    items=line.strip('\n').split('\t') 
    if float(items[20])<=80.0:
	del_freq=1.0
    else:
        del_freq=((float(items[20])/1000.0)**(-0.36)*41.41-2.84)/100.0   ##deletion size(pair_dis): kb
    if otsite=='yes':
	new_line=items[0].split('~~')[0]+'\t'+'\t'.join(items[1:9])+'\t'+items[9].split('//')[0]+'\t'+items[10].split('//')[0]+'\t'+items[11].split('//')[0]+'\t'+items[9].split('//')[1]+'\t'+items[10].split('//')[1]+'\t'+items[11].split('//')[1]+'\t'+'\t'.join(items[12:18])+'\t'+'%.2f' % float(CG1)+'\t'+'%.2f' % float(CG2)+'\t'+items[20]+'\t'+'%.2f' % float(cover)+'\t'+'%.2f' % float(del_freq)+'\t'+items[18]+'\t'+items[19]+'\n'
    else:
        new_line=items[0].split('~~')[0]+'\t'+'\t'.join(items[1:9])+'\t'+items[9].split('//')[0]+'\t'+items[10].split('//')[0]+'\t'+items[11].split('//')[0]+'\t'+items[9].split('//')[1]+'\t'+items[10].split('//')[1]+'\t'+items[11].split('//')[1]+'\t'+'\t'.join(items[12:18])+'\t'+'%.2f' % float(CG1)+'\t'+'%.2f' % float(CG2)+'\t'+items[20]+'\t'+'%.2f' % float(cover)+'\t'+'%.2f' % float(del_freq)+'\n'
    return new_line



if gRNA=='sg':
    infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
    outfile.write('query_info'+'\t'+'gRNA'+'\t'+'pos'+'\t'+'strand'+'\t'+'SSC_score'+'\t'+'on_target'+'\t'+'phastCE'+'\t'+'ot_num'+'\t'+'ot_num_in_exon'+'\t'+'CG'+'\t'+'ot_site'+'\n')   
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
            items=line.strip('\n').split('\t')
	    CG=float(items[0].count('C')+items[0].count('G'))/float(len(items[0]))
	    
	    if float(items[7])<float(ot_num) and float(items[8])<float(ot_num_in_exon):
	        if float(CGmin)<float(CG)<float(CGmax): 
		    if ontarget=='exon':
		        if '_exon' in items[6]:
			    outfile.write(new_line_sg(line,otsite))
		    else:
			outfile.write(new_line_sg(line,otsite))
        else:
	    i=0
    infile.close(); outfile.close()    
 


if gRNA=='pg':
    infile=open(sys.argv[1],"r"); outfile=open(sys.argv[2],"w")
    outfile.write('query_info'+'\t'+'gRNA1'+'\t'+'gRNA2'+'\t'+'pos1'+'\t'+'pos2'+'\t'+'strand1'+'\t'+'strand2'+'\t'+'SSC_score1'+'\t'+'SSC_score2'+'\t')
    outfile.write('on_target1'+'\t'+'on_target2'+'\t'+'on_target_pair'+'\t'+'phastCE1'+'\t'+'phastCE2'+'\t'+'phastCE_pair'+'\t')
    outfile.write('ot_num1'+'\t'+'ot_num2'+'\t'+'ot_num_pair'+'\t'+'ot_num1_in_exon'+'\t'+'ot_num2_in_exon'+'\t'+'ot_num_pair_in_exon'+'\t') 
    outfile.write('CG1'+'\t'+'CG2'+'\t'+'pair_dis'+'\t'+'cover'+'\t'+'del_freq'+'\t'+'ot_site1'+'\t'+'ot_site2'+'\n')
    i=1
    while i!=0:
        line=infile.readline()
        if len(line)>0:
            items=line.strip('\n').split('\t')
	    CG1=float(items[1].count('C')+items[1].count('G'))/float(len(items[1]))
	    CG2=float(items[2].count('C')+items[2].count('G'))/float(len(items[2]))
	    cover=float(items[20])/float(items[0].split('~~')[1])
		
	    if float(items[12])<float(ot_num) and float(items[13])<float(ot_num) and float(items[14])<float(ot_num):
		if float(items[15])<float(ot_num_in_exon) and float(items[16])<float(ot_num_in_exon) and float(items[17])<float(ot_num_in_exon):
		    if float(CGmin)<float(CG1)<float(CGmax) and float(CGmin)<float(CG2)<float(CGmax):
			if float(pair_dis_min)<float(items[20])<float(pair_dis_max):
			    if float(coverage_min)<float(cover)<float(coverage_max):
				
				if ontarget=='exon':
				    if '_exon' in items[9] or '_exon' in items[10] or '_exon' in items[11]:
					if strand=='same':
					    if items[5]==items[6]:
						outfile.write(new_line_pg(line,otsite))
					if strand=='diff':
					    if items[5]!=items[6]:
						outfile.write(new_line_pg(line,otsite))
					if strand=='all':
					    outfile.write(new_line_pg(line,otsite))				
				else:
				    if strand=='same':
					if items[5]==items[6]:
					    outfile.write(new_line_pg(line,otsite))
				    if strand=='diff':
					if items[5]!=items[6]:
					    outfile.write(new_line_pg(line,otsite))
				    if strand=='all':
					outfile.write(new_line_pg(line,otsite))
        else:
	    i=0
    infile.close(); outfile.close()    


  
                


    
