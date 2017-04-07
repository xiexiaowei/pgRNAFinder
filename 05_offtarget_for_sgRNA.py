########################################
#run off-targets for hg19...
#command: python step5_get_offtarget.py 04.ssc.ontarget.phast 05.otsite -g genome -pam PAM -ul gRNA_length -dl downstream_length -mis mismatch -ota number
########################################
import sys,os 
genome=sys.argv[4]
PAM=sys.argv[6]
gRNA_length=sys.argv[8]
downstream_length=sys.argv[10]
mismatch=sys.argv[12]
ot_num_cutoff=sys.argv[14]
Results='/data8/xiexw/Off-Spotter/Results'                                #modifiable
OT_dir='/data8/xiexw/Off-Spotter/'+genome+'/'                             #modifiable
exon_for_offtarget='/data8/xiexw/gRNA_tool/genomes/'+genome+'.off.exon'   #modifiable



OUT=open(sys.argv[2],"w")
#run off-spotter; using w303 to replace the names of other 7 genomes
for info in open(sys.argv[1],"r").read().split('\n'):   ##one file for one sgRNA
    if len(info)>0:
        outfile=open("gRNA_20bp.txt","w"); gRNA_20bp=info.split('\t')[0][int(gRNA_length)-20:int(gRNA_length)]; outfile.write(gRNA_20bp+'\n'); outfile.close()
        if genome=='hg19' or genome=='hg38' or genome=='mm10':
            os.system(Results+' -i gRNA_20bp.txt -f -p '+PAM+' -n '+mismatch+' -g '+genome+' -t '+OT_dir+' -a '+exon_for_offtarget+' -o ot.txt -m 1')
        else:
            os.system(Results+' -i gRNA_20bp.txt -f -p '+PAM+' -n '+mismatch+' -g w303 -t '+OT_dir+' -a '+exon_for_offtarget+' -o ot.txt -m 1')        



        #only leave 20bp-gRNA,chr,pos from off-spotter results + sort by gRNA (must do for the latter step)
        infile=open("ot.txt","r"); outfile=open("ot.sum","w")
        lines=infile.read().split('\n')
        sort=[]
        for line in lines:
            if len(line)>0:
                items=line.split('\t')
                exons=[]   ##for delete repeated exon
                for exon in items[7].strip('~').split('~'):
                    if exon.strip('\r') not in exons:
                        exons.append(exon.strip('\r'))
                if '-' in exons:
                    sort.append([items[4],'chr'+items[0],int(items[2]),items[4]+'\t'+'chr'+items[0]+'\t'+items[2]+'\t'+'non-exonic'])
                else:
                    sort.append([items[4],'chr'+items[0],int(items[2]),items[4]+'\t'+'chr'+items[0]+'\t'+items[2]+'\t'+'~'.join(exons)])
        sort.sort()
        for line in sort:
            outfile.write(line[3]+'\n')        
        infile.close(); outfile.close()


        
        #stat ot sites for each gRNA_20bp; dic: {gRNA,{chr,pos}}      
        on_target_info=[] 
        chr=info.split('\t')[4].split(':')[0]; ontarget=info.split('\t')[6].split('//')[0]
        if info.split('\t')[3]=='+':
            on_target_info.append(gRNA_20bp+'\t'+chr+'\t'+str(int(info.split('\t')[1])+1+int(gRNA_length)-20)+'\t'+ontarget)
        else:
            on_target_info.append(gRNA_20bp+'\t'+chr+'\t'+str(int(info.split('\t')[1])+1+int(downstream_length))+'\t'+ontarget)
        
        infile=open("ot.sum","r")
        lines=infile.read().split('\n') 
        dic={}
        for line in lines:
            if len(line)>0:
                if line not in on_target_info:   ##on_target_info: list with ontarget sites, otherwise offtarget number will be +1
                    items=line.split('\t')
                    if items[0] not in dic:       
                        dic2={}
                        dic[items[0]]=dic2
                        dic2[items[1]]=items[2]+'__'+items[3]
                    else:
                        if items[1] not in dic2:
                            dic2[items[1]]=items[2]+'__'+items[3]
                        else:
                            dic2[items[1]]+=','+items[2]+'__'+items[3]
        
        #combine ot infos with score infos. Here count off-target is just to reduce the output in 05.otsite
        if gRNA_20bp in dic:
            offtarget_dic=dic[gRNA_20bp]
            chr_s=set(offtarget_dic.keys())
            ot_num=0
            for chr in chr_s:
                sites=offtarget_dic[chr].split(',')
                ot_num+=len(sites)                  
            if float(ot_num)<float(ot_num_cutoff): 
                OUT.write(info+'\t'+str(offtarget_dic)+'\n')
        else:
            OUT.write(info+'\t'+'0'+'\n')
        infile.close(); os.remove('gRNA_20bp.txt'); os.remove('ot.txt'); os.remove('ot.sum')
OUT.close()




