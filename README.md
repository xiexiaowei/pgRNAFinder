# pgRNAFinder
design distance independent paired-gRNA
pgRNAFinder: a software package for designing both sgRNA (single guide RNA), paired nickase and pgRNA (paired guide RNA),
as well as ranking the candidates by evaluating gRNA efficiency and potential off-target sites.
##########################################################################################
Anyone can use the source codes or the excutable file of pgRNAFinder free of charge for non-commercial use. 
For commercial use, please contact the author. Please send bug reports to: xiexw3@mail2.sysu.edu.cn
This README file covers the following topics:

1. Prerequest for pgRNAFinder
2. pgRNAFinder frame
3. Prepare pgRNAFinder input files
4. How to run pgRNAFinder and the usage of parameters
5. pgRNAFinder output files explanation





##########################################################################################
1. Prerequest for pgRNAFinder
1). Python 2.7
2). bedtools v2.25.0
    fastaFromBed,intersectBed and groupBy were used.
3). SSC (SSC0.1)
    website: https://sourceforge.net/projects/spacerscoringcrispr/?source=typ_redirect
    install steps:
    download SSC0.1.tar.gz
    tar -zxvf SSC0.1.tar.gz
    run "make" in unzipped directory
    done
4). Off-Spotter_v0.2.2
    website: https://cm.jefferson.edu/downloads/off-spotter-help/
    install steps:
    download Off-Spotter_v0.2.2.zip
    unzip Off-Spotter_v0.2.2.zip
    run "make" in unzipped directory
    done
    Before running, you must build index by "Table_Creation" like this:
    Table_Creation -i fasta -o directory
    Table_Creation -i /data8/xiexw/gRNA_tool/genomes/hg19.fa -o /data8/xiexw/Off-Spotter/hg19/
5). necessary reference files:
    ①. gennome sequence              
    ②. gene position                 
    ③. exon region for on-target     
    ④. fasta for Off-Spotter         
    ⑤. exon region for Off-Spotter
    ⑥. phast conserved element
    ⑦. gene structure
    Corresponding files of ten genomes have been provided in ftp://222.200.187.83/genomes
    Current ten genomes: hg19, hg38, mm10, bosTau, canFam, danRer, galGal, ratNor, sacCer and susScr.
    More infomation about reference can be seen from "reference_README.txt" in the "Links" page of pgRNAFinder website.





##########################################################################################
2. pgRNAFinder frame
pgRNAFinder is a software package, which contains these scripts as follows: 
1). 01_get_sequence_for_search.py
2). 02_search_sgRNA.py
3). 03_get_sgRNA_candidate.py
4). 04_sgRNA_efficiency.py       (It would call ontarget_phast.py)
5). 05_offtarget_for_sgRNA.py 
6). 06_count_offtarget.py        (Itwould call ontarget_phast_pgRNA.py and ontarget_phast.py)
7). 07_select_output.py 
8). run.py 
Notice: run.py is the main script. All you need to do is running run.py, which will automatically call scripts from 01 to 07.

The goals achieved by each step:
script1. Getting fasta sequence from genomic region, gene name and target sequence.
script2: Getting sgRNA followed by four common PAMs such as NGG, NAG, NNGRRT, NNNNACA and any self-definded PAM comprised of A,G,C,T,N,R.
script3: Then, sgRNA having low-complexity sequence and specific restriction enzyme site will be discarded.
script4. The remaining sgRNA will be submit to SSC for evaluating efficiency (restricted to human and mouse). 
       At the same time, judging whether on-target sites locate in exon region and overlap with phast conserved element.
script5. Getting off-target sites based on Off-Spotter for ten species.
script6. Counting off-targets and off-targets located in exon for sgRNA.
       Combining pgRNA from sgRNA and counting off-targets and off-targets located in exon for pgRNA.
       At the same time, judging whether knockout regions of pgRNA locate in exon region and overlap with phast conserved element.
script7. Counting CG content, pair distance, deletion frequency and coverage. Finally output the suitable pgRNA with user-defined parameters. 





##########################################################################################
3. Prepare pgRNAFinder input 
pgRNAFinder permits three kinds of inputs:
1). genomic_region (chr:start-end)
    chr1:69091-70008
    chr1:38605757-38605866
    chr8:128752630-128753214
    ...
2). gene name 
    OR4F5
    tert
    ...
3). target_sequence
    The file must have the following structure(fasta format). 
    Notice: Header that starts with ">" and is followed by the unique sequence name with no space. 
    >seq1
    CTGC...GACGCT
    >seq2
    CTTTCTTAAAG..
    .....AAAGAGGC





##########################################################################################
4. How to run pgRNAFinder and the usage of parameters
Usage: 
python run.py -i input -o output.txt -g genome -f flag -pam PAM -ul gRNA_length -dl downstream_length -enzys enzymes \
-s SSC_score -mis mismatch -gRNA type -ontarget ontarget -ota OTnum -ote OTnum_in_exon -cg CG -dis distance -cover coverage -strand strand -otsite yes,no

Examples:
1). for----sg
python run.py -i genomic_region.txt -o output.txt -g hg19 -f genomic_region -pam NGG -ul 20 -dl 7 -enzys no \
-s 0.5 -mis 2 -gRNA sg -ontarget all -ota 10 -ote 5 -cg 0.2~0.8 -otsite yes
2). for----pg 
python run.py -i genomic_region.txt -o output.txt -g hg19 -f genomic_region -pam NGG -ul 20 -dl 7 -enzys no \
-s 0.5 -mis 2 -gRNA pg -ontarget all -ota 10 -ote 5 -cg 0.2~0.8 -dis 0~1000 -cover 0.0~1.0 -strand all -otsite yes

Options:
-i <str>	      Input file
-o <str>              Output file
-g <str>              Genome: hg19, hg38, mm10, bosTau, canFam, danRer, galGal, ratNor, sacCer or susScr.  
-f <str>              Flag: genomic_region, gene_name or target_sequence; should be consistent with the format of input.     
-pam <str>            PAM sequence: NGG, NAG, NNGRRT, NNNNACA or self-definded sequence comprised of A,G,C,T,N,R
-ul <int>             gRNA length: the number of base in the upstream of PAM; at better >=20bp for off-target analysis
-dl <int>             Downstream length: the number of base in the downstream of PAM
-enzys <str>          Restriction enzyme：no enzyme or GAGACC or GAGACC,GGTCTC or self-defined like e1,e2,...en
-s <float>            SSC score: for evaluating gRNA efficiency, only SSC score higher than the cutoff will be reported.
-mis <int>            Mismatch of Off-Spotter when finding off-targets(at most): 0, 1, 2, 3, 4, 5
-gRNA <str>           gRNA type: sg(for designing single gRNA) or pg(for designing paired gRNA)
-ontarget <str>       gRNA site: exon(>=1 gRNA located in exon) or all(report all gRNAs)
-ota <int>	      Off-targets located in genome(at most): Only those gRNAs whose Off-targets is lower than the cutoff will be reported. 
-ote <int>	      Off-targets located in exon(at most): Lower is better.
-cg <float>	      Minimum~Maximum value of GC content
-dis <int>	      Minimum~Maximum distance between pgRNA
-cover <float>        Minimum~Maximum coverage: equal to distance between pgRNA divided by the length of sequence
-strand <str>         Strand: same(two gRNAs are in the same stand), diff(two gRNAs are in the different stand), both(report all gRNAs)
-otsite <str>         yes (output detailed off-target sites) or no (don't output detailed off-target sites)

Advised options for better gRNAs: (some of which are also the default options of website)
-ul                   >=20bp is recommended
-s                    >0.5 is recommended
-mis                  default is 2
-gRNA                 default is pg
-ontarget             exon is recommended; As for target sequence, here this option must be "all".
-ota     	      default is 10
-ote                  default is 5
-cg                   0.2~0.8 is recommended
-strand               default: diff for paired nick; both for long pair
-dis                  default: -2~32bp for paired nickase; 200~1000bp for paired guide RNA
-cover                default 0.0~1.0





##########################################################################################
5. pgRNAFinder output files explanation
1) pgRNA:
output example: (We advise you to view the output in excel.)
query_info	gRNA1	gRNA2	pos1	pos2	strand1	strand2	SSC_score1	SSC_score2	on_target1	on_target2	on_target_pair	phastCE1	phastCE2	phastCE_pair	ot_num1	ot_num2	ot_num_pair	ot_num1_in_exon	ot_num2_in_exon	ot_num_pair_in_exon	CG1	CG2	pair_dis	cover	del_freq	ot_site1	ot_site2
chr1:67091-70091	GAGTGTGTAGGACTAAGAAATGGGATTCAG	CAGAGTAGTAAAGAGAAAAGTGGAATTTCC	67351-67381	67378-67408	+	+	0.51	0.54	non-exonic	non-exonic	non-exonic	non-phastCE	non-phastCE	non-phastCE	2	5	1	0	0	0	0.43	0.37	3	0.00	1.00	{'chr19': '108940__non-exonic', 'chr15': '102464544__non-exonic'}	{'chr13': '53296344__non-exonic', 'chr7': '24663568__non-exonic', 'chr19': '108967__non-exonic', 'chr3': '157860263__non-exonic', 'chr1': '210265930__non-exonic'}
chr1:67091-70091	GAGTGTGTAGGACTAAGAAATGGGATTCAG	CTATACCTTCATGTCTCCCGTGGAATGTTA	67351-67381	67490-67520	+	+	0.51	0.56	non-exonic	non-exonic	non-exonic	non-phastCE	non-phastCE	non-phastCE	2	2	2	0	0	0	0.43	0.43	109	0.04	0.89	{'chr19': '108940__non-exonic', 'chr15': '102464544__non-exonic'}	{'chr19': '109079__non-exonic', 'chr15': '102464405__non-exonic'}
...

There are 28 columns in this output file. Their meaning are:	
field                  Meaning
query_info             Name of input sequence
gRNA1                  Base sequence of one gRNA  	
gRNA2                  Base sequence of another pairwise gRNA 
pos1                   Position of gRNA1 (the genome position for genomic_region and gene_name; the relative position in input sequence for target sequence)
pos2                   Position of gRNA2
strand1                Strand of gRNA1
strand2                Strand of gRNA2
SSC_score1             Efficiency score for gRNA1 by SSC software (This has two kinds: score or non-SSC. non-SSC means that gRNA length, PAM and species can't meet the request of SSC software.)
SSC_score2             Efficiency score for gRNA2 by SSC software
on_target1             On-target site of gRNA1    (This has three kinds: Ensemblid_genename_exon or non-exonic or unknown. non-exonic means that gRNA locates in non-coding region; unknown appears when flag is target_sequence.)
on_target2             On-target site of gRNA2
on_target_pair         On-target site of knockout region between pgRNA
phastCE1               Phast conserved element overlapping with gRNA1  (This has three kinds: sum,count or non-phastCE or unknown. non-phastCE means no overlap with any phast conserved element; unknown appears when flag is target_sequence.)
phastCE2               Phast conserved element overlapping with gRNA2
phastCE_pair           Phast conserved element overlapping with knockout region between pgRNA
ot_num1_in_exon	       The number of off-targets located in exon for gRNA1
ot_num2_in_exon	       The number of off-targets located in exon for gRNA2
ot_num_pair_in_exon    The number of off-targets located in exon for pgRNA (>=1 off-target located in exon)
ot_num1	               The number of off-targets for gRNA1
ot_num2	               The number of off-targets for gRNA2
ot_num_pair            The number of off-targets for pgRNA (>=1 off-target located in exon)
CG1	               CG content of gRNA1
CG2	               CG content of gRNA2
pair_dis               Distance between pgRNA
cover	               Coverage of pgRNA     
del_freq               Deletion frequency, this value depends on the distance between pgRNA
ot_site1               Detailed off-target sites for gRNA1	
ot_site2               Detailed off-target sites for gRNA2


2) single gRNA:
output example: (We advise you to view the output in excel.)
query_info	gRNA	pos	strand	SSC_score	on_target	phastCE	ot_num	ot_num_in_exon	CG	ot_site
chr1:69091-70008	CCACTGTTATGACAATAAGAAGGTTTCCAA	69187	-	0.59	ENSG00000186092.4_OR4F5_exon	(1466,5)	3	2	0.37	{'chr19': '110783__ENSG00000176695.5_OR4F17_exon', 'chr3': '122986852__non-exonic', 'chr15': '102463136__ENSG00000177693.3_OR4F4_exon'}
chr1:69091-70008	GTACATGGGAGAGTGAAGGTGGGAGTCAGA	69219	-	0.55	ENSG00000186092.4_OR4F5_exon	(746,2)	5	4	0.53	{'chr7': '142760355__ENSG00000179420.10_OR6W1P_exon', 'chr11': '7795190__ENSG00000166408.3_OR5P1P_exon', 'chr19': '110815__ENSG00000176695.5_OR4F17_exon', 'chr15': '102463104__ENSG00000177693.3_OR4F4_exon', 'chr2': '13424__non-exonic'}
...
Explanation is the similar to above.

