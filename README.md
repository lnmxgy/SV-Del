# SV-Del
Detecting genomic deletion with complex clonal structure from next generation sequencing data

Introduction
----------------------------------
    SV-Del aims to detect drift deletions and SNVs around fixed deletion and to attain the number of sub-clones by separating SNVs.
    Contact : Yu Geng <gengyu@jzmu.edu.cn>; Zewen Tian <1196815117@qq.com>  
	June 20, 2018

Installation
------------------------------------
1. SV-Del runs on 64-bit Linux system.

2. Depends: samtools(Version: 0.1.19-44428cd), DELLY (Version: 0.7.7), svprops，picard.jar

3. Usage:
   (1)python runR1.py ref.fa tumor.sam(bam) nor.sam(bam) -k k -s s -e e -n n /path/picard.jar 
        required parameters: ref.fa : reference genome
                             tumor.sam(bam)：the bam file of tumor sample
		      			     nor.sam(bam) ： the bam file of normal sample
							 /path/picard.jar ：picard.jar path
        optional parameters: -k  the length of kmer[10]
		                     -s  the kmer supports of candidate position[5]
                             -e  the extension length of constructing SNV reference[5000]
                             -n  at least supports of each SNV when separating SNVs[4]
   (2)chmod -R 777 run1.sh
   (3)command: ./run1.sh

4.Output file:
    (1) drift deletion：BreakPoint.txt
	(2) SNVs：snvall.txt
	(3) separating result：snv.txt

5.Example：
    (1)cd ./example
    (2)python ../runR1.py chr19_300w.fa s150.bam nor150.bam 10 5 5000 4 /root/pindel/GATK/picard.jar
    (3)chmod -R 777 run1.sh
    (4)./run1.sh
