files="sample_AAACGG.fq
sample_AACGTT.fq
sample_AACTGA.fq
sample_AAGACG.fq
sample_AAGCTA.fq
sample_AATATC.fq
sample_AATGAG.fq
sample_ACAAGA.fq
sample_ACAGCG.fq
sample_ACATAC.fq
sample_ACCATG.fq
sample_ACCCCC.fq
sample_ACTCTT.fq
sample_ACTGGC.fq
sample_AGCCAT.fq
sample_AGCGCA.fq
sample_AGGGTC.fq
sample_AGGTGT.fq
sample_AGTAGG.fq
sample_AGTTAA.fq
sample_ATAGTA.fq
sample_ATCAAA.fq
sample_ATGCAC.fq
sample_ATGTTG.fq
sample_ATTCCG.fq
sample_CAAAAA.fq
sample_CAATCG.fq
sample_CACCTC.fq
sample_CAGGCA.fq
sample_CATACT.fq
sample_CCATTT.fq
sample_CCCGGT.fq
sample_CCCTAA.fq
sample_CCGAGG.fq
sample_CCGCAT.fq
sample_CCTAAC.fq
sample_CGAGGC.fq
sample_CGCAGA.fq
sample_CGCGTG.fq
sample_CGGTCC.fq
sample_CGTCTA.fq
sample_CGTGAT.fq
sample_CTACAG.fq
sample_CTCGCC.fq
sample_CTGCGA.fq
sample_CTGGTT.fq
sample_CTTATG.fq
sample_CTTTGC.fq
sample_TCTTCT.fq
sample_TGAACC.fq
sample_TGACAA.fq
sample_TGCCCG.fq
sample_TGCTTA.fq
sample_TGGGGA.fq
sample_TTATGA.fq
sample_TTCCGT.fq
sample_TTCTAG.fq
sample_TTGAGC.fq
sample_TTTAAT.fq
sample_TTTGTC.fq"

for file in $files
do
bowtie2 --sensitive -x Ssco_bowtieRef.fa \
-S ${file}_2013-24-10_SscoAlign.sam -t \
-U ${file}
--un /home/joneslab/Documents/Sarahs/PstI_Data/alignments/${file}_2013-15-08_unpaired.sam \
-p 6

done
