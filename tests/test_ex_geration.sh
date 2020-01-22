gzhead /home/mk446/bio/pipeline/SIMDATA/HG002.triorandom/thr30_GVCF/JC50_HG002.triorandom_r1.chr1.geno.vcf.gz.mnv.vqsr.vcf.gz 100200 > test_trio_JC50.vcf
python mutanno.py annot -vcf test_trio_JC50.vcf -out test_trio_JC50.annot.vcf -ds datastructure_gnomadonly.json
tabixgz test_trio_JC50.annot.vcf
