zcat /home/mk446/mutanno/PRECALVEP/chr1/0/chr1_1_100000.tsv.gz | bgzip -c > /home/mk446/mutanno/PRECALVEP/vep.hg38.1.tsv.gz;
zcat /home/mk446/mutanno/PRECALVEP/chr1/0/chr1_100001_200000.tsv.gz | grep -v '^#' | bgzip -c >> /home/mk446/mutanno/PRECALVEP/vep.hg38.1.tsv.gz;
zcat /home/mk446/mutanno/PRECALVEP/chr1/0/chr1_200001_300000.tsv.gz | grep -v '^#' | bgzip -c >> /home/mk446/mutanno/PRECALVEP/vep.hg38.1.tsv.gz;
