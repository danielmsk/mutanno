python /home/mk446/mutanno/SRC/scripts/precal_vep/s03_merge_vep_by_chrom.py 9 | bgzip -c > /home/mk446/mutanno/PRECALVEP/vep.98.hg38.9.tsv.gz;sleep 120;tabix -f -p vcf /home/mk446/mutanno/PRECALVEP/vep.98.hg38.9.tsv.gz;