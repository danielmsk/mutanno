[TOC]

# Variant Annotation

1. [Conservation](#conservation)
2. [hg19 coordinates](#hg19_coordinates)

## Make Micro-Annotation Source Files

```
python s04_merge_and_tabixgz.py > s04.sh
./s04.sh
```
**Input:**

1. Wigfile for PhastCons, PhyloP, Gerp
2. SiPhy bed file (tabix-indexed)

**Output:**

1. `[DATASOURCE_PATH]/CONSERVATION/conservation_scores.hg38.chr##CHROM##.tsi.gz`

 ```
#CHROM	POS	ID	PHASTCONS100|PHASTCONS30|PHASTCONS20|PHYLOP100|PHYLOP30|PHYLOP20|GERP|SIPHY_OMEGA|SIPHY_P
5	11014		|0.113|||0.146||||
5	11015		|0.106|||0.158||||
5	11016		|0.099|||0.177||||
```




## Conservation

### Merge Scores for 

```
python preproc_convert_and_merge_conservation_wigfix2tsi.py
```

**Input:**

1. PhastCons
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS100WAY/hg38/chr#CHROM#.phastCons100way.wigFix.gz`
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS30WAY/hg38/chr#CHROM#.phastCons30way.wigFix.gz`
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS20WAY/hg38/chr#CHROM#.phastCons20way.wigFix.gz`
1. PhyloP
1. Gerp
1. SiPhy
 * 

**Output:**

1. `[DATASOURCE_PATH]/CONSERVATION/conservation_scores.hg38.chr#CHROM#.tsi`


### SiPhy

```
python ./conservation/preproc_siphy.py
vcf-sort -c [DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt \
   > [DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.sorted.bed 
vcf-sort -c [DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.tsv \
   > [DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.tsv.sorted.bed
```

**Input:**

1. pi score file: `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg19_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz`
1. omega score file: `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg19_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz`

**Interim Output:**

1. `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.tsv.gz`
1. `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.gz`

**Output:**

1. `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_omega_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.txt.sorted.bed.gz`
1. `[DATASOURCE_PATH]/CONSERVATION/SIPHY/29mammals/hg38liftover_29way_pi_lods_elements_12mers.chr_specific.fdr_0.1_with_scores.tsv.sorted.bed.gz`






# Gene Annotation

## Make Datasource file

**For all genes**

```
mutanno makedata -ds ./tests/data/datastructure_v0.4.3ds_mvp.json \
    -out ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3 \
    -vartype GENE \
    -outtype json
gz ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3.json
```

**For coding genes**

```
mutanno makedata -ds ./tests/data/datastructure_v0.4.3ds_mvp.json \
    -out ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3.coding_gene \
    -vartype CODINGGENE \
    -outtype json
gz ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3.coding_gene.json
```


## EnsembleIDxref

```
python ./gene/preproc_EnsembleIDxref.py
```
**input:**

1. Uniprot Human IDmapping: `[DATASOURCE_PATH]/UNIPROT/HUMAN_9606_idmapping.dat.gz`
 * from: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz)

 ```
$ gzhead HUMAN_9606_idmapping.dat.gz
P31946	UniProtKB-ID	1433B_HUMAN
P31946	Gene_Name	YWHAB
P31946	GI	4507949
P31946	GI	377656702
P31946	GI	67464628
P31946	GI	1345590
P31946	GI	1034625756
P31946	GI	21328448
P31946	GI	377656701
```

1. Uniprot Huamn IDmapping selected tab: `[DATASOURCE_PATH]/UNIPROT/HUMAN_9606_idmapping.dat.gz`
 * from: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz)
1. Ensembl Uniprot ID mapping file: `[DATASOURCE_PATH]/ENSEMBL/hg38/Homo_sapiens.GRCh38.98.uniprot.tsv.gz`
 * from: [ftp://ftp.ensembl.org/pub/release-99/tsv/homo_sapiens/Homo_sapiens.GRCh38.99.uniprot.tsv.gz](ftp://ftp.ensembl.org/pub/release-99/tsv/homo_sapiens/Homo_sapiens.GRCh38.99.uniprot.tsv.gz)

**output:**

* `[DATASOURCE_PATH]/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.tsv.gz`


## EnsembleIDxref_transcriptid
```
python ./gene/preproc_EnsembleIDxref_transcriptid.py
```
**input:**

1. EnsembleIDxref result: `[DATASOURCE_PATH]/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.tsv.gz`

**output:**

`[DATASOURCE_PATH]/ENSEMBL/hg38/idmap_ensembl_uniprot_xref.transcriptid.tsv.gz`




## hg19 coordinates
```
python ./gene/preproc_hg19_coordination.py
```

**input**

1. GRCh37.p13 (last version): `[DATASOURCE_PATH]/ENSEMBL/hg38/Homo_sapiens.GRCh37.75.gtf.gz`
 * from: [ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz](ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)

**output**

1. `[DATASOURCE_PATH]/GENE/ensgid_hg19_coordination.tsv`

 ```
 #chrom_hg19	spos_hg19	epos_hg19	strand_hg19	ensgid
1	11869	14412	+	ENSG00000223972
1	14363	29806	-	ENSG00000227232
1	29554	31109	+	ENSG00000243485
1	34554	36081	-	ENSG00000237613
```