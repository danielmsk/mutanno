<!-- ##TOC## -->
1. [Variant Annotation](#variant-annotation)
	1. [Make Annotation Source Files](#make-annotation-source-files)
		1. [Make micro-annotation source files](#make-micro-annotation-source-files)
		1. [Make main-annotation source files](#make-main-annotation-source-files)
		1. [Compare json data structure file with googlesheet table](#compare-json-data-structure-file-with-googlesheet-table)
	1. [CLINVAR](#clinvar)
	1. [Conservation](#conservation)
		1. [Merge Scores for](#merge-scores-for)
		1. [SiPhy](#siphy)
	1. [LIFTOVER](#liftover)
1. [Gene Annotation](#gene-annotation)
	1. [Make Datasource file](#make-datasource-file)
	1. [EnsembleIDxref](#ensembleidxref)
	1. [EnsembleIDxref_transcriptid](#ensembleidxref_transcriptid)
	1. [hg19 coordinates](#hg19-coordinates)
<!-- ##TOCEND## -->

# Variant Annotation


## Make Annotation Source Files

```
python s01_mk_datasource_script.py [DS_FILE] [OUT_PATH] [chunksize(default:1000000)]
python s04_merge_and_tabixgz.py > s04.sh
./s04.sh
```

### Make micro-annotation source files
```
python s01_mk_datasource_script.py \
	../tests/data/datastructure_microannot_v0.4.2ds.json \
	../../DATASOURCE/MICROANNOT/tmp/ \
	1000000
```

### Make main-annotation source files
```
python s01_mk_datasource_script.py \
	../tests/data/datastructure_v0.4.4ds_mvp.json \
	../../DATASOURCE/MAINANNOT/tmp/ \
	500000
```


### Compare json data structure file with googlesheet table
```
python maptable/check_variant_json_vs_googlesheet.py \
	../tests/data/datastructure_v0.4.4ds.json \
	../tests/data/datastructure_v0.4.4ds.googlesheet.txt
```
## CLINVAR
```
python s02_preproc_clinvar.py 

vcf-sort -c clinvar_20200329_submission.tsv > clinvar_20200329_submission.sorted.tsv
tabixgz clinvar_20200329_submission.sorted.tsv.gz
```

*Input:*

* `clinvarxml = path + "ClinVarFullRelease_2020-03.xml.gz"`
* `clinvarvcf = path + "clinvar_20200329.vcf.gz"`
* download from [ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar)

*Output:*

* tsv file : `clinvar_20200329.tsv`

	```
[0]	#CHROM :	1
[1]	POS :	930248
[2]	VARIATIONID :	789256
[3]	REF :	G
[4]	ALT :	A
[5]	ALLELEID :	707587
[6]	CLNDN :	not_provided
[7]	CLNDNINCL :
[8]	CLNDISDB :	MedGen:CN517202
[9]	CLNDISDBINCL :
[10]	CLNHGVS :	NC_000001.11:g.930248G>A
[11]	CLNREVSTAT :	criteria_provided,_single_submitter
[12]	CLNSIG :	Likely_benign
[13]	CLNSIGCONF :
[14]	CLNSIGINCL :
[15]	CLNVC :	single_nucleotide_variant
[16]	CLNVCSO :	SO:0001483
[17]	CLNVI :
[18]	GENEINFO :	SAMD11:148398
[19]	MC :	SO:0001583|missense_variant
[20]	ORIGIN :	1
[21]	SSR :
```

* submission tsv file: `clinvar_20200329_submission.tsv`

	```
[0]	#CHROM :	1
[1]	POS :	930248
[2]	VARIATIONID :	789256
[3]	REF :	G
[4]	ALT :	A
[5]	ClinVarAccession :	SCV001119512
[6]	Interpretation :	Likely benign
[7]	DateLastEvaluated :	2019-01-02
[8]	ReviewStatus :	criteria provided, single submitter
[9]	Method :	clinical testing
[10]	Condition :	not provided
[11]	AlleleOrigin :	germline
[12]	Submitter :	Invitae
[13]	SubmitterID :	500031
[14]	Citation :
```


## Conservation

### Merge Scores for 

```
python preproc_convert_and_merge_conservation_wigfix2bed.py
```

**Input:**

1. Wigfile for PhastCons, PhyloP, Gerp
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS100WAY/hg38/chr#CHROM#.phastCons100way.wigFix.gz`
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS30WAY/hg38/chr#CHROM#.phastCons30way.wigFix.gz`
 * `[DATASOURCE_PATH]/CONSERVATION/PHASTCONS/PHASTCONS20WAY/hg38/chr#CHROM#.phastCons20way.wigFix.gz`

2. SiPhy bed file (tabix-indexed)
 * 

**Output:**

1. `[DATASOURCE_PATH]/CONSERVATION/conservation_scores.hg38.chr##CHROM##.tsi.gz`

 ```
#CHROM	POS	ID	PHASTCONS100|PHASTCONS30|PHASTCONS20|PHYLOP100|PHYLOP30|PHYLOP20|GERP|SIPHY_OMEGA|SIPHY_P
5	11014		|0.113|||0.146||||
5	11015		|0.106|||0.158||||
5	11016		|0.099|||0.177||||
```


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


## LIFTOVER

```
python ./liftover/liftover.py [b37 vcf file]
```
**Input**

* b37 VCF file

**Output**

* hg38-liftover annotated VCF file




# Gene Annotation

## Make Datasource file

**For coding genes in main chromosomes (1-22,X,Y,MT)**

```
mutanno makedata -ds ./tests/data/datastructure_v0.4.4ds_mvp.json \
    -out ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.4.coding_gene_main_chrom \
    -vartype CODING_GENE_MAIN_CHROM \
    -outtype json
gz ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.4.coding_gene_main_chrom.json
pytest scripts/gene/test_googlesheet_gene_table.py
```


**For all genes**

```
mutanno makedata -ds ./tests/data/datastructure_v0.4.3ds_mvp.json \
    -out ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3 \
    -vartype GENE \
    -outtype json
gz ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3.json
python ./gene/check_googlesheet_gene_table.py 
```

**For coding genes**

```
mutanno makedata -ds ./tests/data/datastructure_v0.4.3ds_mvp.json \
    -out ../DATASOURCE/MUTANOANNOT/mvp_gene_datasource_v0.4.3.coding_gene \
    -vartype CODING_GENE \
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