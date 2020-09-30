Version History
===============

v0.4.x release series
---------------------
0.4.8 (2020.09.30) :
	- check redundancy of VEP fields in full-annotation

0.4.7 (2020.09.24) :
	- add preprocessing functions

0.4.7 (2020.09.23) :
	- debug VEP annotation missing error of some multiallele variants in preproc module

0.4.6 (2020.09.18) :
	- add multiple data source files.

0.4.5 (2020.08.28) :
	- bug fix in download module. add ``-silence`` option

0.4.4 (2020.08.25) :
	- bug fix missing genotype for multiallelic variant when ``-split_multi_allelic_variant`` option is applied. (`issue #5 <https://github.com/dbmi-bgm/mutanno/issues/5>`_)

0.4.3 (2020.08.09) :
	- add ``-fixpl`` option when converting multi-allelic variant to biallelic variants. 
	- change for full-annotation data source v0.4.8
	- update value mapping order in external function.
	- add configration jsonfile (./data/conf.json)
	- add download sub-command
	- add tqdm library for downloading progress bar
0.4.2 (2020.06.25) :
	- add parameter set in data structure json file.
0.4.1 (2020.06.08) :
	- turn to be 'single_source_mode = True' when data structure file has only single source file.
	- allocate all sample names when users use -genoinfo without indicating sample names.
0.4.0 (2020.06.05) :
	- add variant type
	- support multiple datasources for vcf annotation
	- remove block-size option
	- make `MUTANNO` source name for 'samplevariantkey', 'hgvsg', and 'variant_class'
	- add '-variant_class' option
	- refactoring
	- add multi-source mode
	- refactoring makedata mode	


v0.3.x release series
---------------------

0.3.19 (2020.05.25) :
	- split pred and score in pathogenicity score from VEP (SIFT and Polyphen)
	- add 'region_vcf' option in 'makedata'
0.3.18 (2020.05.13) :
	- add chain file option for pyliftover
0.3.17 (2020.05.06) :
	- add is_most_severe_transcript in external function.py
	- add keep_equal_str filter in data structure file
	- update -geno_info option to get sample ID list (ex. -geno_info HG002 HG003 HG004 )
	- update hgvs chromosome.
0.3.16 (2020.04.27) :
	- add GENE_MAIN_CHROM, CODING_GENE_MAIN_CHROM in -vartype option
0.3.15 (2020.03.26) :
	- add 'web' option
0.3.14 (2020.03.21) :
	- implement pypi installation
	- apply pytest
0.3.13 (2020.03.10) :
	- add 'add_genoinfo' option
	- add 'split_multi_allelic_variant' option
	- add 'clean_tag' option
0.3.12 (2020.03.05) :
	- add 'default' value
0.3.11 (2020.02.28) :
	- add SAMPLEGENO tag
	- minor fix in header
0.3.10 (2020.02.26) :
	- minor update in version tag
0.3.9 (2020.02.26) :
	- remove unannotated variant (-remove_unannotated_variant)
	- update the format of multiallele tag
0.3.8 (2020.02.25)
	- debugged first variant missing in annot module
0.3.7
	- update annot module
0.3.6
	- add makedata	
0.3.5
	- merge some dbNSFP fields into transcript table (dbNSFPTranscript)
0.3.4
	- change type of `dbNSFP SiPhy_29way_pi` to list
0.3.3
	- encode 'space' to '%20' (remove blank space)









