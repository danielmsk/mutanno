Annotated VCF file format
=========================

All annotations are inserted in INFO field of VCF file. 


VCF-Specific Restrictions
-------------------------

For the annotated VCF we make use of INFO fields to encapsulate our annotations. This field is part of the VCF structure and has the following restrictions on values within the field (ie: ‘AC=2;VEP=1|2…’ etc).

1. String format 

	a. Conversion to type specified on the mapping table is done later

2. No whitespace

	a. No tabs, spaces or otherwise

3. No semicolon

	a. Semicolon delineates annotation fields within the INFO block

	b. Marks the end

4. No equals (=)

	a. Equals delineates annotation fields within the INFO block in tandem with semicolon

	b. Marks the beginning (ie: ’AC=2;VEP=1|2…’)

5. Commas are only used to separate lists of values

	a. In our case, this applies to VEP annotations which are both multi-valued and can have multiple entries 
	
	b. ie: ‘VEP=1|2,3|2;’ would mean there are TWO VEP annotations with field values 1|2 and 3|2

Requirements
------------

1. Annotation fields that should be processed as such must be marked with a MUTANNO tag in the VCF metadata as below. ##MUTANNO=<ID=VEP,Version="v97.4",Date="07/04/2019">

2. Annotation fields that have MUTANNO tags must also have a corresponding INFO tag. This tag must specify a format if the annotation is multi-valued and must be pipe (|) separated. An example of each is below.

	a. ##INFO=<ID=1000GP,Number=.,Type=String,Description="population AF by 1000GP.

	b. ##INFO=<ID=SNPEFFLOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant by SNPEFF. Format:'Gene_Name|Gene_ID|Number_of_transcripts_in_gene|Percent_of_transcripts_affected' ">


3. If an annotation field can have multiple entries, as is the case with VEP, these entries must be comma separated as consistent with the VCF requirements

	a. VEP=1%3A65565|G|CCDS30547.1|CCDS30547.1|Transcript|upstream_gene_variant|||||||MODIFIER|3526|1||SNV||||protein_coding|YES||||CCDS30547.1|CCDS30547.1|||||||||||||||||||||,1%3A65565|G|ENSG00000186092|ENST00000335137|Transcript|upstream_gene_variant|||||||MODIFIER|3490|1||SNV|OR4F5|HGNC|HGNC%3A14825|protein_coding|YES|||P1|CCDS30547.1|ENSP00000334393|Q8NH21||UPI0000041BC1|||||||||||||||||| … 

4. If an annotation field within a sub-embedded object is an array, such as vep_domains, those entries must be tilde (~) separated

	a. VEP= … |val_1~val_2~val_3| … → process field as [val_1, val_2, val_3]

5. In addition, if an annotation field is a sub-embedded object, the DESCRIPTION field in the INFO tag must have ‘Subembedded=’<val>’ where ‘val’ is the name of the sub_embedded_grouping on the mapping table

6. No further nesting is allowed

	a. You cannot have an annotation sub-field that is an array of arrays and so forth

Separator Summary
-----------------

* Tab separates VCF specific fields and is thus restricted
* Semicolon separates different annotation fields within INFO and is thus restricted
* Comma separates sub-embedded objects within a single INFO field (such as VEP) and cannot be used in any other way
* Pipe separates multi-valued annotation fields and cannot be used in any other way
* Tilde separates sub-embedded objects that are also arrays, such as vep_domain and cannot be used in any other way

Example
-------

Given these restrictions, below is a detailed walk through of how the VCF parses the annotation fields given this specification. A truncated example entry is below. Assume we are able to grab appropriate MUTTANO/INFO header information. New lines are inserted for readability but are not present in the actual file.

.. code::

	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002
	chr1	65565	.	A	G	58.56	VQSRTrancheSNP99.00to99.90	

* Bold above represents the VCF data header. Fields other than INFO are readily accessible. All annotation fields are collapsed into the INFO section. FORMAT and HG002 follow after INFO.
* The fields below are tab separated as consistent with the VCF specification. A tab separates the last part of the data above and the INFO data below.

.. code::

	AC=2;AF=0.500;AN=4;DP=24;ExcessHet=0.7918;FS=0.000;MLEAC=2;MLEAF=0.500;MQ=65.65;NEGATIVE_TRAIN_SITE;QD=29.28;SOR=2.303;VQSLOD=-3.874e+00;culprit=DP;


* These annotations are all single valued and are thus processed directly as strings. Conversion to actual types is done later 

.. code::

	VEP=1%3A65565|G|CCDS30547.1|CCDS30547.1|Transcript|upstream_gene_variant|||||||MODIFIER|3526|1||SNV||||protein_coding|YES||||CCDS30547.1|CCDS30547.1|||||||||||||||||||||,

	1%3A65565|G|ENSG00000186092|ENST00000335137|Transcript|upstream_gene_variant|||||||MODIFIER|3490|1||SNV|OR4F5|HGNC|HGNC%3A14825|protein_coding|YES|||P1|CCDS30547.1|ENSP00000334393|Q8NH21||UPI0000041BC1||||||||||||||||||,

	1%3A65565|G|ENSG00000240361|ENST00000492842|Transcript|downstream_gene_variant|||||||MODIFIER|1678|1||SNV|OR4G11P|HGNC|HGNC%3A31276|transcribed_unprocessed_pseudogene|||||||||||||||||||||||||||;


* Above is a VEP annotation entry that is both multi-valued and has multiple entries
* To parse this we first split on the comma to get the groups. Newlines are inserted to visualize the groups.
* We then split on pipe since the fields are pipe separated. Even if a field is blank a pipe must be present for that field otherwise we will not be able to determine which fields go with which values
* Once we have all the fields, we then go through each one and post-process. If it is an array field (not shown in this example but consistent with point 4 above) then we split again on tilde to determine the array elements, otherwise the field value is cast to the appropriate type.

Item Generation
---------------

1. Once we have processed the VCF a dictionary is created that roughly represents the structure of each VCF record (one per line)

	a. Keys are annotation fields, values are either direct or keyed again (sub-dictionary) on the subfield

	b. Ex: { ‘CHROM’ : ‘chr1’, ‘POS’: 65565 … ‘VEP’ : { ‘Location’ : <val>, ‘Allele’: <val> … } }

2. This dictionary is then converted to the format expected by Elasticsearch (TBD)

	a. The above described dictionary format should thus be considered temporary
