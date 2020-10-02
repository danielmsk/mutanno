

def run_preproc_vep(data):
    from .preproc_vep import PreprocVEP
    pp = PreprocVEP(data)
    pp.run()


def run_preproc_clinvar(data):
    from .preproc_clinvar import PreprocClinVar
    pp = PreprocClinVar(data)
    pp.run()


def run_preproc_gnomad(data):
    from .preproc_gnomad import PreprocGnomAD
    pp = PreprocGnomAD(data)
    pp.run()


def run_preproc_spliceai(data):
    from .preproc_spliceai import PreprocSpliceAI
    pp = PreprocSpliceAI(data)
    pp.run()


def run_preproc_cytoband(data):
    from .preproc_cytoband import PreprocCytoband
    pp = PreprocCytoband(data)
    pp.run()

def run_preproc_ucsc_repeat(data):
    from .preproc_ucsc_repeat import PreprocUCSCRepeat
    pp = PreprocUCSCRepeat(data)
    pp.run()

def run_preproc_uniprot_transmem(data):
    from .preproc_uniprot_transmem import PreprocUniprotTransmem
    pp = PreprocUniprotTransmem(data)
    pp.run()

def run_preproc_dbsnp(data):
    from .preproc_dbsnp import PreprocDbSNP
    pp = PreprocDbSNP(data)
    pp.run()

def run_preproc_cadd(data):
    from .preproc_cadd import PreprocCADD 
    pp = PreprocCADD(data)
    pp.run()

def run_preproc_topmed(data):
    from .preproc_topmed import PreprocTopmed 
    pp = PreprocTopmed(data)
    pp.run()