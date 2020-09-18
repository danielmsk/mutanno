

def run_preproc_vep(data):
    from .preproc_vep import PreprocVEP
    pp = PreprocVEP(data)
    pp.run()


def run_preproc_clinvar(data):
    from .preproc_clinvar import PreprocClinVar
    pp = PreprocClinVar(data)
    pp.run()
