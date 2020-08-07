

def run_preproc_clinvar(data):
    from .preproc_clinvar import PreprocClinVar
    pp = PreprocClinVar(data)
    pp.run()
