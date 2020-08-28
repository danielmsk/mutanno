import sys
import pytest
sys.path.append('..')
from src.mutanno.util import vcf_util


def test_convert_to_metadata():

    exlist = [
        (
            {'MUTANNO': {
                'MUTANNO': {'Version': '0.3.19', 'Date': '2020.05.25'}, 
                'VEP': {'Version': 'v99', 'Date': '11/04/2019'}, 
                'gnomADgenome': {'Version': '3.0', 'Date': '03/06/2019'}, 
                'CLINVAR': {'Version': '20200106', 'Date': '01/06/2020'}, 
                'SpliceAI': {'Version': '20191004', 'Date': '10/04/2019'}
            }},
            '''##MUTANNO=<ID=MUTANNO,Version="0.3.19",Date="2020.05.25">
##MUTANNO=<ID=VEP,Version="v99",Date="11/04/2019">
##MUTANNO=<ID=gnomADgenome,Version="3.0",Date="03/06/2019">
##MUTANNO=<ID=CLINVAR,Version="20200106",Date="01/06/2020">
##MUTANNO=<ID=SpliceAI,Version="20191004",Date="10/04/2019">'''
            )
    ]
    for ex in exlist:
        assert vcf_util.convert_to_metadata(ex[0]) == ex[1]
    
if __name__=='__main__':
    test_convert_to_metadata()