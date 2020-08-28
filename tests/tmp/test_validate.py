import sys
import pytest
sys.path.append('..')
from src.mutanno import validate
from src.mutanno.validate import AnnotVCFValidator


def test_pars_header():
    exlist = [
        ('##source=CombineGVCFs',{'source':'CombineGVCFs'}),
        ('##contig=<ID=chrM,length=16569>',
            {'contig':{'ID':'chrM', 'length':'16569'}}),
        ('##FILTER=<ID=LowQual,Description="Low quality TEST=T2">',
            {'FILTER':{'ID':'LowQual', 'Description':'"Low quality TEST=T2"'}}),
        ('##INFO=<ID=VEP,Number=.,Type=String,Description="Predicted nonsense. Format:\'Gene|Feature|Feature_type\'">',
            {'INFO':{'ID':'VEP','Number':'.','Type':'String','Description':'"Predicted nonsense. Format:\'Gene|Feature|Feature_type\'"'}})
    ]

    for ex in exlist:
        assert validate.pars_header(ex[0]) == ex[1]

def test_read_metadata_and_check_duplicate():
    exlist = [
        "##source=CombineGVCFs\n##source=CombineGVCFs",
        "##INFO=<ID=SpliceAI,Number=.,Type=String>\n##INFO=<ID=SpliceAI,Number=.,Type=String>"
    ]

    for ex in exlist:
        with pytest.raises(Exception):
            validate.read_metadata_and_check_duplicate(ex)


    exlist = [
        ("##source=CombineGVCFs\n##source=CombineGVCFs2",{'source':['CombineGVCFs','CombineGVCFs2']}),
        ("##INFO=<ID=SpliceAI,Type=String>\n##INFO=<ID=SpliceAI2,Type=String>",
            {'INFO':{'SpliceAI':{'ID':'SpliceAI','Type':'String'},'SpliceAI2':{'ID':'SpliceAI2','Type':'String'}}})
    ]
    for ex in exlist:
        assert validate.read_metadata_and_check_duplicate(ex[0]) == ex[1]


def test_run_AnnotVCFValidator():
    for n1 in [10]:
        for r1 in range(1,6):
            vcflist = []
            # vcflist.append("data/test_trio_"+str(n1)+"_"+str(r1)+".vcf")
            # vcflist.append("data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf")
            vcflist.append("out/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf" + '.microannot.vcf.fullannot.vcf')
            # vcflist.append("out/test_trio_"+str(n1)+"_"+str(r1)+".vcf" + '.microannot.vcf')

            for vcf in vcflist:
                print(vcf)
                va = AnnotVCFValidator({'vcf':vcf})
                va.validate()
                break

            break


if __name__=='__main__':
    test_pars_header()
    test_read_metadata_and_check_duplicate()
    test_run_AnnotVCFValidator()