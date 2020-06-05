import sys
import pytest
import test_conf
sys.path.append('..')
from src.mutanno.model.datasource import DataSource
from src.mutanno.model.datastructure import DataSourceListStructure


def test_run_DataSource():
    ds = test_conf.get_ds()
    dss = DataSourceListStructure(ds['fullannot']['ds'])

    dss.__dict__['test'] = "test"

    print(dss.__dict__)
    print(dss.test)
    # print(dss.variant_sources['DBSNP'].__dict__)

    # datsrc = DataSource(dss.variant_sources['DBSNP'])
    # datsrc.set_datastructure()

    # for n1 in [10]:
    #     for r1 in range(1,6):
    #         vcflist = []
    #         vcflist.append("data/test_trio_"+str(n1)+"_"+str(r1)+".vcf")
    #         vcflist.append("data/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf")
    #         vcflist.append("out/test_novocaller_"+str(n1)+"_"+str(r1)+".vcf.gz" + '.fullannot.vcf')
    #         vcflist.append("out/test_trio_"+str(n1)+"_"+str(r1)+".vcf.gz" + '.microannot.vcf')
    #         for vcf in vcflist:
    #             print(vcf)
    #             va = DataSource({'vcf':vcf})
    #             va.validate()
    #         break


if __name__=='__main__':
    
    test_run_DataSource()