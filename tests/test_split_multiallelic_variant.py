import sys
sys.path.append('..')
from src.mutanno.util import vcf_util

# from mutanno.util import vcf_util


def get_gt(rec):
    return rec.split(':')[0]

def get_ad(rec):
    return rec.split(':')[1]

def get_dp(rec):
    return rec.split(':')[2]

def get_pl(rec):
    return rec.split(':')[6]

def test_split_multiallelic_variants():
    vcfline = "chr9	134031310	.	T	TGGGGGGG,TGGG	268.18	."
    vcfline += "	AC=1,1;AF=0.167,0.167;AN=6;BaseQRankSum=-1.213e+00;DP=100;ExcessHet=3.01;FS=113.147;MLEAC=1,1;MLEAF=0.167,0.167;MQ=60.00;MQRankSum=0.00;QD=5.16;ReadPosRankSum=-2.791e+00;SOR=5.130"
    vcfline += "	GT:AD:DP:GQ:PGT:PID:PL:PS"
    vcfline += "	0/0:43,0,0:43:63:.:.:0,63,945,63,945,945	0|1:20,8,0:28:99:0|1:134031310_T_TGGGGGGG:238,0,2316,305,2340,2644:134031310	0|2:20,0,4:24:47:0|1:134031310_T_TGGG:47,111,1727,0,1617,1605:134031310"
    vcfrecord = vcfline.split('\t')
    rstrecord = vcf_util.split_multiallelic_variants(vcfrecord)

    assert len(rstrecord) == len(vcfrecord[4].split(',')), "Now matched number of split records."
    for i in range(len(rstrecord)):
        assert get_gt(rstrecord[i][9]) == "0/0" #first genotype
        assert get_ad(rstrecord[i][9]) == "43,0" #first allele depth
        assert get_dp(rstrecord[i][9]) == "43" #third total depth
        assert get_pl(rstrecord[i][9]) == "0,63,945" #third total depth

        if rstrecord[i][4] == "TGGGGGGG":
            assert get_gt(rstrecord[i][10]) == "0|1" #second genotype
            assert get_ad(rstrecord[i][10]) == "20,8" #second allele depth
            assert get_dp(rstrecord[i][10]) == "28" #second total depth
            assert get_pl(rstrecord[i][10]) == "238,0,2316" #second total depth

            assert get_gt(rstrecord[i][11]) == "0|0" #third genotype
            assert get_ad(rstrecord[i][11]) == "20,0" #third allele depth
            assert get_dp(rstrecord[i][11]) == "20" #third total depth
            assert get_pl(rstrecord[i][11]) == "47,111,1727" #third total depth
        else:
            assert get_gt(rstrecord[i][10]) == "0|0" #second genotype
            assert get_ad(rstrecord[i][10]) == "20,0" #second allele depth
            assert get_dp(rstrecord[i][10]) == "20" #second total depth
            assert get_pl(rstrecord[i][10]) == "238,305,2644" #second total depth

            assert get_gt(rstrecord[i][11]) == "0|1" #third genotype
            assert get_ad(rstrecord[i][11]) == "20,4" #third allele depth
            assert get_dp(rstrecord[i][11]) == "24" #third total depth
            assert get_pl(rstrecord[i][11]) == "47,0,1605" #third total depth
    

if __name__ == "__main__":
    test_split_multiallelic_variants()
    pass
