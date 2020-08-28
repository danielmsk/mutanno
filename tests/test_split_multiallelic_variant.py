import sys
sys.path.append('..')
from src.mutanno.util import vcf_util

# from mutanno.util import vcf_util


def get_gt(rec):
    return rec.split(':')[0]

def get_ad(rec):
    try:
        rst = rec.split(':')[1]
    except IndexError:
        rst = ""
    return rst

def get_dp(rec):
    try:
        rst = rec.split(':')[2]
    except IndexError:
        rst = ""
    return rst

def get_pl(rec):
    try:
        rst = rec.split(':')[6]
    except IndexError:
        rst = ""
    return rst


def get_ac(rec):
    rst = ""
    for r1 in rec.split(';'):
        if "AC=" in r1:
            rst = r1.replace('AC=','')
            break
    return rst


def get_af(rec):
    rst = ""
    for r1 in rec.split(';'):
        if "AF=" in r1:
            rst = r1.replace('AF=','')
            break
    return rst


vcflines = []
vcflines.append("""
chr9	134031310	.	T	TGGGGGGG,TGGG	268.18	.
	AC=1,1;AF=0.167,0.167;AN=6;BaseQRankSum=-1.213e+00;DP=100;ExcessHet=3.01;FS=113.147;MLEAC=1,1;MLEAF=0.167,0.167;MQ=60.00;MQRankSum=0.00;QD=5.16;ReadPosRankSum=-2.791e+00;SOR=5.130
	GT:AD:DP:GQ:PGT:PID:PL:PS
	0/0:43,0,0:43:63:.:.:0,63,945,63,945,945	0|1:20,8,0:28:99:0|1:134031310_T_TGGGGGGG:238,0,2316,305,2340,2644:134031310	0|2:20,0,4:24:47:0|1:134031310_T_TGGG:47,111,1727,0,1617,1605:134031310
""")
vcflines.append("""
chr1	4420293	.	ATGTGTG	*,A	188.97	.
	AC=3,1;AF=0.750,0.250;AN=4;DP=31;ExcessHet=3.01;FS=0.000;MLEAC=3,1;MLEAF=0.750,0.250;MQ=52.57;QD=8.22;SOR=1.179;SAMPLEGENO=1/1|*/*|0/4/0|NA12877_sample,1/2|*/A|0/14/5|NA12878_sample,./.|./.|0/0/0|NA12879_sample
	GT:AD:DP:GQ:PGT:PID:PL:PS
	1|1:0,4,0:4:12:1|1:4420279_GTATATATATATATA_G:129,12,0,129,12,129:4420279
	1/2:0,14,5:19:99:.:.:759,194,165,555,0,551
	./.:0,0,0
""")
vcflines.append("""
chr3	94029349	.	ATTT	TTTT,A	1376.17	.
	AC=3,1;AF=0.750,0.250;AN=4;DP=40;ExcessHet=3.01;FS=0.000;MLEAC=3,1;MLEAF=0.750,0.250;MQ=57.93;QD=32.29;SOR=1.112
	GT:AD:DP:GQ:PGT:PID:PL:PS
	1/2:0,9,17:28:99:.:.:1276,929,901,349,0,290
	1/1:0,4,0:4:12:.:.:123,12,0,127,15,141
	.|.:0,0,0:.:.:1|0:94029301_CATAT_C:.:94029301
""")


rst = []
d = {0: {}, 1:{}}
d[0][9] =  {'gt': '0/0', 'ad': '43,0', 'dp': '43', 'pl': '0,63,945'}
d[0][10] = {'gt': '0|1', 'ad': '20,8', 'dp': '28', 'pl': '238,0,2316'}
d[0][11] = {'gt': '0|0', 'ad': '20,0', 'dp': '20', 'pl': '47,111,1727'}
d[0]['ac'] = "1"
d[0]['af'] = "0.167"
d[1][9] =  {'gt': '0/0', 'ad': '43,0', 'dp': '43', 'pl': '0,63,945'}
d[1][10] = {'gt': '0|0', 'ad': '20,0', 'dp': '20', 'pl': '238,305,2644'}
d[1][11] = {'gt': '0|1', 'ad': '20,4', 'dp': '24', 'pl': '47,0,1605'}
d[1]['ac'] = "1"
d[1]['af'] = "0.167"
rst.append(d)
d = {0: {}, 1: {}}
d[0][9] = {'gt': '1|1', 'ad': '0,4', 'dp': '4', 'pl': '129,12,0'}
d[0][10] = {'gt': '1/1', 'ad': '0,14', 'dp': '14', 'pl': '759,194,165'}
d[0][11] = {'gt': './.', 'ad': '0,0', 'dp': '', 'pl': ''}
d[0]['ac'] = "4"
d[0]['af'] = "0.667"
d[1][9] = {'gt': '0|0', 'ad': '0,0', 'dp': '0', 'pl': '129,129,129'}
d[1][10] = {'gt': '1/1', 'ad': '0,5', 'dp': '5', 'pl': '759,555,551'}
d[1][11] = {'gt': './.', 'ad': '0,0', 'dp': '', 'pl': ''}
d[1]['ac'] = "2"
d[1]['af'] = "0.333"
rst.append(d)
d = {0: {}, 1: {}}
d[0][9] = {'gt': '1/1', 'ad': '0,9', 'dp': '9', 'pl': '1276,929,901'}
d[0][10] = {'gt': '1/1', 'ad': '0,4', 'dp': '4', 'pl': '123,12,0'}
d[0][11] = {'gt': '.|.', 'ad': '0,0', 'dp': '0', 'pl': '.'}
d[0]['ac'] = "4"
d[0]['af'] = "0.667"
d[1][9] = {'gt': '1/1', 'ad': '0,17', 'dp': '17', 'pl': '1276,349,290'}
d[1][10] = {'gt': '0/0', 'ad': '0,0', 'dp': '0', 'pl': '123,127,141'}
d[1][11] = {'gt': '.|.', 'ad': '0,0', 'dp': '0', 'pl': '.'}
d[1]['ac'] = "2"
d[1]['af'] = "0.333"
rst.append(d)


def test_split_multiallelic_variants():
    for k, vcfline in enumerate(vcflines):
        vcfrecord = vcfline.strip().split('\t')
        rstrecord = vcf_util.split_multiallelic_variants(vcfrecord)

        arralt = vcfrecord[4].split(',')
        assert len(rstrecord) == len(vcfrecord[4].split(',')), "Now matched number of split records."
        
        for i in range(len(rstrecord)):
            for j in [0,1]: # alt index
                if rstrecord[i][4] == arralt[j]:
                    for l in [9,10,11]:
                        assert get_gt(rstrecord[i][l]) == rst[k][j][l]['gt'], get_ad(rstrecord[i][l]) + " != " + rst[k][j][l]['gt']
                        assert get_ad(rstrecord[i][l]) == rst[k][j][l]['ad'], get_ad(rstrecord[i][l]) + " != " + rst[k][j][l]['ad']
                        assert get_dp(rstrecord[i][l]) == rst[k][j][l]['dp'], get_ad(rstrecord[i][l]) + " != " + rst[k][j][l]['dp']
                        assert get_pl(rstrecord[i][l]) == rst[k][j][l]['pl'], get_ad(rstrecord[i][l]) + " != " + rst[k][j][l]['pl']

            assert get_ac(rstrecord[i][7]) == rst[k][i]['ac'], get_ac(rstrecord[i][7]) + " != " + rst[k][i]['ac']
            assert get_af(rstrecord[i][7]) == rst[k][i]['af'], get_af(rstrecord[i][7]) + " != " + rst[k][i]['af']

if __name__ == "__main__":
    test_split_multiallelic_variants()
    pass
