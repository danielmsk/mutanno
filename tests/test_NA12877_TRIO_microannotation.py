import sys
import shlex
sys.path.append('..')
from src import mutanno

prog = "mutanno"

cmdlist = []
# cmdlist.append("""
#     annot \
#     -vcf /Users/pcaso/db/MUTANNO/TESTDATA/NA12877_TRIO_GAPFIS8ZSPEO.vcf.gz \
#     -ds ./data_structure_json/datastructure_microannot_v0.4.5.json \
#     -out output.annot.vcf \
#     -sourcefile /Users/pcaso/db/MUTANNO/DATASOURCE/MICROANNOT/v0.4.4_200614/microannot_datasource.v0.4.4_200614.tsi.gz \
#     -split_multi_allelic_variant \
#     -genoinfo \
#     -single_source_mode
# """)

#----------------------------------
# chrom20 only
#----------------------------------
# cmdlist.append("""
#     annot \
#     -vcf /Users/pcaso/db/MUTANNO/TESTDATA/NA12877_TRIO_GAPFIS8ZSPEO_CHR20.vcf.gz \
#     -ds ./data_structure_json/datastructure_microannot_v0.4.5.json \
#     -out output.annot.vcf \
#     -sourcefile /Users/pcaso/db/MUTANNO/DATASOURCE/MICROANNOT/v0.4.4_200614_chrom/microannot_datasource.20.v0.4.4_200614.tsi.gz \
#     -split_multi_allelic_variant \
#     -genoinfo \
#     -single_source_mode
# """)

#----------------------------------
# multiallelic only
#----------------------------------
cmdlist.append("""
    annot \
    -vcf /Users/pcaso/db/MUTANNO/TESTDATA/NA12877_TRIO_GAPFIS8ZSPEO.multiallelic.vcf.gz \
    -ds ./data_structure_json/datastructure_microannot_v0.4.5.json \
    -out output.annot.vcf \
    -sourcefile /Users/pcaso/db/MUTANNO/DATASOURCE/MICROANNOT/v0.4.4_200614_chrom/microannot_datasource.20.v0.4.4_200614.tsi.gz \
    -split_multi_allelic_variant \
    -genoinfo \
    -single_source_mode
""")



def test_run_download():
    for cmd in cmdlist:
        cmd = prog + " " + cmd.strip()
        sys.argv = shlex.split(cmd)
        print(' '.join(sys.argv))
        # print(cmd)
        # print(shlex.quote(sys.argv))
        mutanno.cli()


if __name__ == "__main__":
    test_run_download()
