import sys
import os
import tabix
import time
from mutanno.util import file_util
from mutanno.util import proc_util


def mk_target_sourcefile(ori_source_file, vcf_file):
    out = ori_source_file + '.target.mti'
    f = open(out, 'w')

    # save header
    for line in file_util.gzopen(ori_source_file):
        line = file_util.decodeb(line)
        if line[0] == '#':
            f.write(line)
            break
        
    tb = tabix.open(ori_source_file)
    for line in file_util.gzopen(vcf_file):
        line = file_util.decodeb(line)
        if line[0] != '#':
            arr = line.split('\t')
            chrom = arr[0].replace('chr', '')
            pos = int(arr[1])
            alt = arr[4]
            spos = pos - 1
            epos = pos + len(alt)
            varkey = '_'.join(arr[:2])
            cnt = 0
            recs = tb.querys(chrom + ":"+str(spos)+"-" + str(epos))
            for r1 in recs:
                f.write('\t'.join(r1)+'\n')
                cnt += 1

            print(varkey, cnt)
    f.close()

    time.sleep(1)

    cmd = "vcf-sort -c " + out + ' > ' + out + '.sorted.mti'
    proc_util.run_cmd(cmd, True)

    time.sleep(1)

    prev_line = ""
    f = open(out + '.sorted.uniq.mti', 'w')
    for line in open(out + '.sorted.mti'):
        if line != prev_line:
            f.write(line)
        prev_line = line
    f.close()

    time.sleep(1)

    cmd = "tabixgz " + out + '.sorted.uniq.mti'
    proc_util.run_cmd(cmd, True)

    print('Final source file:', out + '.sorted.uniq.mti.gz')


if __name__ == "__main__":
    print('#USAGE: python mk_target_sorucefile.py [Origianl source file] [VCF]')
    ori_source_file = sys.argv[1]
    vcf_file = sys.argv[2]
    mk_target_sourcefile(ori_source_file, vcf_file)
