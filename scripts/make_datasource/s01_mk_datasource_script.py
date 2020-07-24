#!/usr/bin/env python
# -*- coding: utf-8 -*-
# split_run_datasource_4_microannot.py
# made by Min-Seok Kwon
# 2020-01-14 12:14:40
#########################
import sys
import os
SVRNAME = os.uname()[1]
if "MBI" in SVRNAME.upper():
    sys_path = "/Users/pcaso/bin/python_lib"
elif SVRNAME == "T7":
    sys_path = "/ms1/bin/python_lib"
else:
    sys_path = "/home/mk446/bin/python_lib"
sys.path.append(sys_path)


s1 = "{:0>"+str(5)+"d}"

def split_run_datasource(ds_json_file, path, chunksize=1000000):
    # seq_util.load_refseq_info('b38d')
    
    file_util.check_dir(path + 'aa')

    regionlist = seq_util.get_split_region(chunksize)
    # print(regionlist[:20])
    # print(regionlist)
    print(len(regionlist))
    mcmd = ""
    rcmd = ""
    i = 0
    cmd = ""
    mchrom_cmd = {}
    rmchrom_cmd = {}

    file_save(path[:-1] + '_splitrun.sh', '', 'w')
    file_save(path[:-1] + '_merge.sh', '', 'w')

    prev_chrom = ''
    for r1 in regionlist:
        i += 1

        chrom = r1[0]
        cmd += 'mutanno makedata '
        # cmd += "-ds /home/mk446/mutanno/SRC/tests/datastructure_microannot_v0.2.json "
        cmd += " -ds " + ds_json_file
        out = path + "mc_" + s1.format(r1[3]) + ".tsi"
        cmd += " -out " + out
        cmd += " -vartype SNV "
        cmd += " -region " + chrom + ":" + str(r1[1]) + "-" + str(r1[2]) + " "

        # cmd += "sleep 5;"
        # tsifile = path + "mc_" + str(r1[3]) + ".tsv.tsi"
        # cmd += "python /home/mk446/mutanno/SRC/scripts/microannotation/s02_mk_check_tsi.py " + tsifile
        # cmd += " " + str(r1[1])
        # cmd += " " + str(r1[2])
        # cmd += ";"
        # print(cmd)
        cmd += "\n"

        file_save(path[:-1] + '_splitrun.sh', cmd, 'a')
        cmd = ""

        if i == 1:
            mcmd += "cat " + out + ".tsi > " + path + "all.tsi;\n"
        else:
            mcmd += "cat " + out + ".tsi | grep -v '^#' >> " + path + "all.tsi;\n"

        file_save(path[:-1] + '_merge.sh', mcmd, 'a')
        mcmd = ''

        out_chrom = "allchrom_" + chrom + ".tsi"
        if prev_chrom != chrom:
            mchrom_cmd = "cat " + out + ".tsi | grep -v '^#' >> " + path + out_chrom + ";\n"
        else:
            mchrom_cmd = "cat " + out + ".tsi > " + path + out_chrom + ";\n"

        rmchrom_cmd = "rm " + out + ".tsi;\n"
        
        outsh = path[:-1] + '_rm_by_chr'+chrom+'.sh'
        file_save(outsh, rmchrom_cmd, 'a')

        prev_chrom = chrom

    # for chrom in mchrom_cmd.keys():
    #     out_chrom = "allchrom_" + chrom + ".tsi"
    #     mchrom_cmd[chrom] += "tabixgz " + path + out_chrom + ";\n"
    #     outsh = path[:-1] + '_merge_by_chr'+chrom+'.sh'
    #     # file_save(outsh, mchrom_cmd[chrom])

    #     outsh = path[:-1] + '_rm_by_chr'+chrom+'.sh'
    #     file_save(outsh, rmchrom_cmd[chrom])
        

def file_save(out, cont, opt = 'a'):
    file_util.fileSave(out , cont, opt)    
    proc_util.run_cmd('chmod 755 ' + out)
    print('Saved', out)

# This function needs lots of computing time. 
def split_run_datasource_4_microannot_chrom(ds_json_file, out):
    # seq_util.load_refseq_info('b38d')
    
    regionlist = seq_util.CHROM_LEN['b38d']
    
    # print(regionlist)
    # print(len(regionlist))
    for chrom in regionlist.keys():
        if len(chrom) < 3:
            out2 =out.replace('#CHROM#',chrom)
            cmd = ""
            cmd = 'mutanno makedata '
            cmd += " -ds " + ds_json_file
            cmd += " -out " + out2 + " "
            cmd += " -vartype SNV "
            cmd += " -region " + chrom + ":1-" + str(regionlist[chrom]) 
            cmd += " -blocksize 5000 "
            cmd += "; tabixgz " + out2 + ";"
            print(cmd)
            # cmd += "sleep 5;"
            # tsifile = path + "mc_" + str(r1[3]) + ".tsv.tsi"
            # cmd += "python /home/mk446/mutanno/SRC/scripts/microannotation/s02_mk_check_tsi.py " + tsifile
            # cmd += " " + str(r1[1])
            # cmd += " " + str(r1[2])
            # cmd += ";"

def split_run_datasource_4_fullannot_chrom(ds_json_file, out):
    # seq_util.load_refseq_info('b38d')
    
    regionlist = seq_util.CHROM_LEN['b38d']
    
    # print(regionlist)
    # print(len(regionlist))
    for chrom in regionlist.keys():
        if len(chrom) < 3:
            out2 =out.replace('#CHROM#',chrom)
            cmd = ""
            cmd = 'mutanno makedata '
            cmd += " -ds " + ds_json_file
            cmd += " -out " + out2 + " "
            cmd += " -vartype SNV "
            cmd += " -region " + chrom + ":1-" + str(regionlist[chrom]) 
            cmd += " -blocksize 50000 "
            cmd += "; tabixgz " + out2 + ";"
            print(cmd)
            # cmd += "sleep 5;"
            # tsifile = path + "mc_" + str(r1[3]) + ".tsv.tsi"
            # cmd += "python /home/mk446/mutanno/SRC/scripts/microannotation/s02_mk_check_tsi.py " + tsifile
            # cmd += " " + str(r1[1])
            # cmd += " " + str(r1[2])
            # cmd += ";"


# -vartype SNV -blocksize 10000 -region 1:1-1000000

if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    print('#USAGE:python s01_mk_datasource_script.py [DS_FILE] [OUT_PATH] [chunksize(default:1000000)]')
    path = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/tmp/"
    ds_json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_microannot_v0.4.2ds.json"
    out = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microannot_datasource.#CHROM#.v0.4.2_200504.tsi"

    # ds_json_file = sys.argv[1]
    # path = sys.argv[2]
    chunksize = 1000000
    if len(sys.argv) == 4:
        chunksize = int(sys.argv[3].replace(',',''))
    # split_run_datasource(ds_json_file, path, chunksize)
    
    ds_json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_microannot_v0.4.3dso2.json"
    out = "/home/mk446/mutanno/DATASOURCE/MICROANNOT/microannot_datasource.#CHROM#.v0.4.4_200614.tsi"
    # split_run_datasource_4_microannot_chrom(ds_json_file, out)
    

    chunksize = 3000
    ds_json_file = "/home/mk446/mutanno/SRC/tests/data/datastructure_fullannot_v0.4.6ds.json"
    path = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/tmp/"
    out = "/home/mk446/mutanno/DATASOURCE/MUTANOANNOT/fullannot_datasource.#CHROM#.v0.4.6_200617.tsi"
    # split_run_datasource(ds_json_file, path, chunksize)
    split_run_datasource_4_fullannot_chrom(ds_json_file, out)
