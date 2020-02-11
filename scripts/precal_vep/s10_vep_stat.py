#!/usr/bin/env python
# -*- coding: utf-8 -*-
# s10_vep_stat.py
# made by Daniel Minseok Kwon
# 2020-02-05 11:55:01
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
import tabix

def s10_vep_stat_splitrun(chrom, spos, epos, k):
    cnt = {}
    cnt_each = {}
    cnt_merge = {}
    vep = vepfile.replace('#CHROM#',chrom)
    print(vep)
    tb = tabix.open(vep)
    print(vep)
    i = 0
    recs = tb.query(chrom, spos, epos)
    for arr in recs:
        i += 1
        i_cnt_tags = {}
        i_cnt_tag = {}
        i_total_tags = 0
        i_total_tag = 0
        for transcriptinfo in arr[5].split(','):
            arr2 = transcriptinfo.split('|')
            tags = arr2[1]
            try:
                i_cnt_tags[tags] += 1
            except KeyError:
                i_cnt_tags[tags] = 1
            i_total_tags += 1
            
            for tag in tags.split('&'):
                try:
                    i_cnt_tag[tag] += 1
                except KeyError:
                    i_cnt_tag[tag] = 1
                i_total_tag += 1

        for tags in i_cnt_tags.keys():
            try:
                cnt[tags] += i_cnt_tags[tags] / i_total_tags
            except KeyError:
                cnt[tags] = i_cnt_tags[tags] / i_total_tags

        for tag in i_cnt_tag.keys():
            try:
                cnt_each[tag] += i_cnt_tag[tag] / i_total_tag
            except KeyError:
                cnt_each[tag] = i_cnt_tag[tag] / i_total_tag
        
        arrtag = list(i_cnt_tag.keys())
        mtag = '&'.join(sorted(arrtag))

        try:
            cnt_merge[mtag] += 1
        except KeyError:
            cnt_merge[mtag] = 1

        if i % 10000 == 0:
            print(i, chrom, arr[1], len(cnt.keys()))
            # break
    cnt['snv'] = i
    cnt_each['snv'] = i
    cnt_merge['snv'] = i
    print(cnt)
    print(cnt_each)
    print(cnt_merge)

    file_util.jsonSave(vep + '.' + k + '.stat_tag.json', cnt_each)
    file_util.jsonSave(vep + '.' + k + '.stat_tags.json', cnt)
    file_util.jsonSave(vep + '.' + k + '.stat_mergedtag.json', cnt_merge)
    print('Saved', vep + '.' + k + '.stat_tag.json')


def s10_vep_stat(chrom):
    cnt = {}
    cnt_each = {}

    cnt = {}
    cnt_each = {}
    vep = vepfile.replace('#CHROM#',chrom)
    print(vep)
    i = 0
    for line in file_util.gzopen(vep):
        line = line.decode('UTF-8')
        if line[0] != '#':
            arr = line.split('\t')
            i += 1
            for transcriptinfo in arr[5].split(','):
                arr2 = transcriptinfo.split('|')
                tags = arr2[1]
                try:
                    cnt[tags] += 1
                except KeyError:
                    cnt[tags] = 1

                
                for tag in tags.split('&'):
                    try:
                        cnt_each[tag] += 1
                    except KeyError:
                        cnt_each[tag] = 1
                        
            if i % 10000 == 0:
                print(i, chrom, arr[1], len(cnt.keys()))
                # break
    print(cnt)
    print(cnt_each)

    file_util.jsonSave(vep + '.stat_tag.json', cnt_each)
    file_util.jsonSave(vep + '.stat_tags.json', cnt)


def save_stat(jsontype = 'chrom'):
    cnt_each = {}
    cnt_merge = {}
    cnt = {}
    mtagsmap = {}
    tagsmap = {}
    tagmap = {}

    for chrom in seq_util.MAIN_CHROM_LIST:
        vep = vepfile.replace('#CHROM#',chrom)
        if jsontype == 'chrom':
            cnt_merge[chrom] = file_util.jsonOpen(vep + '.stat_mergedtag.json')
            cnt_each[chrom] = file_util.jsonOpen(vep + '.stat_tag.json')
            cnt[chrom] = file_util.jsonOpen(vep + '.stat_tags.json')
        else:
            cnt_each[chrom] = {}
            cnt_merge[chrom] = {}
            cnt[chrom] = {}

            seq_util.load_refseq_info('b38d')
            chrlen = seq_util.CHROM_LEN['b38d'][chrom]
            flag = True
            spos = 1
            k = 0
            while flag:
                k += 1
                epos = spos + bsize - 1
                if epos > chrlen:
                    epos = chrlen
                
                k_cnt_merge = {}
                k_cnt_each = {}
                k_cnt = {}
                if file_util.is_exist(vep + '.' + str(k) + '.stat_tag.json'):
                    # print(vep + '.' + str(k) + '.stat_tag.json')
                    k_cnt_merge = file_util.jsonOpen(vep + '.' + str(k) + '.stat_mergedtag.json')
                    k_cnt_each = file_util.jsonOpen(vep + '.' + str(k) + '.stat_tag.json')
                    k_cnt = file_util.jsonOpen(vep + '.' + str(k) + '.stat_tags.json')
                else:
                    cmd = "python /home/mk446/bio/mutanno/SRC/scripts/precal_vep/s10_vep_stat.py " + chrom
                    cmd += " " + str(spos)
                    cmd += " " + str(epos)
                    cmd += " " + str(k)
                    print(cmd)

                for f1 in k_cnt_merge.keys():
                    try:
                        cnt_merge[chrom][f1] += k_cnt_merge[f1]
                    except KeyError:
                        cnt_merge[chrom][f1] = k_cnt_merge[f1]

                for f1 in k_cnt_each.keys():
                    try:
                        cnt_each[chrom][f1] += k_cnt_each[f1]
                    except KeyError:
                        cnt_each[chrom][f1] = k_cnt_each[f1]

                for f1 in k_cnt.keys():
                    try:
                        cnt[chrom][f1] += k_cnt[f1]
                    except KeyError:
                        cnt[chrom][f1] = k_cnt[f1]
                spos += bsize
                if epos >= chrlen or spos >= chrlen:
                    break

        for tags in cnt[chrom].keys():
            tagsmap[tags] = 1
        for tag in cnt_each[chrom].keys():
            tagmap[tag] = 1
        for tag in cnt_merge[chrom].keys():
            mtagsmap[tag] = 1

    mtagslist = list(mtagsmap.keys())
    tagslist = list(tagsmap.keys())
    taglist = list(tagmap.keys())

    f = open(statfile, 'w')
    f.write("## VEPmergedtag\n")
    cont = ['VEP_tag']
    cont.append('chr' + '\tchr'.join(seq_util.MAIN_CHROM_LIST))
    cont.append('Total')
    header = '\t'.join(cont) + '\n'
    f.write(header)

    for tag in sorted(mtagslist):
        cont = [tag]
        total = 0
        for chrom in seq_util.MAIN_CHROM_LIST:
            try:
                c1 = cnt_merge[chrom][tag]
            except KeyError:
                c1 = 0
            total += c1
            cont.append(str_util.comma(c1))
        cont.append(str_util.comma(total))
        f.write('\t'.join(cont) + '\n')



    f.write("\n\n\n########################\n")
    f.write("## VEPtags\n")
    cont = ['VEP_tag']
    cont.append('chr' + '\tchr'.join(seq_util.MAIN_CHROM_LIST))
    cont.append('Total')
    header = '\t'.join(cont) + '\n'
    f.write(header)

    for tags in sorted(tagslist):
        cont = [tags]
        total = 0
        for chrom in seq_util.MAIN_CHROM_LIST:
            try:
                c1 = cnt[chrom][tags]
            except KeyError:
                c1 = 0
            total += c1
            cont.append(str_util.comma(c1))
        cont.append(str_util.comma(total))
        f.write('\t'.join(cont) + '\n')

    f.write("\n\n\n########################\n")
    f.write("## VEPtag\n")
    cont = ['VEP_tag']
    cont.append('chr' + '\tchr'.join(seq_util.MAIN_CHROM_LIST))
    cont.append('Total')
    header = '\t'.join(cont) + '\n'
    f.write(header)

    for tag in sorted(taglist):
        cont = [tag]
        total = 0
        for chrom in seq_util.MAIN_CHROM_LIST:
            try:
                c1 = cnt_each[chrom][tag]
            except KeyError:
                c1 = 0
            total += c1
            cont.append(str_util.comma(c1))
        cont.append(str_util.comma(total))
        f.write('\t'.join(cont) + '\n')


    f.close()
    print("Saved",statfile)


def run():
    for chrom in seq_util.MAIN_CHROM_LIST:
        cmd = "python /home/mk446/bio/mutanno/SRC/scripts/precal_vep/s10_vep_stat.py " + chrom
        print(cmd)


def run_more_split():
    seq_util.load_refseq_info('b38d')
    for chrom in seq_util.MAIN_CHROM_LIST:
        chrlen = seq_util.CHROM_LEN['b38d'][chrom]
        flag = True
        spos = 1
        k = 0
        while flag:
            k += 1
            epos = spos + bsize - 1
            if epos > chrlen:
                epos = chrlen
            cmd = "python /home/mk446/bio/mutanno/SRC/scripts/precal_vep/s10_vep_stat.py " + chrom
            cmd += " " + str(spos)
            cmd += " " + str(epos)
            cmd += " " + str(k)
            print(cmd)
            spos += bsize
            if epos >= chrlen or spos >= chrlen:
                break

if __name__ == "__main__":
    import proc_util
    import file_util
    import seq_util
    import str_util
    bsize = 1000000
    vepfile = "/home/mk446/bio/mutanno/DATASOURCE/ANNOT/VEP/hg38/v99_SNV/vep.99.hg38.#CHROM#.tsi.gz"
    statfile = "/home/mk446/bio/mutanno/DATASOURCE/ANNOT/VEP/hg38/v99_SNV/vep.99.hg38.tag_stat"
    if len(sys.argv) == 1:
        # run()
        # save_stat()

        # run_more_split()
        save_stat('split')
    elif len(sys.argv) == 2:
        chrom = sys.argv[1]
        s10_vep_stat(chrom)
    else:
        chrom = sys.argv[1]
        spos = int(sys.argv[2])
        epos = int(sys.argv[3])
        k = sys.argv[4]
        s10_vep_stat_splitrun(chrom, spos, epos, k)
