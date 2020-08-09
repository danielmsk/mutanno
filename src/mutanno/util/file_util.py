#!/usr/bin/env python
# -*- coding: utf-8 -*-
# file_util.py
import os
import sys
import tabix
import json
from . import proc_util
import pprint
import mutanno


def getDataPath(datafile):
    return (os.path.join(mutanno.__path__[0], 'data', datafile))


def print(str):
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(str)


def save_gzip(out):
    cmd = "bgzip -c " + out + " > " + out + ".gz"
    proc_util.run_cmd(cmd)


def strip_gzext(fname):
    rst = fname
    if rst.endswith('.gz'):
        rst = rst[:-3]
    return rst


def line2arr(line, delimiter='\t'):
    arr = line.split('\t')
    arr[-1] = arr[-1].strip()
    return arr


def zero_format(i, size):
    # 00001000 format
    s1 = "{:0>" + str(size) + "d}"
    return s1.format(i)


def load_json(jsonfile):
    ds = ""
    with open(jsonfile) as jfp:
        ds = json.load(jfp)
    return ds


def get_tblist(vcfgzlist):
    tblist = []
    for vcfgz in vcfgzlist:
        tblist.append(tabix.open(vcfgz))
    return tblist


def getRowTabix(tb, chrom, pos, ref="", alt=""):
    if not isinstance(pos, int):
        pos = int(pos)
    recs = tb.query(chrom, pos - 1, pos)
    rec = []
    for r1 in recs:
        if int(r1[1]) == pos:
            if ref != "" and alt != "":
                if r1[3] == ref and r1[4].replace(",<NON_REF>", "") == alt:
                    rec = r1
                    break
            elif ref != "":
                if r1[3] == ref:
                    rec = r1
                    break
            else:
                rec = r1
                break

    return rec


def check_dir(fname):
    if fname[0] != '/' and fname[0] != '.':
        fname = './' + fname
    arr = fname.split("/")
    fpath = arr[0]
    for d in arr[1:-1]:
        fpath += "/" + d
        if not is_exist(fpath):
            mkDir(fpath)
            # print("make directory : " + fpath)


def fileOpenOffset(filepath, start_pos, end_pos):
    f = open(filepath)
    f.seek(start_pos)
    cont = f.read(end_pos - start_pos)
    return cont


def fileOpen2(path, opt):
    try:
        thisFile = open(path, opt)
    except IOError:
        print(sys.stderr, "File could not be opened")
        sys.exit(1)
    else:
        return thisFile


def fileOpen(path):
    cont = ""
    if path.endswith('.gz'):
        for line in gzopen(path):
            cont += decodeb(line)
    else:
        f = open(path, "r")
        cont = f.read()
    return cont


def fileSave(path, cont, opt, gzip_flag="n"):
    if gzip_flag == "gz":
        import gzip
        f = gzip.open(path, opt)
        f.write(cont)
        f.close()
    else:
        f = open(path, opt)
        f.write(cont)
        f.close


def gzsave(path, cont, opt):
    return fileSave(path, cont, opt, "gz")


def listdir(dirPath, ext=""):
    flist = []
    for fname in os.listdir(dirPath):
        if (len(ext) > 0 and fname.endswith(ext)) or len(ext) == 0:
            flist.append(fname)
    flist.sort()
    return flist


def walk(dirPath, ext=""):
    flist = []
    for root, dirs, files in os.walk(dirPath):
        for fname in files:
            if (len(ext) > 0 and fname.endswith(ext)) or len(ext) == 0:
                fullpath = os.path.join(root, fname)
                flist.append(fullpath)
    flist.sort()
    return flist


def walkdir(dirPath, cls="simple"):
    dlist = []
    for root, dirs, files in os.walk(dirPath):
        for dname in dirs:
            if dname != "":
                if cls == "simple":
                    dlist.append(dname)
                else:
                    dlist.append(os.path.join(root, dname))
    dlist.sort()
    return dlist


def read_table_col(filename, delemeter="\t"):
    tab = {}
    i = 0
    header = []
    repeatmap = {}
    for line in gzopen(filename):
        if len(line.strip()) > 0:
            arr = line.split(delemeter)
            arr[-1] = arr[-1].strip()

            if i == 0:
                for j in range(0, len(arr)):
                    h1 = arr[j].strip()
                    if h1 in header:
                        header.append(h1 + "_" + str(repeatmap[h1]))
                        repeatmap[h1] += 1
                    else:
                        header.append(h1)
                        repeatmap[h1] = 1
            else:
                for j in range(0, len(header)):
                    try:
                        tab[header[j]].append(arr[j].strip())
                    except KeyError:
                        tab[header[j]] = []
                        tab[header[j]].append(arr[j].strip())
        i += 1
    return tab


def read_table(filename, delemeter="\t"):
    tab = []
    i = 0
    header = []
    repeatmap = {}
    for line in gzopen(filename):
        if filename.endswith('.gz'):
            line = line.decode('UTF-8')
        if len(line.strip()) > 0:
            arr = line.split(delemeter)
            arr[-1] = arr[-1].strip()

            if i == 0:
                for j in range(0, len(arr)):
                    h1 = arr[j].strip()
                    if h1 in header:
                        header.append(h1 + "_" + str(repeatmap[h1]))
                        repeatmap[h1] += 1
                    else:
                        header.append(h1)
                        repeatmap[h1] = 1
            else:
                t1 = {}
                for j in range(0, len(header)):
                    t1[header[j]] = arr[j].strip()
                tab.append(t1)
        i += 1
    return tab


def read_table_key(filename, key_col_idx=0, delemeter="\t"):
    tab = {}
    i = 0
    header = []
    repeatmap = {}
    for line in open(filename):
        if len(line.strip()) > 0:
            arr = line.split(delemeter)
            arr[-1] = arr[-1].strip()

            if i == 0:
                for j in range(0, len(arr)):
                    h1 = arr[j].strip()
                    if h1 in header:
                        header.append(h1 + "_" + str(repeatmap[h1]))
                        repeatmap[h1] += 1
                    else:
                        header.append(h1)
                        repeatmap[h1] = 1
            else:
                t1 = {}
                for j in range(0, len(header)):
                    t1[header[j]] = arr[j].strip()
                tab[arr[key_col_idx]] = t1
        i += 1
    return tab


def comp_gz(target_file, gzip_file=""):
    import gzip
    if gzip_file == "":
        gzip_file = target_file + ".gz"
    f_in = open(target_file, 'rb')
    f_out = gzip.open(gzip_file, 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()


def mkDir(dpath):
    if not os.path.isdir(dpath):
        os.mkdir(dpath)


def is_exist(fpath):
    return os.path.exists(fpath)


def getFileInfo(fpath):
    st = {}
    (st['mode'], st['ino'], st['dev'], st['st_nlink'], st['uid'], st['gid'],
     st['size'], st['atime'], st['mtime'], st['ctime']) = os.stat(fpath)
    return st


def getFileSize(fpath):
    st = {}
    (st['mode'], st['ino'], st['dev'], st['st_nlink'], st['uid'], st['gid'],
     st['size'], st['atime'], st['mtime'], st['ctime']) = os.stat(fpath)
    return st['size']


def copy(src, dst):
    import shutil
    shutil.copy(src, dst)


def rm(fname):
    if is_exist(fname):
        os.remove(fname)


def fopen(fname):
    if fname.endswith(".gz"):
        import gzip
        f1 = gzip.GzipFile(fname, "r")
    else:
        f1 = open(fname)
    return f1


def gzopen(fname):
    if fname.endswith(".gz"):
        import gzip
        f1 = gzip.GzipFile(fname, "r")
    else:
        f1 = open(fname)
    return f1


def extract_line(fname, findme, column=-9, delimiter="space", outfile="", header=False):
    import str_util
    if fname[-3:] == ".gz":
        f1 = gzopen(fname)
    else:
        f1 = open(fname)

    if outfile != "":
        fileSave(outfile, "", "w")

    i = 0
    for line in f1:
        if header and i == 0:
            cont = line.strip()
        else:
            cont = ""
            if column < 0:
                if findme in line:
                    cont = line.strip()
            else:
                if delimiter == "space":
                    line = str_util.strip_space(line)
                    arr = line.split(" ")
                else:
                    arr = line.strip().split(delimiter)
                if findme in arr[column]:
                    cont = line.strip()

        if cont != "" and outfile != "":
            fileSave(outfile, cont + "\n", "a")
        elif cont != "":
            print(cont)

        i += 1


# Function Name : load_seies_file
# load data file with series format in GEO database
def load_series_file(fname, coat="n"):
    m = {}
    cont = ""
    for line in gzopen(fname):
        if line[0] == "!":
            arr = line[1:].strip().split('"\t"')
            arr_k = arr[0].split("\t")
            k = arr_k[0].strip()
            if len(arr_k) > 1:
                v = [arr_k[1].strip()]
            else:
                v = []
            if len(arr) > 1:
                v.extend(arr[1:])
            for i in range(len(v)):
                v[i] = v[i].replace('"', '')
            if len(v) > 0:
                try:
                    arr2 = m[k]
                    arr2.append(v)
                    '''
                    if len(arr) == 2:
                        arr2.append(v)
                    else:
                        arr2.append(v)
                    '''
                except KeyError:
                    m[k] = []
                    m[k].append(v)
                    '''
                    if len(arr) == 2:
                        m[k].append(v)
                    else:
                        m[k].append(v)
                    '''
            else:
                if len(cont) > 0:
                    k = k.replace("_end", "")
                    try:
                        arr2 = m[k]
                        arr2.append(cont)
                        m[k] = arr2
                    except KeyError:
                        m[k] = []
                        m[k].append(cont)
                    cont = ""
        else:
            cont += line
    return m


def log(msg, logfile="", silence=False):
    import time_util

    if "#DATE#" in logfile:
        logfile = logfile.replace("#DATE#", time_util.getToday())

    cont = "[" + time_util.getNow() + "] " + msg + "\n"

    if not silence:
        print(cont,)

    if logfile != "":
        fileSave(logfile, cont, "a")


def get_header(fname):
    arr = []
    for line in gzopen(fname):
        if fname.endswith(".gz"):
            line = line.decode('UTF-8')
        arr = line.split('\t')
        arr[-1] = arr[-1].strip()
        break
    return arr


def get_header_idx(fname):
    h = {}
    arr = get_header(fname)
    for k in range(len(arr)):
        h[arr[k]] = k
    return h


def decodeb(bstr):
    if isinstance(bstr, bytes):
        rst = bstr.decode('UTF-8')
    else:
        rst = bstr
    return rst


def rmext(fname, ext):
    return fname[:(-1 * len(ext))]


def trim_split(v1, delimeter):
    rst = []
    for a1 in v1.split(delimeter):
        rst.append(a1.strip())
    return rst


def tmp_load_clinvar_idmap():
    clinvar_idmapfile = getDataPath('CLINVAR_hg38_20200329_submission.sorted.tsi.gz_idmap.txt.gz')
    clinvar_idmap = {}
    for line in gzopen(clinvar_idmapfile):
        line = decodeb(line)
        arr = line.strip().split('\t')
        clinvar_idmap[arr[0]] = arr[1]
    return clinvar_idmap
