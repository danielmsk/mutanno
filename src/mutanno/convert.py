#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import tabix
from .util import file_util
from .util import vcf_util
from .util import proc_util


def conv_cgap_filed(f1):
    return f1.replace('+', '').replace('-', '_')
    # return f1.replace('.','_')


class CONVFILE():
    def __init__(self, opt):
        self.opt = opt
        self.datastruct = file_util.load_json(opt['ds'])
        print("Staring... convert.CONVFILE")

    def get_header(self):
        rst = self.get_pars_cont(self.datastruct)
        return '#' + '\t'.join(rst) + "\n"

    def convert(self):
        print("Converting..")
        header = self.get_header()

        if self.opt['chromsplit']:
            fp = {}
            fname = {}
        else:
            fname = self.opt['out'] + '.tsv'
            fpc = open(fname, 'w')
            fpc.write(header)

        i = 0
        prev_chrom = ""

        if 'region' in self.opt.keys() and self.opt['region'] != '':
            tb = tabix.open(self.opt['in'])
            for arr in tb.querys(self.opt['region']):
                i += 1
                line = '\t'.join(arr)
                rst = self.get_pars_cont(self.datastruct, line, [])
                chrom = rst[0]
                if self.opt['chromsplit']:
                    if prev_chrom != chrom:
                        fname[chrom] = self.opt['out'] + '.' + chrom + '.tsv'
                        try:
                            fpc.close()
                        except AttributeError:
                            pass
                        fpc = open(fname[chrom], 'w')
                        fpc.write(header)
                        print("Saving...", fname[chrom])
                fpc.write('\t'.join(rst) + '\n')
                if i % 1000000 == 0:
                    print("...", rst[:2])
                prev_chrom = chrom
        else:
            for line in file_util.gzopen(self.opt['in']):
                if self.opt['in'].endswith('.vcf.gz'):
                    line = line.decode('UTF-8')
                if line[0] != '#':
                    i += 1
                    rst = self.get_pars_cont(self.datastruct, line, [])
                    chrom = rst[0]
                    if self.opt['chromsplit']:
                        if prev_chrom != chrom:
                            fname[chrom] = self.opt['out'] + '.' + chrom + '.tsv'
                            try:
                                fpc.close()
                            except AttributeError:
                                pass
                            fpc = open(fname[chrom], 'w')
                            fpc.write(header)
                            print("Saving...", fname[chrom])
                    fpc.write('\t'.join(rst) + '\n')
                    if i % 1000000 == 0:
                        print("...", rst[:2])
                    prev_chrom = chrom
        fpc.close()

        # TODO: convert this codes using tabix gzip library.
        time.sleep(60)
        if self.opt['chromsplit']:
            for chrom in fp.keys():
                cmd = "tabixgz " + fname[chrom]
                proc_util.run_cmd(cmd, True)
        else:
            cmd = "tabixgz " + fname
            proc_util.run_cmd(cmd, True)

    def conv_map_from_keyvar(self, arrcont):
        m = {}
        for cont in arrcont:
            arr = cont.strip().split('=')
            m[arr[0].strip()] = arr[1].strip()
        return m

    def get_pars_cont(self, ds, cont='', rst=[]):
        if cont != '':
            arrcont = cont.split(ds['delimiter'])
            arrcont[-1] = arrcont[-1].strip()
        for field_ds in ds['fields']:
            if 'fields' not in field_ds.keys():
                if cont != '':
                    if "type" in ds.keys() and ds['type'] == "keyvar":
                        arrcontmap = self.conv_map_from_keyvar(arrcont)
                        if field_ds['name'] in arrcontmap.keys():
                            rst.append(arrcontmap[field_ds['name']].strip())
                        else:
                            rst.append('')
                    else:
                        rst.append(arrcont[field_ds['colidx']].strip())
                else:
                    rst.append(field_ds['name'])
            else:
                if cont != '':
                    rst = self.get_pars_cont(field_ds, arrcont[field_ds['colidx']], rst)
                else:
                    rst = self.get_pars_cont(field_ds, '', rst)
        return rst


class VCF2TSV():
    vcfheader = None

    def __init__(self, opt):
        self.opt = opt
        self.vcfheader = None
        self.datastruct = file_util.load_json(opt['ds'])
        print('Staring... convert.VCF2TSV')

    def pars_vcfline(self, vcfline):
        v1 = vcf_util.pars_vcfline(vcfline, self.vcfheader.collist)
        for infokey in v1['INFO'].keys():
            if infokey in self.vcfheader.mutannofield.keys():
                arrd = []
                for f1 in v1['INFO'][infokey].split(','):
                    d = {}
                    arr = f1.split('|')
                    for i in range(len(self.vcfheader.mutannofield[infokey]['fieldlist'])):
                        d[self.vcfheader.mutannofield[infokey]['fieldlist'][i]] = arr[i].strip()
                    arrd.append(d)
                v1['INFO'][infokey] = arrd
        return v1

    def convert(self):
        out = self.opt['out']
        f = open(out, 'w')

        headercont = ""
        for line in file_util.gzopen(self.opt['in']):
            if self.opt['in'].endswith('.vcf.gz'):
                line = line.decode('UTF-8')

            if line[0] == "#":
                headercont += line
            else:
                if self.vcfheader is None:
                    self.vcfheader = vcf_util.VCFHEADER(headercont)

                v1 = self.pars_vcfline(line)
                # print (v1)
                cont = "### " + v1['CHROM'] + ":" + v1['POS'] + " " + v1['REF'] + ">" + v1['ALT'] + "\n"
                for source in self.vcfheader.mutannofield.keys():
                    if source != "MUTANNO":
                        cont += source + ' '
                        if source in v1['INFO'].keys():
                            cont += "(" + str(len(v1['INFO'][source])) + " fields)\n"
                            for f1 in v1['INFO'][source]:
                                cont += '-----------\n'
                                for fname in f1.keys():
                                    cont += "\t" + fname + " : " + f1[fname] + '\n'
                        else:
                            cont += "(0 fields)\n"

                # print (cont)
                f.write(cont)

                '''
                valuemap = {}
                for field in infofields:
                    if '=' in field:
                        arrf = field.split('=')
                        k1 = arrf[0]
                        if k1 in fieldnames.keys():
                            values = arrf[1].split('|')
                            print (k1, len(values), len(fieldnames[k1])),
                            for k in range(len(fieldnames[k1])):
                                valuemap[k1+"_"+fieldnames[k1][k]] = values[k]
                # print (arr[7])
                valuelist = []
                for k1 in fieldkeylist:
                    for fidx in range(len(fieldnames[k1])):
                        print (k1, fieldnames[k1][fidx])
                        try:
                            valuelist.append( decode_infovalue(valuemap[k1+"_"+fieldnames[k1][fidx]]) )
                        except KeyError:
                            valuelist.append('')
                f.write('\t'.join(valuelist)+'\n')
                '''
        # print (fieldnames)
        # print (colnames)
        f.close()
        print('Saved', out)


class VEP2TAB():
    def __init__(self, opt):
        self.opt = opt
        print('Staring... convert.VEP2TAB')

    def convert(self):
        out = self.opt['out']
        f = open(out, 'w')

        for line in file_util.gzopen(self.opt['in']):
            if self.opt['in'].endswith('.gz'):
                line = line.decode('UTF-8')
            arr = line.split('\t')
            if line[0] == "#":
                if arr[0] == "#Uploaded_variation":
                    cont = ['#CHROM', 'POS', 'REF', 'ALT']
                    arr[0] = arr[0][1:]
                    cont.extend(arr)
                    f.write('\t'.join(cont))
            else:
                arr2 = arr[0].split('_')
                cont = [arr2[0].replace('chr', ''), arr2[1], arr2[2].replace('/', '\t')]
                cont.extend(arr)
                f.write('\t'.join(cont))
        f.close()
        print('Saved ' + out)
        # cmd = "bgzip -c " + out + " > " + out + ".gz;tabix -f -p vcf " + out + ".gz;"
        # proc_util.run_cmd(cmd, True)


if __name__ == "__main__":
    pass
