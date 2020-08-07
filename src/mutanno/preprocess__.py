#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .util import file_util
from .util import proc_util


class MakeDbnsfpTranscript():
    def __init__(self, infile, outfile, ds, opt):
        self.infile = infile
        self.outfile = outfile
        if self.outfile.endswith('.gz'):
            self.outfile2 = outfile[:-3]
        else:
            self.outfile2 = outfile
        self.datastruct = self.get_dbNSFP_Transcipt_struct(ds)

    def get_dbNSFP_Transcipt_struct(self, ds):
        dstruct = file_util.load_json(ds)
        for source in dstruct['source']:
            if source['name'] == "dbNSFP_transcript":
                rst = source
        return rst

    def run(self):
        # print (self.datastruct)

        f = open(self.outfile2, 'w')
        headermap = {}
        i = 0
        for line in file_util.gzopen(self.infile):
            i += 1
            if self.infile.endswith('.gz'):
                line = line.decode('UTF-8')

            arr = line.split('\t')
            arr[-1] = arr[-1].strip()
            if line[0] == '#':
                for k in range(len(arr)):
                    headermap[arr[k].strip()] = k

            enstarr = arr[headermap['Ensembl_transcriptid']].strip().split(';')

            if len(enstarr) > 1:
                # print ("=>", arr)
                for eidx in range(len(enstarr)):
                    cont = []
                    for field in self.datastruct['fields']:
                        arrv1 = arr[headermap[field['name']]].strip().split(';')

                        if len(arrv1) == 1:
                            v1 = arrv1[0]
                        else:
                            # if field['name'] == 'MutationTaster_score':
                            #     print (len(enstarr), field['name'], arrv1, arr[headermap['MutationTaster_pred']],
                            #             arr[headermap['MutationTaster_model']])
                            v1 = arrv1[eidx].strip()
                        cont.append(v1)
                    f.write('\t'.join(cont) + '\n')
                # print()
                # break
            else:
                cont = []
                for field in self.datastruct['fields']:
                    v1 = arr[headermap[field['name']]].strip()
                    cont.append(v1)
                f.write('\t'.join(cont) + '\n')

            if i % 10000 == 0:
                print(i, arr[0], arr[1])
                pass
                # break

        f.close()
        if self.outfile.endswith('.gz'):
            proc_util.run_cmd('tabixgz ' + self.outfile2)


if __name__ == "__main__":
    pass
