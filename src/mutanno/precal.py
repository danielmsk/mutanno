from pyfaidx import Fasta
from os import path
from .util import file_util

VEP_OPTION_CMD = "#VEP# -i #INPUTVCF# -o #OUT#"
VEP_OPTION_CMD += " --hgvs"
VEP_OPTION_CMD += " --fasta #FASTA#"
VEP_OPTION_CMD += " --assembly GRCh38"
VEP_OPTION_CMD += " --use_given_ref"
VEP_OPTION_CMD += " --offline"
VEP_OPTION_CMD += " --cache_version #CACHE_VERSION#"
VEP_OPTION_CMD += " --dir_cache #VEPCACHE#"
VEP_OPTION_CMD += " --everything"
VEP_OPTION_CMD += " --force_overwrite"
VEP_OPTION_CMD += " --vcf"
# PLUGIN
# VEP_OPTION_CMD += "--plugin LoF --plugin LoF,human_ancestor_fa:/home/mk446/BiO/Data/vep/human_ancestor.fa.gz "
VEP_OPTION_CMD += " --plugin MaxEntScan,/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99/fordownload"
VEP_OPTION_CMD += " --plugin TSSDistance"
VEP_OPTION_CMD += " --dir_plugins /home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99"
# v99 added
VEP_OPTION_CMD += " --plugin SpliceRegion,Extended"
# VEP_OPTION_CMD += " --plugin LoF,loftee_path:/home/mk446/bio/mutanno/ANNOT3TOOLS/BIN/VEP_plugins-release-99/loftee"


class PreCalculate():
    def __init__(self, opt):
        self.opt = opt

    def getseq(self, chrom, spos, epos):
        pass

    def get_split_region(self, fastaobj, chunksize=1000000):
        chunklist = []
        k2 = 1
        for chrom in list(fastaobj.keys()):
            if len(chrom) < 6:
                no_split = int(1.0 * len(fastaobj[chrom]) / chunksize)
                for k in range(1, no_split + 2):
                    spos = chunksize * (k - 1) + 1
                    epos = chunksize * (k)
                    if epos > len(fastaobj[chrom]):
                        epos = len(fastaobj[chrom])
                    chunklist.append((chrom, spos, epos, k2))
                    k2 += 1
        return chunklist

    def save_input_vcf(self, fastaobj, inputvcf, chrom, spos, epos):
        file_util.check_dir(inputvcf)
        f = open(inputvcf, 'w')
        header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1'
        f.write(header + '\n')
        for k in range(spos, epos + 1):
            ref = fastaobj[chrom][k - 1]
            print(ref)
            if ref != 'N' and ref.strip() != '':
                cont = chrom
                cont += '\t' + str(k)
                cont += '\t' + '.'
                cont += '\t' + ref
                cont += '\t#ALT#'
                cont += '\t' + '.'
                cont += '\t' + '.'
                cont += '\t' + '.'
                cont += '\t' + 'GT:AD:DP'
                cont += '\t' + '0/1:0,30:60'
                cont += '\n'
                for alt in ['A', 'T', 'G', 'C']:
                    if ref != alt:
                        f.write(cont.replace('#ALT#', alt))
        f.close()

    def save_runcmd_vep(self, inputvcf, chrom, spos, epos):
        cmd = ""
        if not file_util.is_exist(inputvcf):
            cmd += "mutanno precal"
            cmd += " -make_input_vcf "
            cmd += " -out " + inputvcf
            for param in ['fasta']:
                cmd += " -" + param + " " + self.opt[param]
            cmd += " -region " + chrom + ':' + str(spos) + '-' + str(epos)
            cmd += ";"
            cmd += "sleep 60;"
        out = self.get_vepout_path(inputvcf)
        vepcmd = VEP_OPTION_CMD
        vepcmd = vepcmd.replace('#VEP#', self.opt['vep'])
        vepcmd = vepcmd.replace('#INPUTVCF#', inputvcf)
        vepcmd = vepcmd.replace('#OUT#', out)
        vepcmd = vepcmd.replace('#FASTA#', path.abspath(self.opt['fasta']))
        vepcmd = vepcmd.replace('#VEPCACHE#', path.abspath(self.opt['vepcache']))
        vepcmd = vepcmd.replace('#CACHE_VERSION#', self.opt['cache_version'])
        cmd += vepcmd + ";"
        cmd += "touch " + inputvcf + ".vep.txt.done;"

        print(cmd)
        # file_util.fileSave(out, cmd, 'w')
        # proc_util.run_cmd('chmod 755 ' + out)

    def get_inputvcf_path(self, chrom, spos, epos):
        out = path.join(self.opt['out'], chrom, str(int(spos / 1000000)),
                        chrom + '_' + str(spos) + '_' + str(epos) + '.vcf')
        out = path.abspath(out)
        return out

    def get_vepout_path(self, inputvcf):
        return inputvcf + '.vep.txt'

    def mk_run_vep_script(self, chunklist):
        for (chrom, spos, epos, k2) in chunklist:
            inputvcf = self.get_inputvcf_path(chrom, spos, epos)
            self.save_runcmd_vep(inputvcf, chrom, spos, epos)
            # if k2 % 1000 == 0:
            #     print("processing..", k2, chrom, spos, inputvcf)
            # break

    # vep_type: vcf, txt
    def check_vep_result(self, vcf, vep_result, output_type='vcf'):
        print('check_vep_result')
        varmap = {}
        for line in file_util.gzopen(vcf):
            if vcf.endswith('.gz'):
                line = line.decode('UTF-8')
            if line[0] != '#':
                arr = line.split('\t')
                # vid = arr[0] + '_' + arr[1] + '_' + arr[3] + '/' + arr[4]
                # (ref, alt, pos) = vcf_util.remove_nonvariant_base(arr[3], arr[4], int(arr[1]), '-')
                pos = int(arr[1])
                ref = arr[3]
                alt = arr[4]
                if output_type == 'txt':
                    if ref[0] == alt[0]:
                        ref = ref[1:]
                        alt = alt[1:]
                        pos += 1
                vid = arr[0] + '_' + str(pos) + '_' + ref + '/' + alt
                # print(vid, ref, alt, pos)
                varmap[vid] = 0

        cont = "no. of variants in VCF: " + str(len(varmap.keys())) + '\n'
        additional_varno = {}
        mismatch_column = {}
        total_lines = 0
        colno = 0
        colidx = {}
        for line in file_util.gzopen(vep_result):
            if vep_result.endswith('.gz'):
                line = line.decode('UTF-8')
            if line[:2] == "##":
                pass
            elif line[0] == "#":
                arr = line[1:].split('\t')
                for k in range(len(arr)):
                    colidx[arr[k]] = k
                colno = len(arr)
                cont += 'no. of column in VEP: ' + str(colno) + '\n'
            elif line[0] != '#':
                arr = line.split('\t')
                if output_type == "txt":
                    vid = arr[colidx['Uploaded_variation']]
                elif output_type == "vcf":
                    vid = arr[0] + '_' + arr[1] + '_' + arr[3] + '/' + arr[4]
                if colno != len(arr):
                    mismatch_column[vid] = 1
                try:
                    varmap[vid] += 1
                except KeyError:
                    additional_varno[vid] = 1
                    print(vid)
                total_lines += 1

        miss_varno = 0
        annot_varno = 0
        for k1 in varmap.keys():
            if varmap[k1] == 0:
                miss_varno += 1
            else:
                annot_varno += 1
        cont += "total lines in VEP: " + str(total_lines) + '\n'
        cont += "no. of annotated variants in VEP: " + str(annot_varno) + '\n'
        cont += "no. of missing variants in VEP: " + str(miss_varno) + '\n'
        cont += "no. of variants with mismatching column in VEP: " + str(len(mismatch_column.keys())) + '\n'
        cont += "no. of additional variants in VEP: " + str(len(additional_varno.keys())) + '\n'
        print(cont)
        if miss_varno > 0 or len(mismatch_column.keys()) > 0:
            out = vep_result + '.error'
        else:
            out = vep_result + '.checked'
        file_util.fileSave(out, cont, 'w')
        print("Saved", out)

    def mk_shellscript_merge_vep(self, fasta, out, chunklist):
        for (chrom, spos, epos, k2) in chunklist:
            inputvcf = self.get_inputvcf_path(chrom, spos, epos)
            print(inputvcf)
        pass

    def run(self):
        chunksize = 100000
        if self.opt['fasta'] != "":
            fasta = Fasta(self.opt['fasta'], as_raw=True, sequence_always_upper=True)

        if self.opt['run_vep']:
            chunklist = self.get_split_region(fasta, chunksize)
            print('Total chunk files:' + str(len(chunklist)))
            self.mk_run_vep_script(chunklist)
        if self.opt['make_input_vcf']:
            chrom = self.opt['region'].split(':')[0]
            poslist = self.opt['region'].split(':')[1].split('-')
            spos = int(poslist[0])
            epos = int(poslist[1])
            self.save_input_vcf(fasta, self.opt['out'], chrom, spos, epos)
        if self.opt['merge_vep']:
            chunklist = self.get_split_region(fasta, chunksize)
            self.mk_shellscript_merge_vep(fasta, self.opt['out'], chunklist)
        if self.opt['check_vep_result']:
            self.check_vep_result(self.opt['vcf'], self.opt['vep_result'], 'vcf')
