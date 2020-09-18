from ..util import file_util, vcf_util
from ..model.datasource import DataSourceList


class PreprocVEP():
    def __init__(self, data):
        self.selected_fields = []
        self.infile = data['infile']
        self.dslist = None
        if data['ds'] != '':
            self.dslist = DataSourceList(data['ds'])
        self.out = data['out']
        self.set_outfile_extension()

    def set_outfile_extension(self):
        out2 = self.out.replace('.mti.gz', '.mti')
        if out2[-4:] != ".mti":
            out2 += ".mti"
        self.out = out2

    def get_mti_header(self):
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        info_arr = []
        for s1 in self.dslist.available_source_list:
            info_header = []
            for f1 in s1.available_field_list:
                info_header.append(f1.name2)
                if s1.name == "VEP":
                    self.selected_fields.append(f1.name2)
            info_arr.append(s1.name + '=' + '|'.join(info_header))
        header.append(';'.join(info_arr))
        return '\t'.join(header)

    def get_mti_header_nonselection(self, vep_field_list):
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        header.append('VEP=' + '|'.join(vep_field_list))
        return '\t'.join(header)


    def select_fields(self, info_field, vep_field_list, alt, flag_select=False):
        sel_idx_list = []
        for sel_field in self.selected_fields:
            sel_idx_list.append(vep_field_list.index(sel_field))

        rst = "VEP="
        for field in info_field.split(';'):
            if field[:len("CSQ=")] == "CSQ=":
                rst_section = []
                for section in field.replace('CSQ=', '').split(','):
                    arr = section.split('|')
                    allele = arr[0]
                    if (flag_select and alt[1:] == allele) or not flag_select:
                        rst_field = []
                        for idx in sel_idx_list:
                            rst_field.append(self.encode_vep_field(arr[idx], vep_field_list[idx]))
                        rst_section.append('|'.join(rst_field))
                rst += ','.join(rst_section)
        return rst

    def encode_vep_field(self, value, fieldname):
        if fieldname == "Consequence":
            value = value.replace('&', '~')
        else:
            value = vcf_util.encode_value(value)
        return value

    def convert_vep2mti(self):
        fp = open(self.out, 'w')

        for line in file_util.gzopen(self.infile):
            line = file_util.decodeb(line)
            if line[0] == '#':
                if '##INFO=<ID=CSQ' in line:
                    vep_field_list = line.split('Format:')[-1].replace('">', '').strip().split('|')
                    if self.dslist == None:
                        fp.write(self.get_mti_header_nonselection(vep_field_list) + '\n')
                        self.selected_fields = vep_field_list
                    else:
                        fp.write(self.get_mti_header() + '\n')
            else:
                arr = line.split('\t')
                arralt = arr[4].split(',')

                for alt in arralt:
                    rst = []
                    rst.append(arr[0].replace('chr', ''))
                    rst.append(arr[1])
                    rst.append('.')
                    rst.append(arr[3])
                    rst.append(alt)
                    rst.append('')
                    rst.append('')
                    rst.append(self.select_fields(arr[7], vep_field_list, alt, len(arralt) > 1))
                    line = '\t'.join(rst) + '\n'
                    fp.write(line)
        fp.close()
        print("Saved", self.out)


    def run(self):
        print('preprocVEP run..')
        self.convert_vep2mti()

        
        
