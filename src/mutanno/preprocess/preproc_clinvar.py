import os
import xmltodict
from ..util import file_util

VCFCOL = ["CHROM", "POS", "ALLELEID", "REF", "ALT", "VARIATIONID", "CLNDN", "CLNDNINCL", "CLNDISDB",
          "CLNDISDBINCL", "CLNHGVS", "CLNREVSTAT", "CLNSIG", "CLNSIGCONF", "CLNSIGINCL", "CLNVC",
          "CLNVCSO", "GENEINFO", "MC", "ORIGIN"]
# XMLCOL = ["ClinVarAccession","Interpretation","DateLastEvaluated","ReviewStatus","Method","Condition","AlleleOrigin",
# "Submitter","SubmitterID","Citation","Comment"]
XMLCOL = ["ClinVarAccession", "Interpretation", "DateLastEvaluated", "ReviewStatus",
          "AssertionMethod", "AssertionMethodCitation", "Method", "Condition", "ConditionXRef",
          "AlleleOrigin", "Submitter", "SubmitterID", "Comment", "Citation"]
XMLCOL_LIST_TYPE = ["Comment", "Citation", "Condition", "ConditionXRef"]


def conv_list(list1):
    if type(list1) != list:
        if list1 == '':
            list1 = []
        else:
            list1 = [list1]
    return list1


def get_dict_data(dict1, keylist):
    rst = ""
    rst1 = dict1
    for i, key1 in enumerate(keylist):
        if type(rst1) == list:
            rst1 = rst1[0]
        try:
            rst1 = rst1[key1]
            rst = rst1
        except KeyError:
            rst = ""
            break
    return rst


def cnt_field(cnt, f1, v1):
    try:
        cnt[f1]
    except KeyError:
        cnt[f1] = {}
    try:
        cnt[f1][v1] += 1
    except KeyError:
        cnt[f1][v1] = 1
    return cnt


CONV_SIGN_LIST = []
CONV_SIGN_LIST.append((';', "%3B"))
CONV_SIGN_LIST.append(('=', "%3D"))
CONV_SIGN_LIST.append(("|", "%7C"))
CONV_SIGN_LIST.append((",", "%2C"))
CONV_SIGN_LIST.append(('"', "%22"))
CONV_SIGN_LIST.append(('~', "%7E"))
CONV_SIGN_LIST.append((' ', '%20'))


def encode_value(v1, except_list=[]):
    for s1 in CONV_SIGN_LIST:
        if s1[0] not in except_list:
            v1 = v1.replace(s1[0], s1[1])
    # v1 = v1.replace(';', "%3B").replace('=', "%3D").replace("|", "%7C").replace(",", "%2C")
    # v1 = v1.replace('"', "%22").replace('~', "%7E").replace(' ', '%20')
    return v1


class PreprocClinVar():
    def __init__(self, data):
        self.rawfiles = data['rawfiles']
        # self.file_version = data['file_version']
        # self.refversion = data['refversion']
        # self.source_name = data['source_name']
        self.dirpath = data['dirpath']
        self.outfile_title = data['outfile_title']
        self.clinvarxml, self.clinvarvcf = self.get_rawfiles(self.rawfiles)
        self.out_variant, self.out_submission = self.get_out_filename()
        self.ref_trait = {}

    def get_out_filename(self):
        out_variant = os.path.join(self.dirpath, self.outfile_title + '_variant.tsi')
        out_submission = os.path.join(self.dirpath, self.outfile_title + '_submission.tsi')
        return out_variant, out_submission

    def get_rawfiles(self, rawfiles):
        if '.xml' in rawfiles[0]:
            xml = rawfiles[0]
            vcf = rawfiles[1]
        else:
            xml = rawfiles[1]
            vcf = rawfiles[0]
        return xml, vcf

    def save_clinvarvcf(self):
        global VCFCOL

        i = 0
        cvar = {}
        cnt = {}
        f = open(self.out_variant, 'w')
        f.write('#'+'\t'.join(VCFCOL[:5]) + '\tCLINVAR=' + '|'.join(VCFCOL[5:]) + '\n')
        for line in file_util.gzopen(self.clinvarvcf):
            line = file_util.decodeb(line)
            if line[0] != '#':
                i += 1
                arr = line.split('\t')
                m = {}
                for f1 in VCFCOL:
                    m[f1] = ""
                m['CHROM'] = arr[0].strip()
                if m['CHROM'] == "MT":
                    m['CHROM'] = "M"
                m['POS'] = arr[1].strip()
                m['VARIATIONID'] = arr[2].strip()
                m['REF'] = arr[3].strip()
                m['ALT'] = arr[4].strip()
                for f1 in arr[7].strip().split(';'):
                    arr2 = f1.split('=')
                    m[arr2[0]] = arr2[1].strip()

                cont = []
                for f1 in VCFCOL[:5]:
                    cont.append(m[f1])

                info = []
                for idx in range(5, len(VCFCOL)):
                    f1 = VCFCOL[idx]
                    if f1 == "MC":
                        arr2 = m[f1].split(',')
                        arr3 = []
                        for a2 in arr2:
                            arr3.append(a2.split('|')[0])
                        m[f1] = '~'.join(arr3)
                    elif f1 == "CLNDISDB":
                        m[f1] = m[f1].replace('.|', '').replace('.|', '').replace(
                            '.|', '').replace('.|', '').replace('.|', '').replace('.|', '')
                        m[f1] = m[f1].replace('.|', '').replace('.|', '').replace(
                            '.|', '').replace('.|', '').replace('.|', '').replace('.|', '')
                        m[f1] = m[f1].replace('.', '')
                        m[f1] = m[f1].replace('|', '~').replace(',', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "CLNDN":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = m[f1].replace('|', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "CLNDNINCL":
                        m[f1] = m[f1].replace('_', '%20')
                        # m[f1] = m[f1].replace('|', '~')
                        m[f1] = encode_value(m[f1], [])
                    elif f1 == "CLNREVSTAT":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = m[f1].replace(',', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "CLNSIG":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = encode_value(m[f1], [])
                        # print(m[f1])
                    elif f1 == "CLNSIGINCL":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = m[f1].replace('|', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "CLNSIGCONF":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = m[f1].replace(',', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    # elif f1 == "CLNVI":
                    #     m[f1] = m[f1].replace('|', '~')
                    #     m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "CLNVC":
                        m[f1] = m[f1].replace('_', '%20')
                        m[f1] = encode_value(m[f1], [])
                    elif f1 == "CLNDISDBINCL":
                        m[f1] = m[f1].replace('.|', '').replace('.|', '').replace(
                            '.|', '').replace('.|', '').replace('.|', '').replace('.|', '')
                        m[f1] = m[f1].replace('.|', '').replace('.|', '').replace(
                            '.|', '').replace('.|', '').replace('.|', '').replace('.|', '')
                        m[f1] = m[f1].replace('.', '')
                        m[f1] = m[f1].replace('|', '~').replace(',', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    elif f1 == "GENEINFO":
                        m[f1] = m[f1].replace('|', '~')
                        m[f1] = encode_value(m[f1], ['~'])
                    else:
                        m[f1] = encode_value(m[f1])
                    info.append(m[f1])

                    cnt = cnt_field(cnt, f1, m[f1])

                cont.append('CLINVAR='+'|'.join(info))
                f.write('\t'.join(cont) + '\n')
                cvar[m['VARIATIONID']] = [m['CHROM'], m['POS'], m['VARIATIONID'], m['REF'], m['ALT']]

                if i % 100000 == 0:
                    print(i, m)

        f.close()
        print('Saved', self.out_variant, len(cvar.keys()))

        for k1 in VCFCOL[5:]:
            statfile = self.out_variant + "_stat/stat_" + k1.replace(' ', '_') + '.cnt'
            file_util.check_dir(statfile)
            cont = ''
            ks = cnt[k1].keys()
            for k2 in sorted(ks):
                cont += k2 + '\t' + str(cnt[k1][k2]) + '\n'
            file_util.fileSave(statfile, cont, 'w')

        return cvar

    def update_ref_trait(self, xml):
        ref_trait_list = ref_trait_name_list = get_dict_data(xml, ['ClinVarSet', 'ReferenceClinVarAssertion',
                                                                   'TraitSet', 'Trait'])
        for ref_trait_xml in conv_list(ref_trait_list):
            ref_trait_name_list = get_dict_data(ref_trait_xml, ['Name'])
            for ref_trait_name in conv_list(ref_trait_name_list):
                etype = get_dict_data(ref_trait_name, ['ElementValue', '@Type'])
                if etype == "Preferred":
                    ref_trait_name_text = get_dict_data(ref_trait_name, ['ElementValue', '#text'])
                    for ref_trait_xref in conv_list(get_dict_data(ref_trait_xml, ['XRef'])):
                        ref_xref_db = get_dict_data(ref_trait_xref, ['@DB'])
                        ref_xref_id = get_dict_data(ref_trait_xref, ['@ID'])
                        xid = ref_xref_id if ':' in ref_xref_id else ref_xref_db + ":" + ref_xref_id
                        if ref_trait_name_text != '' and ref_xref_db != '':
                            self.ref_trait[ref_trait_name_text.strip()] = xid
                            self.ref_trait[xid] = ref_trait_name_text.strip()

        # print(">>>>>self.ref_trait:", len(self.ref_trait.keys()))

    def get_submission_from_xml(self, setcont, cvar):
        xml = xmltodict.parse(setcont)
        cset = xml['ClinVarSet']
        asetlist = get_dict_data(xml, ['ClinVarSet', 'ClinVarAssertion'])

        if type(asetlist) != list:
            asetlist = [asetlist]

        varid = get_dict_data(cset, ['ReferenceClinVarAssertion', 'MeasureSet', '@ID'])
        try:
            cvar[varid]
        except KeyError:
            return []

        submissionlist = []
        for aset in asetlist:
            m = {}
            for f1 in XMLCOL:
                m[f1] = ''

            m['VARIATIONID'] = varid
            # col1: Interpretation
            m['Interpretation'] = get_dict_data(aset, ['ClinicalSignificance', 'Description'])
            m['DateLastEvaluated'] = get_dict_data(aset, ['ClinicalSignificance', '@DateLastEvaluated'])

            # col2: Review status
            m['ReviewStatus'] = get_dict_data(aset, ['ClinicalSignificance', 'ReviewStatus'])
            m['Method'] = get_dict_data(xml, ['ClinVarSet', 'ReferenceClinVarAssertion',
                                              'ObservedIn', 'Method', 'MethodType'])
            # AssertionMethod
            m['AssertionMethod'] = ""
            m['AssertionMethodCitation'] = ""
            attributeset_list = get_dict_data(aset, ['AttributeSet'])
            for attributeset in conv_list(attributeset_list):
                att_type = get_dict_data(attributeset, ['Attribute', '@Type'])
                if att_type == "AssertionMethod":
                    m['AssertionMethod'] = get_dict_data(attributeset, ['Attribute', '#text'])
                    assert_method_url = get_dict_data(attributeset, ['Citation', 'URL'])
                    assert_method_source = get_dict_data(attributeset, ['Citation', 'ID', '@Source'])
                    assert_method_txt = get_dict_data(attributeset, ['Citation', 'ID', '#text'])
                    if assert_method_url != "":
                        m['AssertionMethodCitation'] = "URL:" + assert_method_url
                    if assert_method_source != "":
                        m['AssertionMethodCitation'] = assert_method_source + ":" + assert_method_txt
                m['AssertionMethod'] = m['AssertionMethod'].replace('_', '%20')

            # col3: Condition
            self.update_ref_trait(xml)

            condition_list = []
            conditionxref_list = []
            for trait_xml in conv_list(get_dict_data(aset, ['TraitSet', 'Trait'])):
                condition = ""
                xid = ""
                for trait_name in conv_list(get_dict_data(trait_xml, ['Name'])):
                    etype = get_dict_data(trait_name, ['ElementValue', '@Type'])
                    if etype == "Preferred":
                        condition = get_dict_data(trait_name, ['ElementValue', '#text'])
                        break

                for xref_xml in conv_list(get_dict_data(trait_xml, ['XRef'])):
                    xref_db = get_dict_data(xref_xml, ['@DB'])
                    xref_id = get_dict_data(xref_xml, ['@ID'])
                    if xref_id != "":
                        xid = xref_id if ':' in xref_id else xref_db + ":" + xref_id
                        break
                if condition == "" and xid != "":
                    try:
                        condition = self.ref_trait[xid]
                    except KeyError:
                        pass
                elif condition != "" and xid == "":
                    try:
                        xid = self.ref_trait[condition]
                    except KeyError:
                        pass

                if condition != "" and xid != "":
                    condition_list.append(condition)
                    conditionxref_list.append(xid)

            m['Condition'] = '~'.join(condition_list)
            m['ConditionXRef'] = '~'.join(conditionxref_list)

            m['AlleleOrigin'] = get_dict_data(
                xml, ['ClinVarSet', 'ReferenceClinVarAssertion', 'ObservedIn', 'Sample', 'Origin'])

            # col4: Submitter
            m['Submitter'] = get_dict_data(aset, ['ClinVarSubmissionID', '@submitter'])
            # https://www.ncbi.nlm.nih.gov/clinvar/submitters/<ID>/
            m['SubmitterID'] = get_dict_data(aset, ['ClinVarAccession', '@OrgID'])
            m['ClinVarAccession'] = get_dict_data(aset, ['ClinVarAccession', '@Acc'])

            # col5: SupportingInfo

            comment = []
            for clin_sig_comment_xml in conv_list(get_dict_data(aset, ['ClinicalSignificance', 'Comment'])):
                try:
                    comment.append(get_dict_data(clin_sig_comment_xml, ['#text']))
                    # comment.append(clin_sig_comment_xml)
                except TypeError:
                    comment.append(clin_sig_comment_xml)
                    # print(clin_sig_comment_xml)

            m['Comment'] = '~'.join(comment)

            supinfo = []
            for clin_sig_citation_xml in conv_list(get_dict_data(aset, ['ClinicalSignificance', 'Citation'])):
                sInfo_DB = get_dict_data(clin_sig_citation_xml, ['ID', '@Source'])  # ---------------
                sInfo_txt = get_dict_data(clin_sig_citation_xml, ['ID', '#text'])  # ---------------
                sInfo_url = get_dict_data(clin_sig_citation_xml, ['URL'])  # ---------------
                if sInfo_DB != '':
                    supinfo.append(sInfo_DB + ':' + sInfo_txt)
                if sInfo_url != '':
                    supinfo.append('URL:' + sInfo_url)
            for clin_citation_xml in conv_list(get_dict_data(aset, ['Citation'])):
                sInfo_DB = get_dict_data(clin_citation_xml, ['ID', '@Source'])  # ---------------
                sInfo_txt = get_dict_data(clin_citation_xml, ['ID', '#text'])  # ---------------
                sInfo_url = get_dict_data(clin_citation_xml, ['URL'])  # ---------------
                if sInfo_DB != '':
                    supinfo.append(sInfo_DB + ':' + sInfo_txt)
                if sInfo_url != '':
                    supinfo.append('URL:' + sInfo_url)
            m['Citation'] = '~'.join(supinfo)
            submissionlist.append(m)
        return submissionlist

    def save_clinvarxml(self, cvar):
        f = open(self.out_submission, 'w')
        cont = ['#CHROM', 'POS', 'VARIATIONID', 'REF', 'ALT', "CLINVAR_SUBMISSION="+'|'.join(XMLCOL)]
        f.write('\t'.join(cont) + '\n')
        flag = False
        i = 0
        cnt = {}
        for line in file_util.gzopen(self.clinvarxml):
            line = file_util.decodeb(line)

            if '<ClinVarSet ' in line:
                flag = True
                setcont = ""

            if flag:
                setcont += line

            if '</ClinVarSet>' in line:
                i += 1

                submissionlist = self.get_submission_from_xml(setcont, cvar)
                for m in submissionlist:
                    try:
                        vpos = cvar[m['VARIATIONID']]
                        xmlfield = []
                        for f1 in XMLCOL:
                            if f1 in XMLCOL_LIST_TYPE:
                                m[f1] = encode_value(m[f1], ['~'])
                            else:
                                m[f1] = encode_value(m[f1])
                            xmlfield.append(m[f1])
                            cnt = cnt_field(cnt, f1, m[f1])

                        cont = '\t'.join(vpos) + '\t' + "CLINVAR_SUBMISSION=" + '|'.join(xmlfield)
                        f.write(cont + '\n')
                        # print(cont)
                    except KeyError:
                        pass

                    flag = False
                    if i % 10000 == 0:
                        print(i, m)
                        # break

        f.close()

        for k1 in XMLCOL:
            statfile = self.out_submission + "_stat/stat_" + k1.replace(' ', '_') + '.cnt'
            file_util.check_dir(statfile)
            cont = ''
            ks = cnt[k1].keys()
            for k2 in sorted(ks):
                # print(k1, k2, cnt[k1][k2])
                cont += k2 + '\t' + str(cnt[k1][k2]) + '\n'
            file_util.fileSave(statfile, cont, 'w')

        print("Saved", self.out_submission)

        # cont = ""
        # for line in open(self.out_submission):
        #     arr = line.split('\t')
        #     cont += arr[-1].strip().replace('CLINVAR_SUBMISSION=', '').replace('|', '\t') + "\n"
        # file_util.fileSave('test.txt', cont, 'w')
        # print(cont)

    def run(self):
        cvar = self.save_clinvarvcf()
        self.save_clinvarxml(cvar)
