from .util import vep_util, file_util


def void():
    pass


def conv_consequence_delimiter(consequence):
    rst = consequence
    if isinstance(consequence, str):
        rst = consequence.replace('&', '~')
    elif isinstance(consequence, list):
        rst = []
        for c1 in consequence:
            rst.append(c1.replace('&', '~'))
    return rst


def get_score_from_vep_pathogenicity(predscore):
    rst = predscore
    if isinstance(predscore, str):
        arr = predscore.replace(')', '').split('(')
        rst = arr[1].strip()
    return rst


def get_pred_from_vep_pathogenicity(predscore):
    rst = predscore
    if isinstance(predscore, str):
        arr = predscore.replace(')', '').split('(')
        rst = arr[0].strip()
    return rst


def add_vep_most_severe(sections):
    most_severe = {}
    most_severe['corder'] = 9999
    for sidx, section in enumerate(sections):
        for consequence in section['Consequence']:
            if (
                (most_severe['corder'] > vep_util.VEP_CONSEQUENCE_ORDER[consequence])
                or (
                    most_severe['corder'] == vep_util.VEP_CONSEQUENCE_ORDER[consequence]
                    and section['CANONICAL'] is not None and section['CANONICAL'] == "1"
                )
            ):
                most_severe['sidx'] = section['Gene']
                most_severe['sidx'] = sidx
                most_severe['canonical'] = section['CANONICAL']
                most_severe['corder'] = vep_util.VEP_CONSEQUENCE_ORDER[consequence]
    for sidx, section in enumerate(sections):
        if sidx == most_severe['sidx']:
            sections[sidx]['MOST_SEVERE'] = '1'
        else:
            sections[sidx]['MOST_SEVERE'] = '0'
    return sections


def vep_select_microannot_add_most_severe(sections, vcf_info_value, select_biotype):
    sections = vep_select_from_microannot(sections, vcf_info_value)
    sections = add_vep_most_severe(sections)
    return sections


def vep_select_biotype_add_most_severe(sections, select_biotype):
    # sections = vep_select_biotype(sections, select_biotype)
    sections = add_vep_most_severe(sections)
    return sections


def vep_select_from_microannot(sections, vcf_info_value):
    selected_feature = {}
    if 'VEP' in vcf_info_value.keys():
        for vcf_info_section in vcf_info_value['VEP']:
            selected_feature[vcf_info_section['Feature']] = vcf_info_section

    selected = []
    for sidx, section in enumerate(sections):
        if section['Feature'] in selected_feature.keys():
            section['Consequence'] = selected_feature[section['Feature']]['Consequence'].split('~')
            selected.append(section)
    return selected


def vep_select_biotype(sections, select_biotype):
    # print (">external_functions.vep_select_biotype():", select_biotype)
    # print(sections)
    selected = []
    for sidx, section in enumerate(sections):
        # print("section['BIOTYPE']:",section['BIOTYPE'])
        if section['BIOTYPE'] in select_biotype:
            selected.append(section)
    return selected


def trim_DIP_ID(v1):
    # ex) DIP-39616N;
    rst = v1
    if isinstance(v1, str):
        rst = v1.split('-')[-1].replace('N', '')
    return rst


def remove_pdb_subversion(id_list):
    rst = id_list
    if isinstance(id_list, str):
        pidmap = {}
        for pid in id_list.split('|'):
            pidmap[pid.split(':')[0]] = 1
        pidlist = list(pidmap.keys())
        rst = '|'.join(pidlist)
    return rst


def convert_rmsk_strand(v1):
    rst = '0'
    if isinstance(v1, str):
        if v1 == '+':
            rst = '1'
    return rst


def convert_NS2blank(v1):
    rst = v1
    if isinstance(v1, str):
        if v1 == 'NS':
            rst = ''
    return rst


def convert_vep_strand(v1):
    rst = v1
    if isinstance(v1, str):
        if v1 == '-1':
            rst = '0'
    return rst


def convert_high_inf_pos(v1):
    rst = v1
    if isinstance(v1, str):
        v1.replace('Y', '1').replace('N', '0')
    return rst


def convert_canonical2boolean(canonical):
    rst = '0'
    # if isinstance(canonical, str):
    if canonical == 'YES':
        rst = '1'
    return rst


def convert_uniprot_transmem(desc_value):
    rst = desc_value
    if isinstance(desc_value, str):
        arr = []
        for v1 in desc_value.split(';'):
            arr.append(v1.strip())
        rst = ";".join(arr)
    return rst


def add_genes_severe_consequence(annotdata, vcfinfo):
    vep_sections = []
    if annotdata is not None and 'VEP' in annotdata.keys():
        vep_sections = annotdata['VEP']
    elif vcfinfo is not None and 'VEP' in vcfinfo.keys():
        vep_sections = vcfinfo['VEP']

    add_item = ['Feature_ncbi', 'HGVSc', 'Amino_acids', 'SIFT_SCORE', 'PolyPhen_SCORE', 'MaxEntScan_diff']

    for attr in vep_sections:
        if attr['MOST_SEVERE'] == '1':
            ensg = attr['Gene']
            most_severe = {}
            most_severe['transcript'] = attr['Feature']
            most_severe['consequence'] = get_most_severe_consequence(attr['Consequence'])
            for ai in add_item:
                most_severe[ai] = attr[ai]
            break
    if len(vep_sections) == 0:
        rst_sections = []
    else:
        d = {}
        d['ensg'] = ensg
        d['most_severe_transcript'] = most_severe['transcript']
        d['most_severe_consequence'] = most_severe['consequence']
        for ai in add_item:
            d['most_severe_' + ai.lower()] = most_severe[ai]

        rst_sections = [d]
    return rst_sections


def get_most_severe_consequence(consequence_list):
    most_severe_order = 100
    most_severe_consequence = ""
    for c1 in consequence_list:
        if most_severe_order > vep_util.VEP_CONSEQUENCE_ORDER[c1]:
            most_severe_order = vep_util.VEP_CONSEQUENCE_ORDER[c1]
            most_severe_consequence = c1
    return most_severe_consequence


def get_value(v1):
    return v1


def tmp_v0_4_8_clinvar_submission(sections, tmp_clinvar_idmap):
    for section in sections:
        section['VariationID'] = tmp_clinvar_idmap[section['ClinVarAccession']]
    return sections
