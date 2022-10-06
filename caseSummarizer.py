#!/usr/bin/env python3

'''
--------------------------------------------------------------------------------
caseSummarizer.py: Tool to generate a summary of relevant findings after\
 running CNAG's diagnostics pipeline
--------------------------------------------------------------------------------
'''

import argparse
from docxtpl import DocxTemplate
import os
import pandas as pd
import sys


__version__ = '1.0.0'


def check_input_files(options):
    options.files = {}
    tmp_file_paths = []
    if 'ShortVariants' in options.variant_types:
        options.files['ShortVariants'] = {}
        short_variants_known_pathogenic_path = '{}/ShortVariants/{}.allchr.n\
orm.annot.snpEff.gnomAD.CADD.Clinvar.REVEL.InterVar.CandidateGenes.SpliceAI.t\
agSpliceAI_0.5.P_and_LP.vcf.gz.full.xlsx'.format(options.input_path,
                                                 options.sample)
        options.files['ShortVariants']['known_pathogenic'] =\
            short_variants_known_pathogenic_path
        tmp_file_paths.append(short_variants_known_pathogenic_path)

        short_variants_splicing_path = '{}/ShortVariants/{}.allchr.norm.anno\
t.snpEff.gnomAD.CADD.Clinvar.REVEL.InterVar.CandidateGenes.SpliceAI.filterSpli\
ceAI_0.8.vcf.gz.full.xlsx'.format(options.input_path, options.sample)
        options.files['ShortVariants']['splicing'] =\
            short_variants_splicing_path
        tmp_file_paths.append(short_variants_splicing_path)

    if 'CNV' in options.variant_types:
        options.files['CNV'] = {}
        for window_size in options.freec_windows:
            options.files['CNV'][window_size] = '{}/CNV/{}/{}.{}_Scoresheet.tx\
t'.format(options.input_path, window_size, options.sample, window_size)
        options.files['CNV']['ClinSV'] = '{}/CNV/ClinSV/{}.RARE_PASS_GENE.ligh\
t.CandidateGenes.xlsx'.format(options.input_path, options.sample)
    if 'SV' in options.variant_types:
        sv_path = '{}/SV/{}.SV.intersection.rehead.annot.\
tsv'.format(options.input_path, options.sample)
        tmp_file_paths.append(sv_path)
        options.files['SV'] = sv_path
    if 'MT' in options.variant_types:
        options.files['MT'] = {
            'MToolbox': '{}/MT/Merged_{}_MT.gnomAD.MitoMAP_Disease.MitoMap_pol\
ymorphisms.MitoTIP.vcf.gz.full.xlsx'.format(options.input_path, options.sample),
            'GATK': '{}/MT/{}.MT_GATK.gnomAD.MitoMAP_Disease.MitoMap_poly\
morphisms.MitoTIP.vcf.gz.full.xlsx'.format(options.input_path, options.sample),
            }
    if 'ROH' in options.variant_types:
        roh_path = '{}/ROH/{}.hom'.format(options.input_path,
                                          options.sample)
        options.files['ROH'] = roh_path
        tmp_file_paths.append(roh_path)
    if 'CYP2D6' in options.variant_types:
        cyp2d6_path = '{}/CYP2D6/{}_CYP2D6.tsv'.format(options.input_path,
                                                       options.sample)
        options.files['CYP2D6'] = cyp2d6_path
        tmp_file_paths.append(cyp2d6_path)
    if 'SMA' in options.variant_types:
        sma_path = '{}/SMA/{}_SMA.tsv'.format(options.input_path,
                                              options.sample)
        options.files['SMA'] = sma_path
        tmp_file_paths.append(sma_path)
    if 'STR' in options.variant_types:
        str_path = '{}/STR/{}_EH4_STRs.tsv'.format(options.input_path,
                                                   options.sample)
        options.files['STR'] = str_path
        tmp_file_paths.append(str_path)
    if 'PGx' in options.variant_types:
        pgx_path = '{}/PGx/{}.merged.stargazer-genotype.\
txt'.format(options.input_path, options.sample)
        options.files['PGx'] = pgx_path
        tmp_file_paths.append(pgx_path)
    if 'MHC' in options.variant_types:
        mhc_path = '{}/MHC/{}_final.result.txt'.format(options.input_path,
                                                       options.sample)
        options.files['MHC'] = mhc_path
        tmp_file_paths.append(mhc_path)
    if 'GBA' in options.variant_types:
        gba_path = '{}/GBA/{}_GBA.tsv'.format(options.input_path,
                                              options.sample)
        options.files['GBA'] = gba_path
        tmp_file_paths.append(gba_path)

    for tfp in tmp_file_paths:
        if not os.path.exists(tfp):
            raise ValueError('ERROR: File {} does not exist'.format(tfp))


def cnv_sniffer_clinsv(clinsv_file, sample, candidate_genes,
                       output_selected_variants=False,
                       selected_variants_folder=None):
    reportable_str = []
    cnv_counter = {}
    for gene_list in candidate_genes:
        cnv_counter[gene_list] = 0
    all_clinsv_pd = pd.read_excel(clinsv_file)
    sample_clinsv_pd = all_clinsv_pd[all_clinsv_pd['SAMPLE'] == sample]
    if len(sample_clinsv_pd['SAMPLE']) == 0:
        raise ValueError('ERROR: The file {} does not contain ClinSV analyses\
 for sample {}'.format(clinsv_file, sample))

    reportable_clinsv_pd = pd.DataFrame(columns=sample_clinsv_pd.columns)

    for gene_list in candidate_genes:
        gene_list_tag = gene_list + '_Tag'
        tmp_pd =\
            sample_clinsv_pd[
                (sample_clinsv_pd['GFEAT'].isin(('exon', 'CDS',
                                                 'start codon'))) &
                (sample_clinsv_pd['SU'] > 0) &
                (~sample_clinsv_pd[gene_list_tag].isna())]
        reportable_clinsv_pd =\
            pd.concat([reportable_clinsv_pd, tmp_pd]).drop_duplicates()

    if len(reportable_clinsv_pd) == 0:
        reportable_str.append('No significant results')
    else:
        for gene_list in candidate_genes:
            gene_list_tag = gene_list + '_Tag'
            tmp_pd = reportable_clinsv_pd[
                ~reportable_clinsv_pd[gene_list_tag].isna()]
            cnv_counter[gene_list] = len(tmp_pd)
        reportable_str.append('{}: {}'.format(gene_list,
                                              cnv_counter[gene_list]))
    returnable_str = ',\n'.join(reportable_str)

    if output_selected_variants is True and\
       len(reportable_clinsv_pd) > 0:
        save_selected_variants(selected_variants_folder,
                               '{}_ClinSV_Autoselected.xlsx'.format(sample),
                               reportable_clinsv_pd)

    return ['CNV', 'CNV', 'ClinSV', returnable_str, '', '']

def cnv_sniffer_freec_simple(freec_file, sample, window_size,
                       output_selected_variants=False,
                       selected_variants_folder=None):
    reportable_str = []
    cnv_counter = {}
    freec_pd = pd.read_csv(freec_file, delimiter='\t')

    pathogenic_pd =\
        freec_pd[(freec_pd['Classification'].str.contains('Pathogenic'))]
    likely_pathogenic_pd =\
        freec_pd[freec_pd['Classification'].str.contains('Likely pathogenic')]

    reportable_freec_pd =\
        pd.concat([pathogenic_pd, likely_pathogenic_pd]).drop_duplicates()

    if len(reportable_freec_pd) == 0:
        reportable_str.append('No significant results')
    else:
        reportable_str.append(
            '{} Pathogenic, {} Likely pathogenic'.format(
                len(pathogenic_pd),
                len(likely_pathogenic_pd)))
    returnable_str = ',\n'.join(reportable_str)

    if output_selected_variants is True and\
       len(reportable_freec_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_Control-Freec_{}_Autoselected.xlsx'.format(sample, window_size),
            reportable_freec_pd)

    return ['CNV', 'CNV', 'Control-FREEC (window size: {})'.format(window_size),
             returnable_str, '', '']


def cnv_sniffer_freec_complex(freec_file, sample, candidate_genes, window_size,
                       output_selected_variants=False,
                       selected_variants_folder=None):
    reportable_str = []
    cnv_counter = {}
    freec_pd = pd.read_csv(freec_file, delimiter='\t')
    freec_['All_protein_coding_genes'] =\
        supu['All protein coding genes'].str.split(', ')
    reportable_freec_pd = pd.DataFrame(columns=freec_pd.columns)

    for gene_list in candidate_genes:
        cnv_counter[gene_list] = 0

        pathogenic_pd =\
            freec_pd[(freec_pd['Classification'].str.contains('Pathogenic'))]
        likely_pathogenic_pd =\
            freec_pd[freec_pd['Classification'].str.contains('Likely pathog\
enic')]
    
        reportable_freec_pd =\
            pd.concat([reportable_freec_pd, pathogenic_pd,
                       likely_pathogenic_pd]).drop_duplicates()

    if len(reportable_freec_pd) == 0:
        reportable_str.append('No significant results')
    else:
        reportable_str = 'No significant results'
        for gene_list in candidate_genes:
            gene_list_tag = gene_list + '_Tag'
            tmp_pd = reportable_clinsv_pd[
                ~reportable_clinsv_pd[gene_list_tag].isna()]
            cnv_counter[gene_list] = len(tmp_pd)
        reportable_str.append('{}: {}'.format(gene_list,
                                              cnv_counter[gene_list]))
    returnable_str = ',\n'.join(reportable_str)

    if output_selected_variants is True and\
       len(reportable_clinsv_pd) > 0:
        save_selected_variants(selected_variants_folder,
                               '{}_ClinSV_Autoselected.xlsx'.format(sample),
                               reportable_clinsv_pd)

    return ['CNV', 'CNV', 'Control-FREEC (window size: {})'.format(window_size),
             returnable_str, '', '']


def cyp2d6_sniffer(cyp2d6_file, sample):
    all_cyp2d6_pd = pd.read_csv(cyp2d6_file, delimiter='\t')
    sample_cyp2d6_pd = all_cyp2d6_pd[all_cyp2d6_pd['Sample'] == sample]
    if len(sample_cyp2d6_pd['Sample']) == 0:
        raise ValueError('ERROR: The file {} does not contain CYP2D6 analyses\
 for sample {}'.format(cyp2d6_file, sample))
    if len(sample_cyp2d6_pd['Sample']) > 1:
        raise ValueError('ERROR: The file {} contains more than 1 CYP2D6\
 analysis for sample {}'.format(cyp2d6_file, sample))

    returnable_str = 'Genotype: {}\nFilter: {}'.format(
        str(sample_cyp2d6_pd.at[0, 'Genotype']),
        str(sample_cyp2d6_pd.at[0, 'Filter']))

    return ['Pseudogenes', 'CYP2D6', 'CYP2D6', returnable_str, '', '']


def gba_sniffer(gba_file, sample):
    all_gba_pd = pd.read_csv(gba_file, delimiter='\t')
    sample_gba_pd = all_gba_pd[all_gba_pd['Sample'] == sample]
    if len(sample_gba_pd['Sample']) == 0:
        raise ValueError('ERROR: The file {} does not contain GBA analyses for\
 sample {}'.format(gba_file, sample))
    if len(sample_gba_pd['Sample']) > 1:
        raise ValueError('ERROR: The file {} contains more than 1 GBA analysis\
 for sample {}'.format(gba_file, sample))
    
    is_biallelic = str(sample_gba_pd.at[0, 'is_biallelic_GBAP1-like_variant_e\
xon9-11'])
    is_carrier = str(sample_gba_pd.at[0, 'is_carrier_GBAP1-like_variant_exon\
9-11'])
    total_cn = str(sample_gba_pd.at[0, 'total_CN'])
    deletion_breakpoint = str(sample_gba_pd.at[0, 'deletion_breakpoint_in_GB\
A_gene'])
    gbap1_like_variant = str(sample_gba_pd.at[0, 'GBAP1-like_variant_exon9-11'])
    other_variants = str(sample_gba_pd.at[0, 'other_variants'])

    if (is_biallelic != 'False') or\
       (is_carrier != 'False') or\
       (total_cn != '4') or\
       (deletion_breakpoint != 'None') or\
       (gbap1_like_variant != 'None') or\
       (other_variants != 'None'):
        returnable_str = '''Is biallelic GBAP1-like variant exon9-11: {}
Is carrier GBAP1-like variant exon9-11: {}
Total CN: {}
Deletion breakpoint in GBA gene: {}
GBAP1-like variant exon9-11: {}
Other variants: {}'''.format(is_biallelic, is_carrier, total_cn,
                             deletion_breakpoint, gbap1_like_variant,
                             other_variants)
    else:
        returnable_str = 'No significant results'

    return ['Pseudogenes', 'GBA', 'GBA', returnable_str, '', '']


def get_cs_options(argv):
    cs_parser = argparse.ArgumentParser(
        description='Tool to generate a summary of relevant findings after\
 running CNAG\'s diagnostics pipeline',
        usage='caseSummarizer.py [options]',
    )
    mandatory_csp = cs_parser.add_argument_group('Mandatory arguments')
    mandatory_csp.add_argument(
        '--input_path',
        '-I',
        help='PATH Folder that stores the results',
        required=True,
    )
    mandatory_csp.add_argument(
        '--output',
        '-O',
        default='-',
        help='FILE Output filename',
        required=True,
    )
    mandatory_csp.add_argument(
        '--candidate_genes_conf',
        '-C',
        help='PATH Config file with the lists of candidate genes',
        required=True,
    )
    cs_parser.add_argument(
        '--reportable_variants_conf',
        '-R',
        default='/caseSummarizer-1.0/DefaultReportable.conf',
        help='PATH Config file with the types of variants to report',
    )
    cs_parser.add_argument(
        '--dp_threshold',
        '-DP',
        type=int,
        default=1,
        help='INT DP threshold to accept a variant',
    )
    cs_parser.add_argument(
        '--gq_threshold',
        '-GQ',
        type=int,
        default=1,
        help='INT GQ threshold to accept a variant. ShortVariants only',
    )
    cs_parser.add_argument(
        '--freec_windows',
        '-FW',
        type=int,
        default=[5000, 20000, 50000],
        help='INT List of window sizes used by Control-FREEC',
    )
    cs_parser.add_argument(
        '--report_template',
        '-RT',
        default='/caseSummarizer-1.0/case_summarizer_default.docx',
        help='FILE Path to the report template',
    )
    cs_parser.add_argument(
        '--dont_output_candidate_variants',
        '-DOCV',
        action='store_true',
        default=False,
        help='BOOL Do not generate a folder with the candidate variants',
    )

    return cs_parser.parse_args(argv)


def load_candidate_genes(options):
    options.candidate_genes = {}
    with open(options.candidate_genes_conf) as cand_genes_conf:
        for cand_genes_list in cand_genes_conf:
            cand_gene_list_path, tag = cand_genes_list.rstrip().split('\t')
            with open(options.candidate_genes_path + '/' + cand_gene_list_path)\
             as cand_gene_list:
                tmp_list = cand_gene_list.readlines()
                cand_list = [x.rstrip() for x in tmp_list if x.startswith('#')
                             is False]
                options.candidate_genes[tag] = cand_list


def make_summary_report(sample, results_records, output_file, template):
    doc = DocxTemplate(template)
    context = {'sample' : sample, 'results_records': results_records}
    doc.render(context)
    doc.save(output_file)


def mt_sniffer_gatk(mt_file_gatk, sample, output_selected_variants=False,
                    selected_variants_folder=None):
    mt_variants_pd = pd.read_excel(mt_file_gatk)

    confirmed_mt_variants_pd =\
        mt_variants_pd[(~mt_variants_pd['PASS'].isna()) &
                       (mt_variants_pd[(sample + '_AF')] >= 0.01) &
                       (mt_variants_pd['Mitomap_Status'] == 'Cfrm')]

    reported_mt_variants_pd =\
        mt_variants_pd[(~mt_variants_pd['PASS'].isna()) &
                       (mt_variants_pd[(sample + '_AF')] >= 0.01) &
                       (mt_variants_pd['Mitomap_Status'] == 'Reported') &
                       (mt_variants_pd['MitoMap_AF'] < 0.5)]
    
    reportable_mt_variants_pd =\
        pd.concat([confirmed_mt_variants_pd,
                   reported_mt_variants_pd]).drop_duplicates()

    returnable_str = '{} Confirmed, {} Reported'.format(
        len(confirmed_mt_variants_pd), len(reported_mt_variants_pd))

    if output_selected_variants is True and\
       len(reportable_mt_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_MT_GATK_Autoselected.xlsx'.format(sample),
            reportable_mt_variants_pd)

    return ['mtDNA', 'MT', 'Known pathogenic (GATK)', returnable_str, '', '']


def mt_sniffer_mtb(mt_file_mtb, sample, output_selected_variants=False,
                    selected_variants_folder=None):

    mt_variants_pd = pd.read_excel(mt_file_mtb)

    confirmed_mt_variants_pd =\
        mt_variants_pd[(~mt_variants_pd['PASS'].isna()) &
                       (mt_variants_pd[(sample + '_AF')] >= 0.01) &
                       (mt_variants_pd['Mitomap_Status'] == 'Cfrm')]

    reported_mt_variants_pd =\
        mt_variants_pd[(~mt_variants_pd['PASS'].isna()) &
                       (mt_variants_pd[(sample + '_AF')] >= 0.01) &
                       (mt_variants_pd['Mitomap_Status'] == 'Reported') &
                       (mt_variants_pd['MitoMap_AF'] < 0.5)]
    
    reportable_mt_variants_pd =\
        pd.concat([confirmed_mt_variants_pd,
                   reported_mt_variants_pd]).drop_duplicates()

    returnable_str = '{} Confirmed, {} Reported'.format(
        len(confirmed_mt_variants_pd), len(reported_mt_variants_pd))

    if output_selected_variants is True and\
       len(reportable_mt_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_MT_MToolbox_Autoselected.xlsx'.format(sample),
            reportable_mt_variants_pd)

    return ['mtDNA', 'MT', 'Known pathogenic (MToolBox)', returnable_str,
            '', '']


def no_flags():
    print('Please re-run the command with "-h" to get usage instructions and a\
 complete list of options\n')
    exit()


def roh_consanguinity_sniffer(roh_file, sample):
    all_roh_pd = pd.read_csv(roh_file, delimiter='\s+')
    sample_roh_pd = all_roh_pd[all_roh_pd['IID'] == sample]
    if len(sample_roh_pd['IID']) == 0:
        raise ValueError('ERROR: The file {} does not contain ROH analyses for\
 sample {}'.format(roh_file, sample))

    tmp_pd = sample_roh_pd[(sample_roh_pd['CHR'] != 23) &
                           (sample_roh_pd['CHR'] != 'X')]
    total_size = tmp_pd['KB'].sum()

    if total_size < 22000:
        roh_classification = 'Non consanguineous'
    elif 22000 <= total_size < 79000:
        roh_classification = 'Uncertain'
    elif 79000 <= total_size < 123000:
        roh_classification = 'Probably consanguineous'
    elif total_size >= 123000:
        roh_classification = 'Consanguineous'

    return ['ROH', 'ROH', 'Consanguinity',
            '{} ({} Mb)'.format(roh_classification, round(total_size/1000, 3)),
            '',
            '']


def roh_uniparental_disomy_sniffer(roh_file, sample):
    big_chunks_catalog = {}
    uniparental_disomy_status = 'Uniparental disomy NOT detected'

    all_roh_pd = pd.read_csv(roh_file, delimiter='\s+')
    sample_roh_pd = all_roh_pd[all_roh_pd['IID'] == sample]
    if len(sample_roh_pd['IID']) == 0:
        raise ValueError('ERROR: The file {} does not contain ROH analyses for\
 sample {}'.format(roh_file, sample))

    for chromosome in sample_roh_pd['CHR'].unique():
        if chromosome != 'X' and chromosome != 23:
            tmp_pd = sample_roh_pd[(sample_roh_pd['CHR'] == chromosome) &
                                   (sample_roh_pd['KB'] >= 5000)]
            if len(tmp_pd) > 0:
                big_chunks_catalog[chromosome] = tmp_pd['KB']
    if len(big_chunks_catalog.keys()) > 0:
        for chrom in big_chunks_catalog:
            if sum(big_chunks_catalog[chrom]) > 10000:
                uniparental_disomy_status = 'Uniparental disomy detected'
    return ['ROH', 'ROH', 'Uniparental disomy', uniparental_disomy_status, '',
            '']


def save_selected_variants(target_folder, target_file, pd_result):
    out_file = target_folder + '/' + target_file
    if out_file.endswith('.xlsx'):
        pd_result.to_excel(out_file, index=False)
    elif out_file.endswith('.tsv'):
        pd_result.to_csv(out_file, sep='\t', index=False)
    else:
        raise ValueError('ERROR: Unknown filetype')


def short_variants_sniffer_known_pathogenic(
        short_variants_file, sample, candidate_genes, dp_threshold=1,
        gq_threshold=1, output_selected_variants=False,
        selected_variants_folder=None):
    short_variants_pd = pd.read_excel(short_variants_file)
    sample_dp = sample + '_DP'
    sample_gq = sample + '_GQ'

    reportable_short_variants_pd = pd.DataFrame(
        columns=short_variants_pd.columns)
    returnable_list = []
    for gene_list in candidate_genes:
        pathogenic_pd =\
            short_variants_pd[(short_variants_pd['InterVar_ClinicalSignifican\
ce'].str.contains('Pathogenic')) & (~short_variants_pd[gene_list].isna()) &
            (short_variants_pd[sample_dp] >= dp_threshold) &
            (short_variants_pd[sample_gq] >= gq_threshold)]
        likely_pathogenic_pd =\
            short_variants_pd[(short_variants_pd['InterVar_ClinicalSignificanc\
e'].str.contains('Likely_pathogenic')) &
            (~short_variants_pd[gene_list].isna()) &
            (short_variants_pd[sample_dp] >= dp_threshold) &
            (short_variants_pd[sample_gq] >= gq_threshold)]
        reportable_short_variants_pd =\
            pd.concat([reportable_short_variants_pd, pathogenic_pd,
                       likely_pathogenic_pd]).drop_duplicates()
        returnable_list.append(
            '{}: {} Pathogenic, {} Likely pathogenic'.format(
                gene_list,
                len(pathogenic_pd),
                len(likely_pathogenic_pd)))

    if output_selected_variants is True and\
       len(reportable_short_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_ShortVariants_Known_pathogenic_Autoselected.xl\
sx'.format(sample),
            reportable_short_variants_pd)

    return ['SNV-InDels', 'ShortVariants', 'Known pathogenic',
            ',\n'.join(returnable_list), '', '']


def short_variants_sniffer_splicing(
        short_variants_file, sample, candidate_genes, dp_threshold=1,
        gq_threshold=1, output_selected_variants=False,
        selected_variants_folder=None):
    short_variants_pd = pd.read_excel(short_variants_file)

    sample_dp = sample + '_DP'
    sample_gq = sample + '_GQ'

    reportable_short_variants_pd = pd.DataFrame(
        columns=short_variants_pd.columns)
    tmp_pd = {}
    returnable_list = []
    for gene_list in candidate_genes:
        tmp_pd[gene_list] =\
            short_variants_pd[(short_variants_pd['gnomAD_WG_AF'] < 0.01) &
                              (~short_variants_pd[gene_list].isna()) &
                              (short_variants_pd[sample_dp] >= dp_threshold) &
                              (short_variants_pd[sample_gq] >= gq_threshold)]
        reportable_short_variants_pd =\
            pd.concat([reportable_short_variants_pd,
                       tmp_pd[gene_list]]).drop_duplicates()
        returnable_list.append('{}: {}'.format(gene_list,
                                               len(tmp_pd[gene_list])))
    if output_selected_variants is True and\
       len(reportable_short_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_ShortVariants_Splicing_Autoselected.xlsx'.format(sample),
            reportable_short_variants_pd)

    return ['SNV-InDels', 'ShortVariants', 'Splicing',
            ',\n'.join(returnable_list), '', '']


def sma_sniffer(sma_file, sample, output_selected_variants=False,
                selected_variants_folder=None):
    all_sma_pd = pd.read_csv(sma_file, delimiter='\t')
    sample_sma_pd = all_sma_pd[all_sma_pd['Sample'] == sample]
    if len(sample_sma_pd['Sample']) == 0:
        raise ValueError('ERROR: The file {} does not contain SMA analyses\
 for sample {}'.format(sma_file, sample))
    if len(sample_sma_pd['Sample']) > 1:
        raise ValueError('ERROR: The file {} contains more than 1 SMA\
 analysis for sample {}'.format(sma_file, sample))

    is_sma = str(sample_sma_pd.at[0, 'isSMA'])
    is_carrier = str(sample_sma_pd.at[0, 'isCarrier'])

    reportable_sma_pd = sample_sma_pd[sample_sma_pd['isSMA'] |
                                      sample_sma_pd['isCarrier']]

    return ['Pseudogenes', 'SMA', 'SMA',
            'Is SMA: {}\nIs Carrier: {}'.format(is_sma, is_carrier), '', '']


def str_sniffer(str_file, sample, output_selected_variants=False,
                selected_variants_folder=None):
    all_str_pd = pd.read_csv(str_file, delimiter='\t')
    sample_str_pd = all_str_pd[all_str_pd['SAMPLE'] == sample]
    reportable_str = []
    reportable_str_pd = sample_str_pd[sample_str_pd['SO'] !=
                                      'SPANNING/SPANNING']
    if len(reportable_str_pd) == 0:
        reportable_str.append('No significant results')
    else:
        reportable_str = list(reportable_str_pd['INFO_REPID'])
    returnable_str = ',\n'.join(reportable_str)

    if output_selected_variants is True:
        save_selected_variants(selected_variants_folder,
                               '{}_STR_Autoselected.tsv'.format(sample),
                               reportable_str_pd)

    return ['STR', 'STR', 'Known pathogenic', returnable_str, '', '']


def summarize_dx_results(options):
    res_records = []

    # ShortVariants
    if 'ShortVariants' in options.variant_types:
        res_records.append(short_variants_sniffer_known_pathogenic(
            options.files['ShortVariants']['known_pathogenic'],
            options.sample,
            options.candidate_genes,
            output_selected_variants=options.output_candidate_variants,
            selected_variants_folder=options.auto_selected_variants_folder,
            dp_threshold=options.dp_threshold,
            gq_threshold=options.gq_threshold))

        res_records.append(short_variants_sniffer_splicing(
            options.files['ShortVariants']['splicing'],
            options.sample,
            options.candidate_genes,
            output_selected_variants=options.output_candidate_variants,
            selected_variants_folder=options.auto_selected_variants_folder,
            dp_threshold=options.dp_threshold,
            gq_threshold=options.gq_threshold))

    # MT
    if 'MT' in options.variant_types:
        res_records.append(mt_sniffer_gatk(
            options.files['MT']['GATK'],
            options.sample,
            output_selected_variants=options.output_candidate_variants,
            selected_variants_folder=options.auto_selected_variants_folder))
    
        # res_records.append(mt_sniffer_gatk(
        #     options.files['MT']['MToolbox'],
        #     options.sample,
        #     output_selected_variants=options.output_candidate_variants,
        #     selected_variants_folder=options.auto_selected_variants_folder))

    # CNV
    if 'CNV' in options.variant_types:
        res_records.append(
            cnv_sniffer_clinsv(
                options.files['CNV']['ClinSV'],
                options.sample,
                options.candidate_genes,
                output_selected_variants=options.output_candidate_variants,
                selected_variants_folder=options.auto_selected_variants_folder))

        for window_size in options.freec_windows:
            res_records.append(
                cnv_sniffer_freec_simple(
                    options.files['CNV'][window_size],
                    options.sample,
                    window_size=window_size,
                    output_selected_variants=options.output_candidate_variants,
                    selected_variants_folder=\
                        options.auto_selected_variants_folder))

    
    # SV
    if 'SV' in options.variant_types:
        res_records.append(
            sv_sniffer_survivor_simple(
                options.files['SV'],
                options.sample,
                output_selected_variants=options.output_candidate_variants,
                selected_variants_folder=options.auto_selected_variants_folder))

    # STR
    if 'STR' in options.variant_types:
        res_records.append(str_sniffer(
            options.files['STR'],
            options.sample,
            output_selected_variants=options.output_candidate_variants,
            selected_variants_folder=options.auto_selected_variants_folder))

    # GBA
    if 'GBA' in options.variant_types:
        res_records.append(gba_sniffer(options.files['GBA'],
                                       options.sample))

    # SMA
    if 'SMA' in options.variant_types:
        res_records.append(sma_sniffer(
            options.files['SMA'],
            options.sample,
            output_selected_variants=options.output_candidate_variants,
            selected_variants_folder=options.auto_selected_variants_folder))

    # CYP2D6
    if 'CYP2D6' in options.variant_types:
        res_records.append(
            cyp2d6_sniffer(options.files['CYP2D6'], options.sample))

    # ROH
    if 'ROH' in options.variant_types:
        res_records.append(
            roh_consanguinity_sniffer(options.files['ROH'],
                                      options.sample))
        res_records.append(
            roh_uniparental_disomy_sniffer(options.files['ROH'],
                                           options.sample))

    # PGx

    make_summary_report(options.sample, res_records, options.output,
                        options.report_template)


def sv_sniffer_survivor_simple(structural_variants_file, sample,
                        output_selected_variants=False,
                        selected_variants_folder=None):

    structural_variants_pd = pd.read_csv(structural_variants_file,
                                         delimiter='\t')

    returnable_list = []
    pathogenic_pd =\
        structural_variants_pd[structural_variants_pd['ACMG_class'] == 'full=5']
    likely_pathogenic_pd =\
        structural_variants_pd[structural_variants_pd['ACMG_class'] == 'full=4']
    reportable_structural_variants_pd = pd.concat(
        [pathogenic_pd, likely_pathogenic_pd]).drop_duplicates()
    returnable_list.append(
        '{} Pathogenic, {} Likely pathogenic'.format(
            len(pathogenic_pd),
            len(likely_pathogenic_pd)))

    if output_selected_variants is True and\
       len(reportable_structural_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_Structural_Variants_Known_pathogenic_Autoselected.xl\
sx'.format(sample),
            reportable_structural_variants_pd)

    return ['SV', 'SV', 'Known pathogenic',
            ',\n'.join(returnable_list), '', '']


def sv_sniffer_survivor_complex(structural_variants_file, sample,
                                candidate_genes, output_selected_variants=False,
                        selected_variants_folder=None):
    structural_variants_pd = pd.read_csv(structural_variants_file,
                                         delimiter='\t')

    reportable_structural_variants_pd = pd.DataFrame(
        columns=structural_variants_pd.columns)
    returnable_list = []
    for gene_list in candidate_genes:
        pathogenic_pd =\
            structural_variants_pd[(structural_variants_pd['InterVar_Clinical\
Significance'].str.contains('Pathogenic')) &
                              (~short_variants_pd[gene_list].isna())]
        likely_pathogenic_pd =\
            short_variants_pd[(short_variants_pd['InterVar_Clinica\
lSignificance'].str.contains('Likely_pathogenic')) &
                              (~short_variants_pd[gene_list].isna())]
        reportable_short_variants_pd =\
            pd.concat([reportable_short_variants_pd, pathogenic_pd,
                       likely_pathogenic_pd]).drop_duplicates()
        returnable_list.append(
            '{}: {} Pathogenic, {} Likely pathogenic'.format(
                gene_list,
                len(pathogenic_pd),
                len(likely_pathogenic_pd)))

    if output_selected_variants is True and\
       len(reportable_short_variants_pd) > 0:
        save_selected_variants(
            selected_variants_folder,
            '{}_ShortVariants_Known_pathogenic_Autoselected.xls\
x'.format(sample),
            reportable_short_variants_pd)

    return ['SV', 'SV', 'Known pathogenic',
            ',\n'.join(returnable_list), '', '']


def validate_options(options):
    options.input_path = os.path.abspath(options.input_path)
    path_list = options.input_path.split('/')
    options.sample = path_list[-1]
    options.pedigree = path_list[-2]
    options.candidate_genes_conf = os.path.abspath(options.candidate_genes_conf)
    options.candidate_genes_path = os.path.dirname(options.candidate_genes_conf)
    load_candidate_genes(options)
    validate_reportable_results(options)
    check_input_files(options)
    auto_selected_variants_folder =\
        os.getcwd() + '/{}_AutoSelectedVariants'.format(options.sample)
    options.auto_selected_variants_folder = auto_selected_variants_folder
    if options.dont_output_candidate_variants is False:
        options.output_candidate_variants = True
        os.makedirs(auto_selected_variants_folder)
    else:
        options.output_candidate_variants = False


def validate_reportable_results(options):
    known_variant_types = ['ShortVariants', 'CNV', 'SV', 'MT', 'ROH', 'STR',
                           'PGx', 'CYP2D6', 'GBA', 'SMA']
    with open(options.reportable_variants_conf) as reportable_variants_conf:
        tmp_variant_types = reportable_variants_conf.read().splitlines()
        if len(tmp_variant_types) == 0:
            raise ValueError('ERROR: At least a valid type of variant must be\
 requested. Known types of variants are {}'.format(', '.join(known_variant_types)))
        for variant_type in tmp_variant_types:
            if variant_type not in known_variant_types:
                raise ValueError('ERROR: {} is an unrecognised variant type.\
 Known types are {}'.format(variant_type, ', '.join(known_variant_types)))
        options.variant_types = tmp_variant_types


def main():
    if len(sys.argv) == 1:
        no_flags()
    elif len(sys.argv) > 1:
        cs_opts = get_cs_options(sys.argv[1:])
        validate_options(cs_opts)
        summarize_dx_results(cs_opts)
        sys.exit(0)


if __name__ == '__main__':
    exit(main())
