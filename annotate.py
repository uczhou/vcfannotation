import utils.vcfio as vcfio
import utils.variantquery as variantquery
import pandas
from functools import reduce
import sys


class Load:
    """
    Load class used to create rank of consequence of variants
    """
    def __init__(self):
        filepath = "consequence_of_variants.csv"
        df = pandas.read_csv(filepath)
        self.rank = dict()
        for i in range(df.shape[0]):
            self.rank.update({'_'.join(df['SO term'][i].lower().split(' ')): i})


consequencerank = Load()


def get_consequence(query):
    if query['consequence'] is None:
        return 'null'
    ret = ""
    for key in query['consequence']:
        if len(ret) == 0:
            ret = key
        else:
            if consequencerank.rank.get(key.lower()) < consequencerank.rank.get(ret.lower()):
                ret = key
    return ret


def get_frequency(query):
    """
    Get allele frequency from requested result. Return 0 if not exist.
    :param query: dict type
    :return: float type
    """
    if 'allele_freq' not in query['variant']:
        return 0
    return query['variant']['allele_freq']


def get_info_byid(info, id):
    """
    Parse the info and return value associated with input id
    :param info: vcf INFO
    :param id: vcf INFO ID, case insenstive
    :return: value associated with INFO ID if found, otherwise return None
    """
    info_array = info.split(';')
    for s in info_array:
        ss = s.split('=')
        if ss[0].upper() == id.upper():
            return ss[1]
    return None


def convert_dict_to_str(dict_data):
    """
    Convert dict data into string with ';' as delimiter.
    :param dict_data: dict data
    :return: string representation of dict data
    """
    ss = []
    for kv in dict_data.items():
        if type(kv[1]) is list:
            val = ','.join([str(s) for s in kv[1]])
        else:
            val = str(kv[1])
        ss.append('='.join([kv[0], val]))
    return ';'.join(ss)


def main():
    """
    The main function to start the annotation.
    :return:
    """
    if len(sys.argv) == 1:
        print("Please specify the input file.")
        return
    vcf_filepath = sys.argv[1]
    baseurl = "http://exac.hms.harvard.edu/rest/variant"
    # 1. Open cvf file and read cvf file content into dataframe
    df = vcfio.read_vcf(vcf_filepath)
    annotation = []
    # 2. Extract variants information: CHROM, POS, REF, ALT, QUAL, FILTER, INFO
    for i in range(df.shape[0]):
        ch = df['CHROM'].iloc[i]
        pos = df['POS'].iloc[i]
        ref = df['REF'].iloc[i]
        alt = df['ALT'].iloc[i].split(',')
        # Create annotation information

        # Type of variation. If there are multiple possibilities, most deleterious is recorded.
        consequence = []
        # Allele frequency of variant from Broad Institute ExAC Project API
        frequency = []

        info = df['INFO'].iloc[i].strip()

        # Depth of sequence coverage at the site of variation
        dp = float(get_info_byid(info, 'DP'))

        # Number of reads supporting the variant
        ao = [float(n) for n in get_info_byid(info, 'AO').split(',')]

        # Number of reads supporting the reference variant
        ro = dp - reduce(lambda x, y: x+y, ao)

        # Percentage of reads supporting the variant versus those supporting reference reads
        ao_ro = [float('inf') if ro == 0 else x/ro for x in ao]
        for j in range(len(alt)):
            url = variantquery.createquery(baseurl, ch, pos, ref, alt[j])
            print(url)
            query = variantquery.restquery(url)
            consequence.append(get_consequence(query))
            frequency.append(get_frequency(query))
        info_dict = {'DP': dp, 'AO': ao, 'AO_RO': ao_ro, 'Consequence': consequence, 'ExACFrequency': frequency}
        annotation.append(convert_dict_to_str(info_dict))
    df['INFO'] = annotation

    info_meta = ['##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">',
                 '##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, '
                 'with partial observations recorded fractionally">',
                 '##INFO=<ID=AO_RO,Number=A,Type=Float,Description="Percentage of reads supporting the '
                 'variant versus those supporting reference reads">',
                 '##INFO=<ID=Consequence,Number=A,Type=String,Description="Type of variation. '
                 'If there are multiple possibilities, most deleterious is recorded">',
                 '##INFO=<ID=ExACFrequency,Number=A,Type=Float,Description="Allele frequency of variant '
                 'from Broad Institute ExAC Project API">']

    with open(vcf_filepath) as f:
        other_meta = []
        format_meta = []
        for line in f:
            if line[0] == '#' and line[1] == '#':
                if line[2:8] == 'FORMAT':
                    format_meta.append(line)
                elif line[2:6] != 'INFO':
                    other_meta.append(line)

    vcfio.write_vcf("annotated.vcf", other_meta, info_meta, format_meta, df)


if __name__ == "__main__":
    main()

