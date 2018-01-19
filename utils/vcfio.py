import pandas
import re


def read_vcf(filepath):
    """
    Read vcf file variants information, return as pandas DataFrame object.
    :param filepath:
    :return: pandas DataFrame object
    """

    with open(filepath) as f:
        variant_info = []
        labels = []
        for line in f:
            if line[0] == '#' and line[1] != '#':
                labels = re.split(r'[\t\s]\s*', line[1:])
            elif line[0] != '#':
                fields = re.split(r'[\t\s]\s*', line)
                variant_info.append(tuple(fields))

    df = pandas.DataFrame.from_records(variant_info, columns=labels)
    return df


def read_vcf_header(filepath):
    """
    Read vcf file header information, return as dict.
    :param filepath: file path
    :return: header information with list
    """
    vcf_header = []
    with open(filepath) as f:
        for line in f:
            if line[0] == '#' and line[1] == '#':
                vcf_header.append(line)
    return vcf_header


def write_vcf(filename, other_meta, info_meta, format_meta, df):
    """
    Write file into filename in vcf format.
    :param filename:
    :param other_meta: meta data such as fileformat, filedate, source, reference and contig.
    :param info_meta: INFO meta data
    :param format_meta: FORMAT meta data
    :param df: pandas DataFrame contains detailed annotated variant information
    :return:
    """
    names = list(df)
    names[0] = "{}{}".format('#', names[0])
    name = '\t'.join(names)
    with open(filename, 'w') as f:
        for line in other_meta:
            f.write(line)
        for line in info_meta:
            f.write(line+'\n')
        for line in format_meta:
            f.write(line)
        f.write(name+'\n')
        for i in range(df.shape[0]):
            f.write('\t'.join([str(x) for x in df.loc[i]])+'\n')


def getkvpair(line):
    """
    Convert string to key, value pair
    :param line: string
    :return: tuple
    Example:
        input: 'INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">'
        output: (INFO, 'ID=NS,Number=1,Type=Integer,Description="Number of samples with data"')
    """
    if '<' in line:
        ss = line.split('<')
        key = ss[0]
        value = ss[1][0:len(ss[1])]
    else:
        ss = line.split('=')
        key = ss[0]
        value = ss[1]
    return key, value

