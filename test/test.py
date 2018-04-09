import unittest
from ..utils import annotate


class TestAnnotate(unittest.TestCase):
    """
    Test get_consequence method
    """
    def test_get_consequence(self):
        query1 = {"consequence": None}
        val1 = 'null'

        query2 = {"consequence": {'Nearest_gene_five_prime_end': 1, 'frameshift_variant': 2}}
        val2 = 'frameshift_variant'

        query3 = {"consequence": {'transcript_ablation': 1, 'splice_donor_variant': 2}}
        val3 = 'transcript_ablation'

        self.assertEqual(annotate.get_consequence(query1), val1)
        self.assertEqual(annotate.get_consequence(query2), val2)
        self.assertEqual(annotate.get_consequence(query3), val3)

    def test_get_frequency(self):
        query1 = {"variant": {"ref": "G", "xpos": 23155003886, "chrom": "X", "pos": 155003886, "alt": "A"}}
        val1 = 0

        query2 = {"variant": {"ref": "G", "allele_freq": 0.1234,
                              "xpos": 23155003886, "chrom": "X", "pos": 155003886, "alt": "A"}}
        val2 = 0.1234

        self.assertEqual(annotate.get_frequency(query1), val1)
        self.assertEqual(annotate.get_frequency(query2), val2)

    def test_get_info_byid(self):
        info1 = 'AB=0;ABP=0;AC=0;AF=0;AN=6;AO=95;CIGAR=1X;' \
                'DP=4124;DPB=4124;DPRA=0.999031;EPP=9.61615;EPPR=316.776;' \
                'GTI=0;LEN=1;MEANALT=1;MQM=59.7579;MQMR=65.2274;NS=2;NUMALT=1;' \
                'ODDS=591.29;PAIRED=0.989474;PAIREDR=0.966741;PAO=0;PQA=0;PQR=0;' \
                'PRO=0;QA=3774;QR=160284;RO=4029;RPL=51;RPP=4.13032;RPPR=101.278;' \
                'RPR=44;RUN=1;SAF=40;SAP=8.15326;SAR=55;SRF=1663;SRP=269.369;SRR=2366;' \
                'TYPE=snp'
        id1 = 'DP'
        val1 = '4124'

        info2 = 'AB=0;ABP=0;AC=0;AF=0;AN=6;AO=95;CIGAR=1X;' \
                'DP=4124;DPB=4124;DPRA=0.999031;EPP=9.61615;EPPR=316.776;' \
                'GTI=0;LEN=1;MEANALT=1;MQM=59.7579;MQMR=65.2274;NS=2;NUMALT=1;' \
                'ODDS=591.29;PAIRED=0.989474;PAIREDR=0.966741;PAO=0;PQA=0;PQR=0;' \
                'PRO=0;QA=3774;QR=160284;RO=4029;RPL=51;RPP=4.13032;RPPR=101.278;' \
                'RPR=44;RUN=1;SAF=40;SAP=8.15326;SAR=55;SRF=1663;SRP=269.369;SRR=2366;' \
                'TYPE=snp'
        id2 = 'CLA'
        val2 = None

        self.assertEqual(annotate.get_info_byid(info1, id1), val1)
        self.assertEqual(annotate.get_info_byid(info2, id2), val2)

    def test_convert_dict_to_str(self):
        dict_data = {"variant": ['A','T'], "xpos": 23155003886, "chrom": ["X", 1]}
        val = "variant=A,T;xpos=23155003886;chrom=X,1"
        self.assertEqual(annotate.convert_dict_to_str(dict_data), val)


if __name__ == "__main__":
    unittest.main()
