import sys
from io import StringIO
import unittest
sys.path.append(".")
from libs.readmapperstats import ReadMapperStats

class TestReadMapperStats(unittest.TestCase):

    def setUp(self):
        self.read_mapper_stats = ReadMapperStats()
        self.example_data = ExampleData()
    
    def test_count_fasta_fh_entries(self):
        fasta_fh_1 = StringIO(self.example_data.fasta_seqs_1)
        self.assertEqual(
            3, 
            self.read_mapper_stats._count_fasta_fh_entries(fasta_fh_1))
        fasta_fh_2 = StringIO(self.example_data.fasta_seqs_2)
        self.assertEqual(
            5, 
            self.read_mapper_stats._count_fasta_fh_entries(fasta_fh_2))

    def test_count_fasta_entry_set(self):
        test_dict = {}
        file_names = ["foo", "bar", "space"]
        file_paths = ["mope/foo", "mope/bar", "mope/space"]
        self.read_mapper_stats._count_fasta_entries = _mock_count_fasta_entries
        self.read_mapper_stats._count_fasta_entry_set(
            file_names, file_paths, test_dict)
        self.assertDictEqual(
            {"bar": "the number of entries in mope/bar",
             "foo": "the number of entries in mope/foo",
             "space": "the number of entries in mope/space"}, test_dict)
        
    def test_count_raw_reads(self):
        file_names = ["lo", "bo", "ro"]
        file_paths = ["lala/lo", "lala/bo", "lala/ro"]
        self.read_mapper_stats._count_fasta_entries = _mock_count_fasta_entries
        self.read_mapper_stats.count_raw_reads(file_names, file_paths)
        self.assertDictEqual(
            {"lo": "the number of entries in lala/lo",
             "bo": "the number of entries in lala/bo",
             "ro": "the number of entries in lala/ro"},
            self.read_mapper_stats.raw_read_countings)

    def test_too_small_clipped_reads(self):
        file_names = ["no", "to", "so"]
        file_paths = ["fufu/no", "fufu/to", "fufu/so"]
        self.read_mapper_stats._count_fasta_entries = _mock_count_fasta_entries
        self.read_mapper_stats.count_too_small_clipped_reads(
            file_names, file_paths)
        self.assertDictEqual(
            {"no": "the number of entries in fufu/no",
             "to": "the number of entries in fufu/to",
             "so": "the number of entries in fufu/so"},
            self.read_mapper_stats.too_small_clipped_reads)

    def test_count_long_enough_clipped_reads(self):
        file_names = ["tick", "trick", "track"]
        file_paths = ["ducks/tick", "ducks/trick", "ducks/track"]
        self.read_mapper_stats._count_fasta_entries = _mock_count_fasta_entries
        self.read_mapper_stats.count_long_enough_clipped_reads(
            file_names, file_paths)
        self.assertDictEqual(
            {"tick": "the number of entries in ducks/tick",
             "trick": "the number of entries in ducks/trick",
             "track": "the number of entries in ducks/track"},
            self.read_mapper_stats.long_enough_clipped_reads)

    def test_write_stats_to_fh(self):
        read_file_name = ["xxx", "yyy", "zzz"]
        self.read_mapper_stats.raw_read_countings = {
            "zzz" : 5, "yyy" : 1, "xxx" : 8}
        self.read_mapper_stats.long_enough_clipped_reads = {
            "zzz" : 42, "yyy" : 23, "xxx" : 5}
        self.read_mapper_stats.too_small_clipped_reads = {
            "zzz" : 0, "yyy" : 1, "xxx" : 2}
        self.read_mapper_stats.no_of_mapped_reads = {
            "zzz" : {"genome1" : 3, "genome2" : 4}, 
            "yyy" : {"genome1" : 1, "genome2" : 10}, 
            "xxx" : {"genome1" : 8,"genome2" : 10}}
        self.read_mapper_stats.no_of_mappings = {
            "zzz" : {"genome1" : 1, "genome2" : 12}, 
            "yyy" : {"genome1" : 2, "genome2" : 9}, 
            "xxx" : {"genome1" : 4,"genome2" : 18}}
        stat_fh = StringIO()
        self.read_mapper_stats._write_stats_to_fh(read_file_name, stat_fh)
        self.assertEqual(
            self.example_data.stat_file_content,
            stat_fh.getvalue())
        
    def test_head_line(self):
        read_file_name = ["xxx", "yyy", "zzz"]
        self.assertEqual(
            "\txxx\tyyy\tzzz", 
            self.read_mapper_stats._head_line(read_file_name))

    def test_value_line(self):
        read_file_name = ["xxx", "yyy", "zzz"]
        value_dict = {"zzz" : 5, "yyy" : 1, "xxx" : 8}
        self.assertEqual(
            "boing\t8\t1\t5",
            self.read_mapper_stats._value_line(
                "boing", value_dict, read_file_name))

    def test_dict_value_sum_line(self):
        read_file_name = ["xxx", "yyy", "zzz"]
        value_dict = {"zzz" : {"a" : 2, "b" : 4},
                      "yyy" : {"a" : 3, "b" : 1}, 
                      "xxx" : {"a" : 5, "b" : 9}}
        self.assertEqual(
            "boing\t14\t4\t6",
            self.read_mapper_stats._dict_value_sum_line(
                "boing", value_dict, read_file_name))

    def test_count_mappings(self):
        sam_fh_1 = StringIO(self.example_data.sam_content_1)
        no_of_mappings, no_of_mapped_reads = (
            self.read_mapper_stats._count_mappings(sam_fh_1))
        self.assertDictEqual(
            {'SL1344': 6, 'SL1344_plasmid1': 3, 'SL1344_plasmid2': 0},
            no_of_mappings)
        self.assertDictEqual(
            {'SL1344': 3.0, 'SL1344_plasmid1': 3.0, 'SL1344_plasmid2': 0},
            no_of_mapped_reads)

def _mock_count_fasta_entries(file_path):
    return("the number of entries in %s" % file_path) 
    

class ExampleData(object):

    fasta_seqs_1 = """>test_1 a random sequence
TTTAGAAATTACACA
>test_2 another random sequence
ACGAGAAATTAAATTAAATT
>test_3 another random sequence
TAGAGACATTGGATTTTATT
"""

    fasta_seqs_2 = """>test_1 a random sequence
TTTAGAAATTACACA
>test_2 another random sequence
ACGAGAAATTAAATTAAATT
>test_3 another random sequence
TAGAGACATTGGATTTTATT
>test_4 another random sequence
TAGAGACATTGGATTTTATT
>test_5 another random sequence
TAGAGACATTGGATTTTATT
"""

    sam_content_1 = """@HD	VN:1.0
@SQ	SN:SL1344	LN:10000000
@SQ	SN:SL1344_plasmid1	LN:5000
@SQ	SN:SL1344_plasmid2	LN:4000
@PG	ID:segemehl	VN:0.9.4-$Rev: 316 $ ($Date: 2011-08-18 16:37:19 +0200 (Thu, 18 Aug 2011) $)
read_01	0	SL1344	10	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:2	XA:Z:Q
read_01	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:2	XA:Z:Q
read_02	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_02	0	SL1344	1	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_02	0	SL1344	1000	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:3	XA:Z:Q
read_03	16	SL1344	1500	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XA:Z:Q
read_04	16	SL1344_plasmid1	1500	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XA:Z:Q
read_05	16	SL1344_plasmid1	1500	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XA:Z:Q
read_06	16	SL1344_plasmid1	1500	255	60M	*	0	0	ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT	*	NM:i:0	MD:Z:60	NH:i:1	XA:Z:Q
"""

    stat_file_content = """\txxx\tyyy\tzzz
Number of raw reads\t8\t1\t5
Reads long enough after clipping\t5\t23\t42
Reads too short after clipping\t2\t1\t0
Total number of mapped reads\t18\t11\t7
Total number of mappings\t22\t11\t13
"""

if __name__ == "__main__":
    unittest.main()
