import unittest
import sys
from io import StringIO
sys.path.append(".")
from libs.readclipper import ReadClipper

class TestReadClipper(unittest.TestCase):

    def setUp(self):
        self.read_clipper = ReadClipper()
        self.example_data = ExampleData()

    def tearDown(self):
        pass

    def test_clip_with_wrong_file_list_sizes(self):
        """If the number of input file is not equal the number of
        output file an exception must be raised."""
        self.assertRaises(
            Exception, self.read_clipper.clip, ["a"], ["a", "b"])

    def test_clip_entries_in_file(self):
        input_fh = StringIO(self.example_data.fasta_input)
        output_fh = StringIO()
        self.read_clipper._clip_entries_in_file(input_fh, output_fh)
        self.assertEqual(output_fh.getvalue(), self.example_data.fasta_output)

class ExampleData(object):

    fasta_input = """>test_1 a seq without poly a tale
CTGAAATATACAGGGACACAAAAACCGGGTAGAGTACG
>test_2 a seq with an internal poly A tail
CTGAAATATACAGGGACACAAAAAAAAAAAAAAAATTTTATTACCCGCACAT
>test_3 a seq with a long terminal poly A strech
AACCAGACCACTAAACCGGATTTAGTAAAAAAAAAAAAAAAAAAAA
>test_3 a seq with a short terminal poly A strech
TTAAATTGGTATGATTTAGTAAAAAA
>test_4 a seq with a long terminal poly A strech and line break
GTTTGTAGTATGGAGATAGGGAAAGAGTT
CCAGTAAAAAAAAAAAAAAAAAAAAAAAA
"""
    fasta_output = """>test_1 a seq without poly a tale
CTGAAATATACAGGGACACAAAAACCGGGTAGAGTACG
>test_2 a seq with an internal poly A tail
CTGAAATATACAGGGACAC
>test_3 a seq with a long terminal poly A strech
AACCAGACCACTAAACCGGATTTAGT
>test_3 a seq with a short terminal poly A strech
TTAAATTGGTATGATTTAGT
>test_4 a seq with a long terminal poly A strech and line break
GTTTGTAGTATGGAGATAGGGAAAGAGTTCCAGT
"""
        
if __name__ == "__main__":
    unittest.main()
