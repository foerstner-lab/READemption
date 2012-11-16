import sys
import unittest
sys.path.append(".")
from libs.parameters import Parameters

class TestParameters(unittest.TestCase):

    def setUp(self):
        self.parameters = Parameters()

    def test_int_parameters(self):
        for parameter in [
            self.parameters.segemehl_hit_strategy,
            self.parameters.python_number_of_threads,
            self.parameters.min_seq_length,
            self.parameters.min_overlap,
            self.parameters.uniquely_mapped_reads_only]:
            assert(isinstance(parameter, int))

    def test_float_parameters(self):
        for parameter in [
            self.parameters.min_read_overlap_perc]:
            assert(isinstance(parameter, float))

    def test_boolean_parameters(self):
        for parameter in [
            self.parameters.uniquely_mapped_reads_only]:
            assert(isinstance(parameter, bool))

    def test_string_parameters(self):
        for parameter in [
            self.parameters.exception_handling]:
            assert(isinstance(parameter, str))

if __name__ == "__main__":
    unittest.main()
