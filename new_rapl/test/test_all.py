import unittest
from test_controller import TestController
from test_fasta import TestFasta
from test_parameters import TestParameters
from test_paths import TestPaths
from test_polyaclipper import TestPolyAClipper
from test_projectcreator import TestProjectCreator
from test_readclipper import TestReadClipper
from test_readmapper import TestReadMapper
from test_segemehl import TestSegemehlIndexBuilding, TestSegemehlMapping
from test_seqsizefilter import TestSeqSizeFilter

test_cases = [
    TestController,
    TestFasta,
    TestParameters,
    TestPaths,
    TestPolyAClipper,
    TestProjectCreator,
    TestReadClipper,
    TestReadMapper,
    TestSegemehlIndexBuilding,
    TestSegemehlMapping,
    TestSeqSizeFilter]

suites = [unittest.TestLoader().loadTestsFromTestCase(test_case) 
          for test_case in test_cases]
combined_test = unittest.TestSuite(suites)
unittest.TextTestRunner(verbosity=2).run(combined_test)

