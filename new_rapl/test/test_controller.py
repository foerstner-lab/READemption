import os
import sys
import unittest
import shutil
sys.path.append(".")
from libs.controller import Controller

class ArgMock(object):
    project_name = None

class TestController(unittest.TestCase):

    def setUp(self):
        self.controller = Controller()
        self.test_project_name = "a_test_project"

    def tearDown(self):
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)
    
    def test_start_project(self):
        arg_mock = ArgMock()
        arg_mock.project_name = self.test_project_name
        self.controller.start_project(arg_mock)
        self.assertEqual(
            list(os.listdir(self.test_project_name)), 
            ['rapl.config', 'input', 'output'])
        if os.path.exists(self.test_project_name):
            shutil.rmtree(self.test_project_name)

if __name__ == "__main__":
    unittest.main()

