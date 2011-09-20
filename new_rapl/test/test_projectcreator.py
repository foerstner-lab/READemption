import unittest
import os
import sys
import shutil
sys.path.append(".")
from libs.projectcreator import ProjectCreator

class TestProjectCreator(unittest.TestCase):

    def setUp(self):
        self.projectcreator = ProjectCreator()
        self.root_folder_name = "a_test_project"
        
    def tearDown(self):
        if os.path.exists(self.root_folder_name):
            shutil.rmtree(self.root_folder_name)

    def test_create_root_folder(self):
        self.projectcreator.create_root_folder(self.root_folder_name)
        assert(os.path.exists(self.root_folder_name))
        shutil.rmtree(self.root_folder_name)
    
    def test_create_subfolders(self):
        self.projectcreator.create_root_folder(self.root_folder_name)
        subfolders = ["test_a", "test_b", "test_c"]
        self.projectcreator.create_subfolders(self.root_folder_name, subfolders)
        for subfolder in subfolders:
            assert(os.path.exists(self.root_folder_name + "/" + subfolder))

    def test_create_config_file(self):
        file_name = "test_rapl_file.json"
        self.projectcreator.create_root_folder(self.root_folder_name)
        self.projectcreator.create_config_file(
            self.root_folder_name, file_name)
        fh = open("%s/%s" % (self.root_folder_name, file_name))
        content = fh.read()
        fh.close()
        self.assertEqual(content, '{"annotation_and_genomes_files": {}}')

if __name__ == "__main__":
    unittest.main()
