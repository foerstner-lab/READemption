import unittest
import os
import sys
import shutil

sys.path.append(".")
from reademptionlib.projectcreator import ProjectCreator


class TestProjectCreator(unittest.TestCase):
    def setUp(self):
        self.root_folder_name = "a_test_project"
        self.species_file_path = ".species.json"
        self.projectcreator = ProjectCreator()

    def tearDown(self):
        if os.path.exists(self.root_folder_name):
            shutil.rmtree(self.root_folder_name)
        if os.path.exists(self.species_file_path):
            os.remove(self.species_file_path)

    def test_create_root_folder(self):
        self.projectcreator.create_root_folder(self.root_folder_name)
        assert os.path.exists(self.root_folder_name)
        shutil.rmtree(self.root_folder_name)

    def test_create_subfolders(self):
        self.projectcreator.create_root_folder(self.root_folder_name)
        subfolders = ["test_a", "test_b", "test_c"]
        subfolders = [
            self.root_folder_name + "/" + subfolder for subfolder in subfolders
        ]
        self.projectcreator.create_subfolders(subfolders)
        for subfolder in subfolders:
            assert os.path.exists(subfolder)

    def test_create_species_file(self):
        species_json_content = (
            "{\"human\": \"Homo sapiens\", "
            "\"staph\": \"Staphylococcus aureus\"}")
        self.projectcreator.create_species_file(self.species_file_path,
            ['human=Homo sapiens',
             'staph=Staphylococcus aureus'])
        with open(self.species_file_path ,"r") as species_file:
            self.assertEqual(species_file.read(), species_json_content)

if __name__ == "__main__":
    unittest.main()
