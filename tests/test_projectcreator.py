import unittest
import os
import sys
import shutil
import subprocess
import Bio
import matplotlib
import pandas as pd
import pysam

sys.path.append(".")
from reademptionlib.projectcreator import ProjectCreator


class TestProjectCreator(unittest.TestCase):
    def setUp(self):
        self.root_folder_name = "a_test_project"
        self.config_file_path = "config.json"
        self.version_file_path = ".version_log.txt"
        self.projectcreator = ProjectCreator()

    def tearDown(self):
        if os.path.exists(self.root_folder_name):
            shutil.rmtree(self.root_folder_name)
        if os.path.exists(self.config_file_path):
            os.remove(self.config_file_path)
        if os.path.exists(self.version_file_path):
            os.remove(self.version_file_path)

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
        config_json_content_expected = (
            '{"species": {"human": "Homo sapiens", '
            '"staph": "Staphylococcus aureus"}}'
        )
        self.projectcreator.create_config_file(
            self.config_file_path,
            {"human": "Homo sapiens", "staph": "Staphylococcus aureus"},
        )
        with open(self.config_file_path, "r") as config_file:
            read_in_config_file_content = config_file.read()
            self.assertEqual(
                read_in_config_file_content, config_json_content_expected
            )

    def test_create_version_file(self):
        # Get READemption version
        reademption_output = subprocess.run(
            ["reademption", "--version"], stdout=subprocess.PIPE
        ).stdout.decode("utf-8")
        reademption_version_not_cleaned = reademption_output.split(" ")[-1]
        reademption_version = reademption_version_not_cleaned.replace("\n", " ")
        # Get Python version
        python_version = sys.version.replace("\n", " ")

        # Create the version file
        self.projectcreator.create_version_file(
            self.version_file_path, reademption_version
        )
        # The expected content of the version file. The version number of the
        # modules are retrieved from the actual system where testing happens
        expected_version_file_content = (
            f"READemption version: {reademption_version}\n"
            f"Python version: {python_version}\n"
            f"Biopython version: {Bio.__version__}\n"
            f"pysam version: {pysam.__version__}\n"
            f"matplotlib version: {matplotlib.__version__}\n"
            f"pandas version: {pd.__version__}\n"
        )
        # Compare the version file with the expected content
        with open(self.version_file_path, "r") as version_file:
            version_file_content = version_file.read()
            self.assertMultiLineEqual(
                version_file_content, expected_version_file_content
            )


if __name__ == "__main__":
    unittest.main()
