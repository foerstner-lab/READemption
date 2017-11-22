import os
import shutil
import sys
sys.path.append("./tests")
from reademptionlib.projectcreator import ProjectCreator


def test_project_creator():
    # Import Object & define dummy project name
    root_folder_name = "a_test_project"
    projectcreator = ProjectCreator()

    # Test root folder creation
    projectcreator.create_root_folder(root_folder_name)
    assert os.path.exists(root_folder_name)
    shutil.rmtree(root_folder_name)

    # Test subfolder creation
    projectcreator.create_root_folder(root_folder_name)
    subfolders = ["test_a", "test_b", "test_c"]
    subfolders = [root_folder_name + "/" + subfolder for
                  subfolder in subfolders]
    projectcreator.create_subfolders(subfolders)
    for subfolder in subfolders:
        assert os.path.exists(subfolder)

    # Check if dummy project exists & delete it
    if os.path.exists(root_folder_name):
        shutil.rmtree(root_folder_name)
    
