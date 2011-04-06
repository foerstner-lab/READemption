import os
import sys
import json

from rapl.raplpathes import RaplPathes

class RaplCreator(object):

    def __init__(self):
        self.pathes = RaplPathes()

    def create_root_folder(self, project_name):
        """Create the root folder of a new project with the given name.
        
        Arguments:
        - `project_name`: Name of the project root folder

        """
        if not os.path.exists(project_name):
            os.mkdir(project_name)
        else:
            sys.stderr.write("Cannot create folder \"%s\"! File/folder with "
                             "the same name exists already.\n" % project_name)
            sys.exit(2)

    def create_subfolders(self, project_name):
        """Create required subfolders in the given folder.
        
        Arguments:
        - `project_name`: Name of the project root folder

        """
        for folder in self.pathes.required_folders():
            folder_in_root_folder = "%s/%s" % (project_name, folder)
            if not os.path.exists(folder_in_root_folder):
                os.mkdir(folder_in_root_folder)
    
    def create_config_file(self, project_name):
        """Create a JSON config file.
        
        Arguments:
        - `project_name`: Name of the project root folder
        """
        config_fh = open("%s/%s" % (project_name, self.pathes.config_file), "w")
        config_fh.write(json.dumps({"annotation_and_genomes_files" : {}}))
        config_fh.close()
