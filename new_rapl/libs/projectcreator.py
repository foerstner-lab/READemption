import json
import os
import sys

class ProjectCreator(object):

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

    def create_subfolders(self, subfolders):
        """Create required subfolders in the given folder.
        
        Arguments:
        - `project_name`: Name of the project root folder

        """
        for folder in subfolders:
            if not os.path.exists(folder):
                os.mkdir(folder)
    
    def create_config_file(self, config_file_path):
        """Create a JSON config file.
        
        Arguments:
        - `project_name`: Name of the project root folder
        """
        config_fh = open(config_file_path, "w")
        config_fh.write(json.dumps({"annotation_and_genomes_files" : {}}))
        config_fh.close()
