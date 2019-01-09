import sys
from reademptionlib.paths import Paths
from reademptionlib.projectcreator import ProjectCreator


class ProjectCreateController(object):
    """Create a READemption project including the necessary
folders in order to perform ongoing procession steps."""

    def __init__(self, args):
        """Create an instance."""
        self._args = args
        self._paths = Paths(args)

    def create_project(self, version):
        """Create a new project."""
        sys.stdout.write(
            "   ___  _______   ___                 __  _\n"
            "  / _ \/ __/ _ | / _ \___ __ _  ___  / /_(_)__  ___\n"
            " / , _/ _// __ |/ // / -_)  ' \/ _ \/ __/ / _ \/ _ \\\n"
            "/_/|_/___/_/ |_/____/\__/_/_/_/ .__/\__/_/\___/_//_/\n"
            "                             / /\n"
            "====================================================\n"
            "========================================\n"
            "=======================\n"
            "==============\n\n"
            "[http://pythonhosted.org/READemption/]\n\n")
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders())
        project_creator.create_version_file(self._paths.version_path, version)
        sys.stdout.write("Created folder \"%s\" and required subfolders.\n" % (
            self._args.project_path))
        print(self._paths.base_path)
        sys.stdout.write(
            "Please copy read files into folder \"%s\" and "
            "reference sequences files into folder \"%s\".\n" % (
                self._paths.read_fasta_folder, self._paths.ref_seq_folder))
