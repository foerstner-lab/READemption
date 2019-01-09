import concurrent.futures
import os
import sys
from reademptionlib.paths import Paths


class Helpers(object):

    def __init__(self, args):
        self._paths = Paths(args)
        self._args = args

    def test_folder_existance(self, task_specific_folders):
        """Test the existance of required folders."""
        for folder in (
                self._paths.required_base_folders() + task_specific_folders):
            if not os.path.exists(folder):
                self.write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder)
        
    def file_needs_to_be_created(self, file_path, quiet=False):
        """Test if a file exists of need to be created."""
        if not self._args.check_for_existing_files:
            return True
        if os.path.exists(file_path):
            if not quiet:
                sys.stderr.write(
                    "File %s exists. Skipping its generation.\n" % file_path)
            return False
        return True

    def check_job_completeness(self, jobs):
        """Check the completness of each job in a list"""
        for job in concurrent.futures.as_completed(jobs):
            if job.exception():
                raise(job.exception())

    def write_err_msg_and_quit(self, msg):
        """Write error message and close the program gracefully."""
        sys.stderr.write(msg)
        sys.exit(1)

    def was_paired_end_alignment(self, lib_names):
        """Check if the mapping was done in paired- or single-end mode"""
        if len(lib_names) * 2 == len(self._paths.get_read_files()):
            return True
        return False
        

