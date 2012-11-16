import sys
from subprocess import call
sys.path.append(".")

class AnnotationOverlap(object):

    def __init__(self, intersectbed_bin="intersectBed"):
        self.intersectbed_bin = intersectbed_bin
    
    def find_overlaps(self, read_mapping_path, annotation_file_path, 
                      output_file_path):
        output_fh = open(output_file_path, "w")
        call([self.intersectbed_bin, "-bed", "-wo", "-abam", read_mapping_path,
              "-b", annotation_file_path], stdout=output_fh)
