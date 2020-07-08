import os
import sys
import pysam


class BamMerger(object):
    def merge(self, output_bam, input_bam_1, input_bam_2):
        # Check input files
        for input_bam in [input_bam_1, input_bam_2]:
            if not os.path.exists(input_bam) is True:
                sys.stderr.write(
                    "Input file %s does exist. Merging not possible.\n"
                    % input_bam
                )
                return
        # Check output file
        if os.path.exists(output_bam) is True:
            sys.stderr.write(
                "Output file %s already exists. Merging not possible.\n"
                % output_bam
            )
            return
        # Merge and generate index
        pysam.merge(output_bam, input_bam_1, input_bam_2)
        pysam.index(output_bam)
