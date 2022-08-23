import os
import pysam


class FragmentBuilder(object):
    def __init__(
        self,
        max_fragment_length
    ):
        self._max_fragment_length = max_fragment_length
        self.sam_flag_property_to_value = {
            "read paired": 1,
            "read mapped in proper pair": 2,
            "read unmapped": 4,
            "mate unmapped": 8,
            "read reverse strand": 10,
            "mate reverse strand": 20,
            "first in pair": 40,
            "second in pair": 80,
            "not primary alignment": 100,
            "read fails platform/vendor quality checks": 200,
            "read is PCR or optical duplicate": 400,
            "supplementary alignment": 800,
            # fragment = pysam.AlignedSegment()
            # 1 QNAME
            # fragment.query_name = alignment.query_name
            # 2 FLAG
            # fragment.flag = alignment.flag
            # 3 RNAME
            # fragment.reference_id = alignment.reference_id
            # 4 POS
            # fragment.reference_start = alignment.reference_start
            # 5 MAPQ
            # fragment.mapping_quality = alignment.mapping_quality
            # 6 CIGAR
            # fragment.cigarstring = alignment.cigarstring
            # 7 RNEXT
            # fragment.next_reference_id = alignment.next_reference_id
            # 8 PNEXT
            # fragment.next_reference_start = alignment.next_reference_start
            # 9 TLEN
            # fragment.template_length = alignment.template_length
            # 10 SEQ
            # fragment.query_sequence = alignment.query_sequence
            # 11 QUAL
            # fragment.query_qualities = alignment.query_qualities
            # 12 TAGS
            # fragment.tags = alignment.tags
        }

    def build_bam_file_with_fragments(
        self, read_alignment_path, fragment_alignment_path
    ):
        with pysam.Samfile(read_alignment_path) as input_bam:
            with pysam.Samfile(
                fragment_alignment_path, "wb", header=input_bam.header
            ) as output_bam:
                # collect the two reads of a pair for one alignment.
                # This is only possible because the bam file has been ordered
                # by read name and alignment hit index
                pair = {}
                for alignment in input_bam.fetch(until_eof=True):
                    # Only build a fragment if the pair is proper paired and the
                    # current read is the second in pair, to avoid building the same
                    # fragement twice (one time for read one and a second time for
                    # read two)

                    if alignment.is_proper_pair and alignment.is_read2 is False:
                        # Don't build a fragment or add the read if the read is
                        # read one and proper paired. The read is cached in the
                        # pair dictionary
                        pair["read1"] = alignment

                    elif alignment.is_proper_pair is False:
                        # if the read is not mapped in proper pair, add the
                        # single read
                        fragment_length = abs(alignment.template_length)
                        if self._max_fragment_length:
                            if fragment_length > self._max_fragment_length:
                                continue
                        output_bam.write(alignment)

                    else:
                        # Add read2 to the cached pair dictionary and build the
                        # fragment
                        pair["read2"] = alignment
                        orientation = self._determine_orientation(pair["read1"])
                        if pair["read1"].is_read2:
                            print(f"read 1: {pair['read1']} is actually read 2: {pair['read2']}")
                        if pair["read2"].is_read1:
                            print(f"read 2: {pair['read2']} is actually read 1: {pair['read1']}")
                        # build fragment for read one and read two
                        start, end = self._build_fragment(
                            pair["read1"], pair["read2"], orientation, read_alignment_path
                        )
                        # no need of adding + 1 to the fragment length, since
                        # start is 0-based and end is 1-based
                        fragment_length = end - start
                        pair["read1"].reference_start = start
                        pair["read1"].query_sequence = fragment_length * "N"
                        pair["read1"].cigarstring = f"{fragment_length}="
                        if self._max_fragment_length:
                            if fragment_length > self._max_fragment_length:
                                continue
                        output_bam.write(pair["read1"])
        # Sort and index the resulting alignment file, since building new
        # fragments results in the new file being out of order
        tmp_sorted_outfile = f"{fragment_alignment_path}_sorted"
        pysam.sort("-o", tmp_sorted_outfile, fragment_alignment_path)
        os.rename(tmp_sorted_outfile, fragment_alignment_path)
        pysam.index(fragment_alignment_path)

    def _build_fragment(self, read1, read2, orientation, read_alignment_path):
        tlen = abs(read1.template_length)
        # start is 0-based
        start = read1.reference_start
        # end is 1-based
        end = read1.reference_end
        start_of_mate = read1.next_reference_start
        # check if the reads are in order
        if orientation == "forward":
            if start <= start_of_mate:
                # Forward read pair in order
                start_of_fragment = start
                end_of_fragment = start + tlen
            else:
                # Forward read pair not in order
                mate = read2
                mate_start = mate.pos
                mate_end = mate.aend
                # check if overlap or exceed. One needs to be added to the
                # start because the start is 0-based and the end is 1-based
                if mate_end < (start + 1):
                    # exceed
                    start_of_fragment = mate_start
                    end_of_fragment = end
                else:
                    start_of_fragment = start
                    end_of_fragment = mate_end

        elif orientation == "reverse":
            if start_of_mate <= start:
                # Reverse read pair in order
                start_of_fragment = start_of_mate
                end_of_fragment = start_of_mate + tlen

            else:
                # Reverse read pair not in order
                mate = read2
                mate_start = mate.pos
                mate_end = mate.aend
                # check if overlap or exceed. One needs to be added to the
                # start because the start is 0-based and the end is 1-based
                if end < (mate_start + 1):
                    # exceed
                    start_of_fragment = start
                    end_of_fragment = mate_end
                else:
                    # overlap
                    start_of_fragment = mate_start
                    end_of_fragment = end
        return start_of_fragment, end_of_fragment

    def _determine_orientation(self, alignment):
        # Determine the orientation of a read:
        # For paired end reads the read one of a pair determines the
        # orientation. If read one has the flag forward it maps to the
        # forward strand and if it has the flag reverse it maps to the
        # reverse strand. If read two has the flag forward it maps to the
        # reverse strand and if it has the flag reverse it maps to the
        # forward strand.
        if (alignment.is_reverse is False and alignment.is_read2 is False) or (
            alignment.is_reverse is True and alignment.is_read2 is True
        ):
            orientation = "forward"
        else:
            orientation = "reverse"
        return orientation
