import pysam


class FragmentBuilder(object):
    def __init__(
        self,
    ):
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
        }

    def build_sam_file_with_fragments(self, read_alignment_path, fragment_alignment_path):
        with pysam.Samfile(read_alignment_path) as input_bam:
            with pysam.Samfile(
                    fragment_alignment_path, "wb", header=input_bam.header
            ) as output_bam:
                for alignment in input_bam.fetch():
                    orientation = self.determine_orientation(alignment)
                    # Only build a fragment if the pair is proper paired and the
                    # current read is the first in pair, to avoid building the same
                    # fragement twice (one time for read one and a second time for
                    # read two)

                    #fragment = pysam.AlignedSegment()
                    #1 QNAME
                    #fragment.query_name = alignment.query_name
                    #2 FLAG
                    #fragment.flag = alignment.flag
                    #3 RNAME
                    #fragment.reference_id = alignment.reference_id
                    #4 POS
                    #fragment.reference_start = alignment.reference_start
                    #5 MAPQ
                    #fragment.mapping_quality = alignment.mapping_quality
                    #6 CIGAR
                    #fragment.cigarstring = alignment.cigarstring
                    #7 RNEXT
                    #fragment.next_reference_id = alignment.next_reference_id
                    #8 PNEXT
                    #fragment.next_reference_start = alignment.next_reference_start
                    #9 TLEN
                    #fragment.template_length = alignment.template_length
                    #10 SEQ
                    #fragment.query_sequence = alignment.query_sequence
                    #11 QUAL
                    #fragment.query_qualities = alignment.query_qualities
                    #12 TAGS
                    #fragment.tags = alignment.tags


                    if (alignment.is_proper_pair and alignment.is_read2 is False):
                        print("build fragment")
                        start, end = self._build_fragment(alignment, orientation, read_alignment_path)
                        # no need of adding + 1 to the fragment length, since
                        # start is 0-based and end is 1-based
                        fragment_length = end - start
                        alignment.query_sequence = fragment_length * "N"
                        alignment.cigarstring = f"{fragment_length}="
                        output_bam.write(alignment)

                    elif alignment.is_proper_pair is False:
                        print("no proper pair: add single read")
                        output_bam.write(alignment)
                    # Don't build a fragment if the read is read two and proper
                    # paired
                    else:
                        print("don't build a fragment or read")
                        continue


                    #output_bam.write(alignment)
        #tmp_sorted_outfile = f"{fragment_alignment_path}_sorted"
        #pysam.sort("-o", tmp_sorted_outfile, fragment_alignment_path)
        #os.rename(tmp_sorted_outfile, fragment_alignment_path)
        pysam.index(fragment_alignment_path)


    def _build_fragment(self, entry, orientation, read_alignment_path):
        tlen = abs(entry.template_length)
        # start is 0-based
        start = entry.reference_start
        # end is 1-based
        end = entry.reference_end
        start_of_mate = entry.next_reference_start
        # check if the reads are in order
        print(f"start: {start}")
        print(f"end: {end}")
        print(f"start of mate from read one: {start_of_mate}")
        print(f"tlen: {tlen}")
        if orientation == "forward":
            if start <= start_of_mate:
                print("fragment is in order")
                # Forward read pair in order
                start_of_fragment = start
                end_of_fragment = start + tlen
            else:
                # Forward read pair not in order
                print("fragment not in order")
                with pysam.Samfile(read_alignment_path) as bam:
                    mate = bam.mate(entry)
                mate_start = mate.pos
                print(f"mate_start: {mate_start}")
                mate_end = mate.aend
                print(f"mate_end: {mate_end}")
                # check if overlap or exceed. One needs to be added to the
                # start because the start is 0-based and the end is 1-based
                if mate_end < (start + 1):
                    # exceed
                    print("exceed")
                    start_of_fragment = mate_start
                    end_of_fragment = end
                else:
                    print("overlap")
                    start_of_fragment = start
                    end_of_fragment = mate_end

        elif orientation == "reverse":
            if start_of_mate <= start:
                # Reverse read pair in order
                print("fragment is in order")
                start_of_fragment = start_of_mate
                end_of_fragment = start_of_mate + tlen

            else:
                # Reverse read pair not in order
                print("fragment is not in order")
                with pysam.Samfile(read_alignment_path) as bam:
                    mate = bam.mate(entry)
                mate_start = mate.pos
                print(f"mate_start: {mate_start}")
                mate_end = mate.aend
                print(f"mate_end: {mate_end}")
                # check if overlap or exceed. One needs to be added to the
                # start because the start is 0-based and the end is 1-based
                if end < (mate_start + 1):
                    print("exceed")
                    start_of_fragment = start
                    end_of_fragment = mate_end
                else:
                    print("overlap")
                    start_of_fragment = mate_start
                    end_of_fragment = end

        print(f"start of fragment: {start_of_fragment}")
        print(f"end of fragment: {end_of_fragment}")
        return start_of_fragment, end_of_fragment

    def determine_orientation(self, alignment):
        # Determine the orientation of a read:
        # For paired end reads the read one of a pair determines the
        # orientation. If read one has the flag forward it maps to the
        # forward strand and if it has the flag reverse it maps to the
        # reverse strand. If read two has the flag forward it maps to the
        # reverse strand and if it has the flag reverse it maps to the
        # forward strand.
        if (alignment.is_reverse is False and alignment.is_read2 is False) or (
                alignment.is_reverse is True and alignment.is_read2 is True):
            orientation = "forward"
        else:
            orientation = "reverse"
        return orientation
