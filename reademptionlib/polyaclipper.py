class PolyAClipper(object):
    """Uses very simple heuristics to detect and remove polyA tails."""

    def clip_poly_a_strech(self, sequence):
        """Search for a longer block of As and clip before it.

        Search for a block of four As. If the seven nucleotides after
        that are at least six A clip before the block of four As.
        """
        sequence = sequence.upper()
        length = 11
        if "AAAA" in sequence:
            for subseq, start_pos in self._aaaa_starting_substrings(
                sequence, length
            ):
                # Tolerate one mismatch
                if subseq.count("A") < length - 1:
                    continue
                else:
                    # Use sequence only to the start of the poly-A
                    # tail
                    sequence = sequence[:start_pos]
                    break
        return sequence

    def _aaaa_starting_substrings(self, sequence, length):
        """Find all AAAA starting substring of the given length"""
        start_pos = 0
        while start_pos != -1:
            start_pos = sequence.find("AAAA", start_pos)
            if start_pos == -1:
                break
            if start_pos + length > len(sequence):
                start_pos = -1
            else:
                cur_start_pos = start_pos
                start_pos = start_pos + 1
                yield (
                    [
                        sequence[cur_start_pos : cur_start_pos + length],
                        cur_start_pos,
                    ]
                )

    def remove_3_prime_a(self, sequence):
        """Remove 3' terminal As"""
        if sequence == "":
            return sequence
        elif sequence[-1] == "A":
            sequence = self.remove_3_prime_a(sequence[:-1])
        return sequence
