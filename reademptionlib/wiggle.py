class WiggleWriter(object):
    def __init__(self, track_str, fh):
        self._fh = fh
        self._fh.write(('track type=wiggle_0 name="%s"\n' % (track_str)))

    def write_replicons_coverages(
        self, replicon_str, coverages, discard_zeros=True, factor=1.0
    ):
        self._fh.write("variableStep chrom=%s span=1\n" % (replicon_str))
        # Filter values of 0 and multiply the remaining ones by
        # the given factor. pos is increased by 1 as a translation
        # from a 0-based sysem (Python list) to a 1 based system
        # (wiggle) takes place.
        self._fh.write(
            "\n".join(
                [
                    "%s %s" % (pos + 1, coverage * factor)
                    for pos, coverage in filter(
                        lambda pos_and_cov: pos_and_cov[1] != 0.0,
                        enumerate(coverages),
                    )
                ]
            )
            + "\n"
        )

    def close_file(self):
        self._fh.close()
