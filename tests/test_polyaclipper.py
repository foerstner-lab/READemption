import sys
sys.path.append("./tests")
from reademptionlib.polyaclipper import PolyAClipper


def test_poly_a_clipper():
    # Define Object
    poly_a_clipper = PolyAClipper()

    # Test: If there is no poly a stretch the sequence is not changed
    seq_no_stretch = result_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAC"
    "AATGGTTAGGTACAGATAG"
    assert poly_a_clipper.clip_poly_a_stretch(seq_no_stretch) == result_seq

    # Test: If the sequence is empty the sequence is not changed
    empty_seq = result_seq = ""
    assert poly_a_clipper.clip_poly_a_stretch(empty_seq) == result_seq

    # Test: Clip a terminal 10 fold A stretch
    test_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAAATTTAGACGACG"
    result_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACAC"
    assert poly_a_clipper.clip_poly_a_stretch(test_seq) == result_seq

    # Test: If there is less than a 10 fold terminal stretch don't clip
    test_seq = result_seq = "ATAGTAGGAGATTTAGACCAGATGACGATGACACAAAAAAAAA"
    "TTTAGACGACG"
    assert poly_a_clipper.clip_poly_a_stretch(test_seq) == result_seq

    # Test: All A sequence
    test_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    result_seq = ""
    assert poly_a_clipper.clip_poly_a_stretch(test_seq) == result_seq

    # Test: Sequence with an starting 'AAAA' substring
    test_seq = "TTTAAAATTTTTTTTAAAACCCCCCCCCCAAAAC"
    assert list(poly_a_clipper._aaaa_starting_substrings(
        test_seq, 11)) == [['AAAATTTTTTT', 3], ['AAAACCCCCCC', 15]]

    # Test: If there is no terminal A stretch, there is no clipping
    test_seq = result_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
    assert poly_a_clipper.remove_3_prime_a(test_seq) == result_seq

    # Test: Removal of one terminal A
    test_seq = "AAAAATTTTCCGCCCGGGAAATTTTA"
    result_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
    assert poly_a_clipper.remove_3_prime_a(test_seq) == result_seq

    # Test: Removal of multiple terminal As
    test_seq = "AAAAATTTTCCGCCCGGGAAATTTTAAAAAA"
    result_seq = "AAAAATTTTCCGCCCGGGAAATTTT"
    assert poly_a_clipper.remove_3_prime_a(test_seq) == result_seq
