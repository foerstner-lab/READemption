import hashlib
import os
import sys
sys.path.append("./tests")
import segemehl_data as sd


def create_tmp_fasta_file(fasta_file_path, content):
    fasta_fh = open(fasta_file_path, "w")
    fasta_fh.write(content)
    fasta_fh.close()


def sha1_of_file(file_path):
    fh = open(file_path, "rb")
    content = fh.read()
    fh.close()
    return hashlib.sha1(content).hexdigest()


def setup_function(function):
    sd.data_segemehl()
    create_tmp_fasta_file(sd.fasta_file_path, sd.genome_fasta_upper)
    sd.segemehl.build_index([sd.fasta_file_path], sd.index_file_path)


def teardown_function(function):
    for file_path in [sd.fasta_file_path, sd.index_file_path]:
        if os.path.exists(file_path):
            os.remove(file_path)


# Test segemehl index building
def test_build_index_lower_letters():
    create_tmp_fasta_file(sd.fasta_file_path, sd.genome_fasta_lower)
    sd.segemehl.build_index([sd.fasta_file_path], sd.index_file_path)
    assert sha1_of_file(
        sd.index_file_path) == "78668505720e53735f807bb5485b0b38cc3cbc22"

    
def test_build_index_upper_letters():
    # create_tmp_fasta_file(sd.fasta_file_path, sd.genome_fasta_upper)
    sd.segemehl.build_index([sd.fasta_file_path], sd.index_file_path)
    assert sha1_of_file(
        sd.index_file_path) == "78668505720e53735f807bb5485b0b38cc3cbc22"


# Test segemehl aligning
def test_align_reads_single_read_perfect_match():
    """
    ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
    """
    read_file_content = (
        ">read_01\n" +
        "ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT"
        + "\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_aligned_1


def test_map_reads_single_read_not_matching():
    """
    ACAACATCCATGAACCGCATCAGCACCACCACCATTACCACCATCACCATTACCACAGGT
    |       | |||     ||    |            |   ||  |            ||
    ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT
    """
    read_file_content = (
        ">read_02\n" +
        "ATGTACCACATGAGAGAGATAGAGAGAGATTGACAACCACACACGAGAGAGAGAGAGAGT"
        + "\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_no_aligned


def test_map_reads_single_read_one_mismatch():
    """A 20 nt long read with 1 mismatch at 95% accu should be
    mapped.
    GCTTTTTTTTCGACCAGAGA
    |||||||||||||||||| |
    GCTTTTTTTTCGACCAGACA
    """
    read_file_content = (
        ">read_03\nGCTTTTTTTTCGACCAGACA\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_aligned_2


def test_map_reads_single_read_two_mismatches_95():
    """A 20 nt long read with 2 mismatches at 95% accuracy should
    not be mapped.
    GCTTTTTTTTCGACCAGAGA
    |||||||||||||||||  |
    GCTTTTTTTTCGACCAGTCA
    """
    read_file_content = (
        ">read_04\nGCTTTTTTTTCGACCAGTCA\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_no_aligned


def test_map_reads_single_read_two_mismatches_90():
    """A 20 nt long read with 2 mismatches at 90% accuracy should
    be mapped.
    GCTTTTTTTTCGACCAGAGA
    |||||||||||||||||  |
    GCTTTTTTTTCGACCAGTCA
    """
    read_file_content = (
        ">read_05\nGCTTTTTTTTCGACCAGTCA\n")
    assert align_reads_and_give_result(
        read_file_content, accuracy=90) == sd.sam_result_aligned_3


def test_map_reads_single_read_three_mismatches():
    """A 20 nt long read with 3 mismatches at 90% accuracy should
    not be mapped.
    GCTTTTTTTTCGACCAGAGA
    ||||| |||||||||||  |
    GCTTTATTTTCGACCAGTCA
    """
    read_file_content = (
        ">read_06\nGCTTTTTTTTCGACCAGTCA\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_no_aligned


def test_map_reads_single_too_short_read():
    """Reads that are too short should be mapped
    """
    read_file_content = (
        ">read_07\nGCTTTTTTT\n")
    assert align_reads_and_give_result(
        read_file_content) == sd.sam_result_no_aligned


def align_reads_and_give_result(read_file_content, **kwargs):
    """
    - read_file_content: the content of a read file (in fasta format)
    - **kwargs: are directly given to map_reads()
    """
    create_tmp_fasta_file(
        sd.read_fasta_file_path, read_file_content)
    sd.segemehl.run_alignment(sd.read_fasta_file_path, sd.index_file_path,
                              [sd.fasta_file_path], sd.aligning_result_path,
                              sd.unmapped_reads_path, **kwargs)
    result_fh = open(sd.aligning_result_path)
    result = result_fh.read()
    result_fh.close()
    if '--silent' in result:
        map_spec, reads = result.split("--silent", 1)
    return reads
    
