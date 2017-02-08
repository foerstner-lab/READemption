import os
import shutil
import sys
sys.path.append("./tests")
import controller_data as cd


def setup_function(function):
    cd.data_controllers()
    cd.project_creator.create_project(cd.version)

    
def teardown_function(function):
    if os.path.exists(cd.test_project_name):
        shutil.rmtree(cd.test_project_name)
        

def test_create_project():
    assert set(list(os.listdir(cd.test_project_name))) == set([
        'input', 'output'])

    
def test_read_aligning():
    # generate input fasta files
    genome_fh = open("{}/{}".format(cd.controller_align._paths.ref_seq_folder,
                                    "agenome.fa"), "w")
    read_fh_1 = open("{}/{}".format(
        cd.controller_align._paths.read_fasta_folder, "libfoo.fa"), "w")
    read_fh_2 = open("{}/{}".format(
        cd.controller_align._paths.read_fasta_folder, "libbar.fa"), "w")
    
    genome_fh.write(cd.genome_fasta)
    genome_fh.close()
    # File containing overlap sequences with reference but not annotation
    read_fh_1.write(cd.read_fasta_1)
    read_fh_1.close()
    # file containing overlap sequences with the reference and annotation
    read_fh_2.write(cd.read_fasta_2)
    read_fh_2.close()
    
    # generate annotation files
    annotation_fh = open("{}/some_annos.gff".format(
        cd.controller_align._paths.annotation_folder), "w")
    annotation_fh.write(cd.gff_content_1)
    annotation_fh.close()

    # Perform the mapping
    cd.controller_align._paths._set_folder_names()
    cd.controller_align.align_reads()

    # Create coverage files
    cd.controller_coverage._paths.required_coverage_folders()
    cd.controller_coverage.create_coverage_files()

    # Perform gene wise quantification
    cd.controller_genequanti.quantify_gene_wise()
    
    # Perform DESeq
    cd.controller_deseq.compare_with_deseq()
