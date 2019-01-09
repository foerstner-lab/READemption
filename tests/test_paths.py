import os
import shutil
import sys
sys.path.append("./tests")
import paths_data as pd


def setup_function():
    pd.data_paths()


def teardown_function(function):
    remove_folder_if_exists(pd.test_folder)


def remove_folder_if_exists(folder):
        if os.path.exists(folder):
            shutil.rmtree(folder)


def test_set_folder_names():
    pd.paths._set_folder_names()
    for folder_name in pd.folder_names:
        assert folder_name != ''
        assert pd.folder_names.count(folder_name) == 1


def test_set_folder_names_with_base_path():
    pd.paths._set_folder_names()
    for folder_name in pd.folder_names:
        assert folder_name != ''
        assert pd.folder_names.count(folder_name) == 1


def test_set_files():
    pd.paths._set_static_files
    for file_name in pd.static_files:
        assert file_name != ''
        assert pd.static_files.count(file_name) == 1


def test_get_sorted_folder_content():
    create_folder_with_files(pd.test_folder, pd.test_lib_names)
    assert pd.paths._get_sorted_folder_content(
        pd.test_folder) == sorted(pd.test_lib_names)


def test_required_folders():
    assert len(pd.paths.required_folders()) == 25


def test_get_read_files():
    folder = pd.test_folder + '/input/reference_sequences'
    create_folder_with_files(folder, pd.test_lib_names)
    assert pd.paths.get_ref_seq_files() == sorted(pd.test_lib_names)


def test_get_ref_seq_files():
    folder = pd.test_folder + '/input/reference_sequences'
    create_folder_with_files(folder, pd.test_lib_names)
    assert pd.paths.get_ref_seq_files() == sorted(pd.test_lib_names)


def test_get_annotation_files():
    folder = pd.test_folder + '/input/annotations'
    create_folder_with_files(folder, pd.test_lib_names)
    assert pd.paths.get_annotation_files() == sorted(pd.test_lib_names)


def test_set_read_files_dep_file_lists():
    pd.paths.set_read_files_dep_file_lists_single_end(
        pd.test_files, pd.test_lib_names)
    for path_list in [
            pd.paths.read_paths, pd.paths.processed_read_paths,
            pd.paths.primary_read_aligner_sam_paths,
            pd.paths.unaligned_reads_paths]:
        assert isinstance(path_list, list)


def test_path_list_without_appendix():
    assert pd.paths._path_list(pd.test_folder, pd.test_lib_names) == [
        '/tmp/test/foo', '/tmp/test/bar']


def test_path_list_with_appendix():
    assert pd.paths._path_list(
        pd.test_folder, pd.test_lib_names, appendix=".boing") == [
            '/tmp/test/foo.boing', '/tmp/test/bar.boing']


def test_set_ref_seq_paths():
    pd.paths.set_ref_seq_paths(pd.test_files)
    assert pd.paths.ref_seq_paths == ["{}/{}".format(
        pd.paths.ref_seq_folder, file) for file in pd.test_files]

    
def create_folder_with_files(folder, file_names):
    os.makedirs(folder)
    for file_name in file_names:
        open("{}/{}".format(folder, file_name), "w").close()
