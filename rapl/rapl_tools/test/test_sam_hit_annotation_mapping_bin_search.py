import unittest
import sys
sys.path.append('.')
sys.path.append('..')
from segemehl_hit_annotation_mapping_bin_search import SegemehlHitAnnotationMapping

class OptionsMock(object):
    min_overlap = 1
    output_file = "x"

class TestSegemehlHitAnnotationMapping(unittest.TestCase):

    def setUp(self):
        options = OptionsMock()
        args = ["humpy-dumpy.map", "boing.ppt", ]
        self.segemehl_hit_annotation_mapping = SegemehlHitAnnotationMapping(
            args, options)

    def test__sort_annotation_entries_by_start(self):
        self.segemehl_hit_annotation_mapping.annotations = [
            {"start" : 50, "end": 150}, {"start" : 10, "end": 50},
            {"start" : 10, "end": 20},  {"start" : 80, "end": 120},
            {"start" : 10, "end": 40}]
        self.segemehl_hit_annotation_mapping._sort_annotation_entries_by_start()
        self.assertEqual(
            [{'start': 10, 'end': 20}, {'start': 10, 'end': 40}, 
             {'start': 10, 'end': 50}, {'start': 50, 'end': 150}, 
             {'start': 80, 'end': 120}], 
            self.segemehl_hit_annotation_mapping.sorted_annotations)

    def test___annotation_index_1(self):
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            1, 3, 3, 3, 4, 6, 9]
        self.assertEqual(
            0,
            self.segemehl_hit_annotation_mapping._annotation_index({"end" : 3})
            )
        self.assertEqual(
            0,
            self.segemehl_hit_annotation_mapping._annotation_index({"end" : 1})
            )
        self.assertEqual(
            4,
            self.segemehl_hit_annotation_mapping._annotation_index({"end" : 6})
            )

    def test___annotation_index_1(self):
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            38, 669, 2287, 3543, 4926]
        self.assertEqual(
            1,
            self.segemehl_hit_annotation_mapping._annotation_index({"end" : 2352, "start" : 2257}) #TODO
            )


    def test__search_in_list_1(self):
        """
            ------------ Annotation
        ======           Read    

        """
        segemehl_hit = {"start" : 180, "end" : 210}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])

    def test__search_in_list_2(self):
        """
            ------------ Annotation
             ======      Read    

        """
        segemehl_hit = {"start" : 220, "end" : 260}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])

    def test__search_in_list_3(self):
        """
            ------------     Annotation
                     ======  Read    

        """
        segemehl_hit = {"start" : 280, "end" : 310}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])

    def test__search_in_list_4(self):
        """
            ------------ Annotation
            ======       Read    

        """
        segemehl_hit = {"start" : 200, "end" : 220}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])

    def test__search_in_list_6(self):
        """
            ------------ Annotation
                  ====== Read    

        """
        segemehl_hit = {"start" : 250, "end" : 300}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])


    def test__search_in_list_7(self):
        """
             ------------    Annotation 1
             --------------- Annotation 2
            ======          Read    

        """
        segemehl_hit = {"start" : 480, "end" : 550}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 500, 'end': 600}, {'start': 500, 'end': 700}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])


    def test__search_in_list_8(self):
        """
          --------    Annotation
         ===========  Read    

        """
        segemehl_hit = {"start" : 180, "end" : 320}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 100, "end": 150},
            {"start" : 200, "end": 300},
            {"start" : 500, "end": 600},
            {"start" : 500, "end": 700},
            {"start" : 600, "end": 700},
            {"start" : 800, "end": 900}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 200, 'end': 300}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])


    def test__search_in_list_9(self):
        """
                 --------    Annotation
         ===========  Read    

        """
        segemehl_hit = {"start" : 1032188, "end" : 1033081}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 1031052, "end": 1031234},
            {"start" : 1031231, "end": 1031977},
            {"start" : 1033073, "end": 1033169},
            {"start" : 1033061, "end": 1033837}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 1033073, 'end': 1033169}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])


    def test__search_in_list_99(self):
        """
669,               2297
  ------------------

                ------------------------
               2287                   3546

                                   3543                  4733
                                    -----------------------
        -----------------
       2257          2352

        """
        segemehl_hit = {"start" : 2257, "end" : 2352}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 669, "end": 2297},
            {"start" : 2287, "end": 3546}, 
            {"start" : 3543, "end": 4733},
            {"start" : 4926, "end": 5825}]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'start': 669, 'end': 2297}, {'start': 2287, 'end': 3546}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])

    def test__search_in_list_98(self):
        segemehl_hit = {"start" : 1268, "end" : 1363}
        self.segemehl_hit_annotation_mapping.sorted_annotations = [
            {"start" : 1030, "end": 2100},
            {"start" : 2112, "end": 2654},
            ]
        self.segemehl_hit_annotation_mapping.sorted_annotation_starts = [
            annotation["start"] for annotation in 
            self.segemehl_hit_annotation_mapping.sorted_annotations]
        self.assertEqual(
            [{'end': 2100, 'start': 1030}],
            [x for x in self.segemehl_hit_annotation_mapping._search_in_list(segemehl_hit)])




        
if __name__ == "__main__": 
    unittest.main()
