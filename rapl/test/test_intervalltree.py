import unittest
import sys
sys.path.append(".")
from libs.intervaltree import Interval, TreeNode, IntervalTree

class TestInterval(unittest.TestCase):
    
    def test_interval(self):
        interval = Interval(1, 20)
        self.assertEqual(interval.start, 1)
        self.assertEqual(interval.end, 20)

    def test_start_higher_than_end(self):
        self.assertRaises(Exception, Interval, 20, 1)

class TestTreeNode(unittest.TestCase):
    
    def test_tree_node(self):
        tree_node = TreeNode(
            100, ["fake_interval"], "fake_node_1", "fake_node_2")
        self.assertEqual(tree_node.center_pos, 100)
        self.assertEqual(tree_node.intervals_center, ["fake_interval"])
        self.assertEqual(tree_node.node_left, "fake_node_1")
        self.assertEqual(tree_node.node_right, "fake_node_2")

class TestIntervalTree(unittest.TestCase):
    
    def setUp(self):
        self.interval_tree = IntervalTree()

    def test_fill(self):
        intervals = [
            Interval(2, 6), Interval(3, 9), Interval(14, 20),   
            Interval(22, 60), Interval(10, 99)]
        self.interval_tree.fill(intervals)
        self.assertEqual(type(self.interval_tree.root_node), TreeNode)
        self.assertEqual(self.interval_tree.root_node.center_pos, 10)

    def test_sort_intervals(self):
        intervals = [
            Interval(3, 4), Interval(2, 10),
            Interval(1, 3), Interval(15, 20)]
        intervals = self.interval_tree._sort_intervals(intervals)
        self.assertEqual(intervals[0].start, 1)
        self.assertEqual(intervals[0].end, 3)
        self.assertEqual(intervals[1].start, 2)
        self.assertEqual(intervals[1].end, 10)
        self.assertEqual(intervals[2].start, 3)
        self.assertEqual(intervals[2].end, 4)
        self.assertEqual(intervals[3].start, 15)
        self.assertEqual(intervals[3].end, 20)

    def test_center_of_intervals_equal_length(self):
        intervals = [
            Interval(1, 4), Interval(4, 8), Interval(11, 24), Interval(14, 20)]
        self.assertEqual(self.interval_tree._center_pos(intervals), 11)
        
    def test_center_of_intervals_odd_length(self):
        intervals = [
            Interval(1, 4), Interval(4, 8), Interval(3, 8), Interval(11, 24), 
            Interval(14, 20)]
        self.assertEqual(self.interval_tree._center_pos(intervals), 4)

    def test_generate_node_single_intervale(self):
        intervals = [Interval(1, 10)]
        tree_node = self.interval_tree._generate_node(intervals)
        self.assertEqual(type(tree_node.intervals_center), list)
        self.assertIsNone(tree_node.node_left)
        self.assertIsNone(tree_node.node_right)

    def test_generate_node_not_overlapping(self):
        intervals = [
            Interval(1, 10), Interval(5, 10), Interval(11, 20),   
            Interval(25, 50), Interval(80, 100)]
        tree_node = self.interval_tree._generate_node(intervals)
        self.assertEqual(tree_node.intervals_center[0].start, 11)
        self.assertEqual(tree_node.intervals_center[0].end, 20)
        self.assertEqual(tree_node.node_left.center_pos, 5)
        self.assertEqual(tree_node.node_right.center_pos, 80)
        self.assertEqual(len(tree_node.intervals_center), 1)

    def test_generate_node_with_two_identical_intervals(self):
        intervals = [
            Interval(1, 10), Interval(11, 15), Interval(11, 15), 
            Interval(25, 50)]
        tree_node = self.interval_tree._generate_node(intervals)
        self.assertEqual(tree_node.intervals_center[0].start, 11)
        self.assertEqual(tree_node.intervals_center[0].end, 15)
        self.assertEqual(tree_node.node_left.center_pos, 1)
        self.assertEqual(tree_node.node_right.center_pos, 25)
        # Same intervals will be same multiple time (but later after
        # querying only returned once)
        self.assertEqual(len(tree_node.intervals_center), 2)
        self.assertEqual(tree_node.intervals_center[0].start, 11)
        self.assertEqual(tree_node.intervals_center[1].start, 11)

    def test_search_point_overlaps(self):
        intervals = [Interval(1, 10)]
        self.interval_tree.fill(intervals)
        overlapping_intervals = self.interval_tree.search_point_overlaps(
            5, self.interval_tree.root_node)
        self.assertEqual(len(overlapping_intervals), 1)
        self.assertEqual(type(overlapping_intervals[0]), Interval)

    def test_search_point_overlaps(self):
        intervals = [Interval(1, 10), Interval(5, 10), 
                     Interval(11, 20), Interval(25, 50)]
        self.interval_tree.fill(intervals)
        self.assertEqual(len(self.interval_tree.search_point_overlaps(
                1, self.interval_tree.root_node)), 1)
        hit_intervals = self.interval_tree.search_point_overlaps(
            5, self.interval_tree.root_node)
        self.assertEqual(len(hit_intervals), 2)
        self.assertEqual(hit_intervals[0].start, 1)
        self.assertEqual(hit_intervals[0].end, 10)
        self.assertEqual(hit_intervals[1].start, 5)
        self.assertEqual(hit_intervals[1].end, 10)

    def test_search_point_overlaps(self):
        intervals = [Interval(1, 10), Interval(5, 10), 
                     Interval(11, 20), Interval(25, 50)]
        self.interval_tree.fill(intervals)
        self.assertEqual(len(self.interval_tree.search_point_overlaps(
                1, self.interval_tree.root_node)), 1)
        hit_intervals = self.interval_tree.search_point_overlaps(
            5, self.interval_tree.root_node)
        self.assertEqual(len(hit_intervals), 2)
        self.assertListEqual(
            sorted([(interval.start, interval.end) 
             for interval in hit_intervals]), 
            [(1, 10), (5, 10)])

    def test_search_interval_overlaps(self):
        intervals = [Interval(1, 10), Interval(5, 10), 
                     Interval(11, 20), Interval(25, 50)]
        self.interval_tree.fill(intervals)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(1, 2))), 1)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(1, 10))), 2)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(1, 20))), 3)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(1, 50))), 4)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(1, 24))), 3)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(11, 24))), 1)
        self.assertEqual(
            len(self.interval_tree.search_interval_overlaps(Interval(11, 25))), 2)

if __name__ == "__main__":
	unittest.main()
