class Interval(object):

    def __init__(self, start, end):
        if start > end:
            raise(Exception(
                     "Error: Start position (%s) is higher than end position"
                     " (%s)" % (start, end)))
        self.start = start
        self.end = end

class TreeNode(object):
    
    def __init__(self, center_pos, intervals_center, node_left, node_right):
        self.center_pos = center_pos
        self.intervals_center = intervals_center
        self.node_left = node_left
        self.node_right = node_right

class IntervalTree(object):
    """
    https://secure.wikimedia.org/wikipedia/en/wiki/Interval_tree
    https://github.com/misshie/interval-tree
    
    Notes:

    - Intervals with the same coordinates will be save multiple times
      but when queried only one of will be returned.

    """

    def fill(self, intervals):
        
        """
        The input intervals don't have to be sorted.
        """

        self.root_node = self._generate_node(intervals)

    def _generate_node(self, intervals):
        """Generate a node by splitting the interval list.
        """
        if not intervals:
            return(None)
        center_pos = self._center_pos(intervals)
        # Initiate the three sets:
        intervals_left = []
        intervals_center = []
        intervals_right = []
        # Append the intervals
        for interval in intervals:
            if interval.end < center_pos:
                intervals_left.append(interval)
            elif interval.start > center_pos:
                intervals_right.append(interval)
            else:
                intervals_center.append(interval)
        return(TreeNode(
                center_pos, intervals_center, 
                self._generate_node(intervals_left),
                self._generate_node(intervals_right)))

    def search_interval_overlaps(self, interval):
        hit_intervals = []
        for point in range(interval.start, interval.end + 1):
            hit_intervals.extend(
            self.search_point_overlaps(point, self.root_node))
        hit_intervals = list(set(hit_intervals))
        return(hit_intervals)

    def search_point_overlaps(self, point, node):
        hit_intervals = []
        for interval in node.intervals_center:
            if point >= interval.start and point <= interval.end:
                hit_intervals.append(interval)
        if point < node.center_pos and node.node_left:
            hit_intervals.extend(
                self.search_point_overlaps(point, node.node_left))
        elif point > node.center_pos and node.node_right:
            hit_intervals.extend(
                self.search_point_overlaps(point, node.node_right))
        return(list(set(hit_intervals)))
    
    def _center_pos(self, intervals):
        center_element = int(len(intervals)/2)
        return(self._sort_intervals(intervals)[center_element].start)

    def _sort_intervals(self, intervals):
        return(sorted(intervals, key=lambda interval: interval.start))

