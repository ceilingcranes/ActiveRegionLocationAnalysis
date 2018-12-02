import unittest
import distance_calculation as dc
import numpy as np

class test_distance_calcuation(unittest.TestCase):

    def test_get_connected_components(self):
        test_graph = dict()
        test_graph[0] = [float('inf'), 2, 3, 0, 0, 0]
        test_graph[1] = [2, float('inf'), 0, 0, 0, 0]
        test_graph[2] = [3, 0, float('inf'), 0, 0, 1]
        test_graph[3] = [0,0,0,float('inf'),3, 0]
        test_graph[4] = [0,0,0,3, float('inf'),0]
        test_graph[5] = [0,0,1,0,0,float('inf')]

        comp = dc.get_connected_components(test_graph)
        print(comp)
        print(len(comp))
        self.assertEqual(len(comp),2)
        self.assertEqual(len(comp[0]), 4)
        self.assertEqual(len(comp[1]), 2)
        self.assertListEqual(sorted(comp[0]), [0,1,2,5])
        self.assertListEqual(sorted(comp[1]), [3, 4])
        # self.assertEqual(comp, 2)


    def test_thresholded_distance_graph(self):
        test_graph = dict()
        for i in range(5):
            test_graph[i] = [i + dist_val for dist_val in range(5)]
            test_graph[i][i] = float('inf')

        thresh = dc.thresholded_distance_graph(test_graph, 4)

        total_edges = 0
        for key in thresh.keys():
            for ind in range(key, len(thresh[key])):
                val = thresh[key][ind]
                if val != 0 and val != float('inf'):
                    total_edges += 1

        self.assertEqual(total_edges, 6)
        self.assertListEqual(thresh[0].tolist(), [float('inf'), 1, 2, 3, 4])
        self.assertListEqual(thresh[4].tolist(), [4, 0, 0, 0, float('inf')])