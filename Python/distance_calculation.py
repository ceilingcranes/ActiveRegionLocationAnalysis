import pandas as pd
import sunpy.coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
import pickle
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


def calculate_distance_graph(series):
    # get all indices that match, get skycoord objects
    skycoord_list = []
    count=series.shape[0]
    # Error because this line only returns 1 row.
    for ind in range(count):
        row = series.loc[series.index[ind]]
        skycoord_list.append(SkyCoord(row["LON_FWT"]*u.deg,
                                      row["LAT_FWT"]*u.deg,
                                      obstime = row["T_REC"],
                                      frame = "heliographic_stonyhurst"))
    dist_graph = dict()
    # Initialize dict
    for i in range(count):
        dist_graph[i] = np.zeros(count)

    for i in range(count):
        for j in range(i, series.shape[0]):
            if i==j:
                dist_graph[i][j] = float('inf')
                dist_graph[j][i] = float('inf')
            else:
                # Calculate the Great-Circle distance
                dist = skycoord_list[i].separation(skycoord_list[j])
                dist_graph[i][j] = dist.degree
                dist_graph[j][i] = dist.degree

    # print(skycoord_list)
    return dist_graph


def thresholded_distance_graph(dist_graph, threshold):
    # Given a distance graph and a threshold, remove all connections with distance greater than threshold and return
    # a new graph
    thresholded_graph = dict()
    # initialize graph

    for key in dist_graph.keys():
        thresholded_graph[key] = np.zeros(len(dist_graph[key]))
        for ind, val in enumerate(dist_graph[key]):
            if val <= threshold:
                thresholded_graph[key][ind] = val
            else:
                thresholded_graph[key][ind] = float('inf')
        thresholded_graph[key][key] = float('inf')
    return thresholded_graph


def get_connected_components(graph):
    # Given a graph in the form of an adjacency list with a dictionary, return a 2d list containing the indices of each
    # component. I.e. if there were 2 components in a 6 node graph, with 3 nodes in eachcomponent, this would return a
    # 2x3 2d array.
    def _visit_val(index):
        visited[index] = True
        for i in range(len(graph[index])):
            if graph[index][i] != 0 and i != index and not visited[i]:
                component_list.append(i)
                _visit_val(i)

    visited = [False for i in graph.keys()]
    components = []
    component_list = []
    for ind in graph.keys():
        if not visited[ind]:
            component_list = [ind]
            _visit_val(ind)
            components.append(component_list)

    return components


def flare_in_component(component_list, rows):
    # Returns a list of indices that are part of a component that had at least one node flare at least once in the 24hrs
    # following t_rec, as given by the dataframe rows in rows.
    flare_indices = []
    for component in component_list:
        has_flare = False
        for ind in component:
            if rows.loc[rows.index[ind],"any_flare_in_24h"] != 0:
                has_flare = True
        if has_flare:
            flare_indices = np.append(flare_indices, component)
    return flare_indices


def get_mean_dist(dist_graph):
    sum = 0
    print("dist graph: ", dist_graph)
    counter = 0
    graph = dist_graph.copy()
    for key in graph.keys():
        distances = np.delete(graph[key],key)

        for d in distances[key:]:
            if d != 0:
                sum += d
                counter += 1
        # if len(distances) != 0:
        #     sum += np.sum(distances)
    return sum / counter

def get_min_max_dist(dists):
    if len(dists.keys()) == 1:
        return 0,0
    maxd = 0
    mind = float('inf')

    for key in dists.keys():
        vals = dists[key][dists[key] != np.inf]
        max_k = np.nanmax(vals)
        min_k = np.min(vals)
        if np.isfinite(max_k) and max_k > maxd:
            maxd = max_k
        if np.isfinite(min_k) and min_k < mind:
            mind = min_k

    return (mind, maxd)

def create_nxgraph(distance_graph):
    # Given a distance graph that uses a dictionary, as is created by calculate_distance_graph, convert that
    # Dictionary to a networkX graph object and return it.
    graph = nx.Graph()

    # Add nodes - Each node is just an index
    graph.add_nodes_from(range(len(distance_graph[0])))

    for node in distance_graph:
        for ind, dist in enumerate(distance_graph[node]):
            if np.isfinite(dist):
                graph.add_edge(node, ind, distance=dist)


    # nx.draw(graph)
    # plt.show()
    return graph

def add_components(dists, thresh, thresh_name, vals):
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("Creating values from {} using threshold {}".format(dists, thresh))
    thresh_graph = thresholded_distance_graph(dists, thresh)
    print("thresholded graph: {}".format(thresh_graph))
    nxgraph = create_nxgraph(thresh_graph)
    centrality = calc_centrality(nxgraph)
    components = get_connected_components(thresh_graph)
    print("Components: {}".format(components))
    print("++++++++++++++++++++++++++++++++++++++++++++++++")
    flare_components = flare_in_component(components, vals)
    for key, dist_ind in enumerate(vals.index):
        # print(key, ind)
        # print(vals.loc[vals.index[0]])
        # print("Distance to 0: ", dists[key][0])

        nodes_in_component = len([i for i in thresh_graph[key] if np.isfinite(i)])
        # print("+++++++++++++++++++++++++++++++++=")
        # print(nodes_in_component)
        # print("-------")
        # print(mean_thresh)
        # print("+++++++++++++++++++++++++++++++++=")

        # print("Thresholded at {}, graph: {}".format(thresh, mean_thresh))
        # Create feature for number of components with a mean threshold at this time step.

        df.loc[df.index[dist_ind], "{}_nodes_in_component".format(thresh_name)] = nodes_in_component
        df.loc[df.index[dist_ind], "{}_connected_components_count".format(thresh_name)] = len(components)
        df.loc[df.index[dist_ind], "{}_total_nodes".format(thresh_name)] = len(dists.keys())
        df.loc[df.index[dist_ind], "{}_mean_distance".format(thresh_name)] = thresh
        df.loc[df.index[dist_ind], "{}_flare_in_component".format(thresh_name)] = (dist_ind in flare_components)
        df.loc[df.index[dist_ind], "{}_eig_centrality".format(thresh_name)] = centrality[key]


def calc_centrality(nxgraph):
    # use networkx centrality algorithms
    eig_cent = nx.eigenvector_centrality(nxgraph)
    return eig_cent


def add_component_data(df, output_filename = "data/distance_data.csv", graph_filename = "data/dist_networks.obj"):
    # added_columns = ["closest_dist", "connected_components_mean_thresh", "total_nodes", "mean_distance",
    #                  "nodes_in_component", "flare_in_component", "eig_centrality"]
    save_data = True
    batch_size = 3000

    # for c in added_columns:
    #     if c not in df.columns:
    #         df[c] = float('NaN')

    try:
        with open(graph_filename, 'rb') as filename:
            all_distance_graphs = pickle.load(filename)
    except:
        all_distance_graphs = dict()
        print("Unable to load stored files.")

    out_file = open(graph_filename, 'wb')
    batch_count = 1
    ind = 226950

    while ind < df.shape[0]:
        # For each row, check to see if feature values have been created. If not, update all rows for that time
        # With the distance calculation.

        # Will expect a dict of index -> list of feature updates
        t_rec = df.loc[df.index[ind], "T_REC"]
        vals = df.loc[df["T_REC"] == t_rec]
        # print("COUNT OF ROWS WITH T_REC {}: {}".format(t_rec, len(vals)))
        # print(df.loc[df["T_REC"]==t_rec,"closest_dist"])
        # Skipping already designated lines is not working, it skips a lot more than that.
        # if np.isnan(df.loc[df.index[ind], added_columns[-1]]):
        if True:
            if not t_rec in all_distance_graphs.keys():
                # print("Values: [{}]={}".format(vals["T_REC"], vals["HARPNUM"]))
                dists = calculate_distance_graph(vals)
                # print("Full graph: {}".format(dists))
                all_distance_graphs[t_rec] = dists
            else:
                dists = all_distance_graphs[t_rec]

            # Calculate values that won't change for individual indices
            # print(vals)

            # Mean threshold
            if len(dists.keys()) > 1:
                mean_thresh = get_mean_dist(dists)
                add_components(dists, mean_thresh, "mean", vals)
            # try:
            # # Steps of 10
            #     minmax_list=[]
            #     minmax_dist = get_min_max_dist(dists)
            #     if minmax_dist[0] != minmax_dist[1]:
            #         minmax_list = [i*minmax_dist[1]/10 for i in range(0,10)]
            #         for i, thresh in enumerate(minmax_list):
            #             add_components(dists, thresh, str(i), vals)
            #
            #         for key, dist_ind in enumerate(vals.index):
            #             closest = min(dists[key])
            #             df.loc[df.index[dist_ind], "closest_dist"] = closest
            # except:
            #     print('error in thresholding')
            #     print(minmax_list)
            #     print(minmax_dist)
        if ind % batch_size == 0 and save_data:
            # if added_graph:
            #     pickle.dump(all_distance_graphs, out_file)
            #     print("Added distance graph that was not found in pickle!")
            df.to_csv(output_filename, index=False)
            print("Batch {} Finished".format(batch_count))
            batch_count += 1
        ind += 1

    out_file.close()
    df.to_csv(output_filename, index=False)


if __name__ == "__main__":
    df = pd.read_csv("data/distance_data.csv")
    # print(df.columns)
    # print(df.drop(columns = ['Unnamed: 0', 'Unnamed: 0.1', 'Unnamed: 0.1.1', 'Unnamed: 0.1.1.1',
    #    'Unnamed: 0.1.1.1.1', 'Unnamed: 0.1.1.1.1.1'], inplace = True))
    print(df.head())
    print(df.shape)
    add_component_data(df, output_filename="__thresh_data.csv")
    print(df.head())