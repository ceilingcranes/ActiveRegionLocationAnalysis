import pandas as pd
import sunpy.coordinates
from astropy.coordinates import SkyCoord
import astropy.units as u
import pickle
import numpy as np


# def calculate_distance_graph(coord_objs):

# def calc_dist_sphere(x1,y1,z1,x2,y2,z2):


def calculate_distance_graph(series):
    # get all indices that match, get skycoord objects
    skycoord_list = []
    count=series.shape[0]
    # Error because this line only returns 1 row.
    for ind in range(count):
        row = series.loc[series.index[ind]]
        # [SkyCoord(row["FWT_LON"]*u.deg, row["FWT_LAT"]*u.deg, obstime = row["T_REC"], frame = "heliographic_stonyhurst")\
        #              for row in df.loc[df["T_REC"] == t_rec]]
        # skycoord_list.append(SkyCoord(row["hc_x"]*u.deg,
        #                           row["hc_y"]*u.deg,
        #                               row["hc_z"]*u.deg,
        #                           obstime = row["T_REC"],
        #                           frame = "heliocentric"))
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

    # print(dist_graph)
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

    return thresholded_graph

# def get_connected_components(graph):
    # Given a graph in the form of an adjacency list with a dictionary, get the total number of components


if __name__ == "__main__":
    df = pd.read_csv("distance_data.csv")
    # print(df.columns)
    # print(df.drop(columns = ['Unnamed: 0', 'Unnamed: 0.1', 'Unnamed: 0.1.1', 'Unnamed: 0.1.1.1',
    #    'Unnamed: 0.1.1.1.1', 'Unnamed: 0.1.1.1.1.1'], inplace = True))
    # print(df.columns)

    batch_size = 1000
    graph_filename = "dist_networks.obj"
    output_filename = "distance_data.csv"
    added_columns = ["closest_dist"]
    save_data = False
    for c in added_columns:
        if c not in df.columns:
            df[c] = float('NaN')
    try:
        with open(graph_filename, 'r') as filename:
            all_distance_graphs = pickle.load(filename)
    except:
        all_distance_graphs = dict()

    out_file = open(graph_filename, 'w')

    ind = 1

    while ind < df.shape[0]:
        # For each row, check to see if feature values have been created. If not, update all rows for that time
        # With the distance calculation.

        # Will expect a dict of index -> list of feature updates
        t_rec = df.loc[df.index[ind], "T_REC"]
        vals = df.loc[df["T_REC"]==t_rec]
        # print(df.loc[df["T_REC"]==t_rec,"closest_dist"])
        if np.isnan(df.loc[df.index[ind],"closest_dist"]):
            if not t_rec in all_distance_graphs.keys():
                dists = calculate_distance_graph(vals)

                all_distance_graphs[t_rec] = dists
            else:
                dists = all_distance_graphs[t_rec]

            for key, ind in enumerate(vals.index):
                # print(key, ind)
                # print(vals.loc[vals.index[0]])
                # print("Distance to 0: ", dists[key][0])
                closest = min(dists[key])
                df.loc[df.index[ind], "closest_dist"] = closest

        if ind % batch_size == 0 and save_data:
            pickle.dump(all_distance_graphs, out_file)
            df.to_csv(output_filename)
            print("Batch {} Finished".format(ind/batch_size))
        ind += 1

    out_file.close()
    df.to_csv(output_filename)