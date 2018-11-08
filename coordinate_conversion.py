import pandas as pd
import sunpy.coordinates
from astropy.coordinates import SkyCoord
from math import isnan
import astropy.units as u 


def calc_hc_coord(row):
    """
    Given a row containing at least "LAT_FWT" and "LON_FWT", the coordinates in heliographic stonyhurst, and
    "T_REC", the time of observation, return a SkyCoord object containing data in heliocentric coordinates.
    Extract the values from the returned value through ret.x.value, ret.y.value, and ret.z.value (center of sun is
    (0,0,0))
    :param row: Pandas series containing at least "LAT_FWT", "LON_FWT" and "T_REC"
    :return: SkyCoord object in heliocentric coordinates
    """
    lat, lon = row["LAT_FWT"], row["LON_FWT"]
    try:
        hgs = SkyCoord(lon*u.deg, lat*u.deg, obstime = row["T_REC"], frame = "heliographic_stonyhurst")
        ret = hgs.transform_to("heliocentric")
    except:
        ret = None

    return ret


def update_full_df(df, out_file):
    """
    Update dataframe in batches. Will check to see if the last value in a batch is defined before running
    coordinate update on full batch. Saves to out_file after each batch, and at the end.
    :param df: Dataframe containing SHARPs data
    :param out_file: Filename to output modified dataframe
    """
    df_size = df.shape[0]

    batch_count = 1
    batch_size = 10000

    while batch_size * batch_count+1 < df_size:
        end = (batch_count * batch_size) - 1
        if end > df_size:
            print("Last batch")
            end = df_size - 1
        start = (batch_count * batch_size) - batch_size
        print("Starting Batch {}, {} rows remaining".format(batch_count, df_size - start))

        # print("Start: {}, End: {}".format(start, end))
        # If last row in batch hasn't been updated, update whole batch
        # print("Data: {}".format(df.loc[df.index[end]],"hc_x"))
        if isnan(df.loc[df.index[end], "hc_x"]):
            # hc_data = calc_hc_coord(df.loc[df.index[end]])
            # df.loc[df.index[end], ["hc_x", "hc_y", "hc_z"]] = [float(hc_data.x.value), hc_data.y.value,
            #                                                        hc_data.z.value]
            ind = start
            while ind <= end:
                row = df.loc[df.index[ind]]
                hc_data = calc_hc_coord(row)

                df.loc[df.index[ind], ["hc_x", "hc_y", "hc_z"]] = [float(hc_data.x.value), hc_data.y.value, hc_data.z.value]

                ind +=1
        batch_count += 1
        if df.shape[0] == df_size:
            df.to_csv(out_file)
        # print("Shape of dataframe: {}".format(df.shape))

    # last = (batch_count * batch_size) - batch_size
    #
    # # TODO: FIX THIS
    # df.at[df.index[last]:, ["hc_x", "hc_y", "hc_z"]] = df.loc[df.index[last]:].apply(calc_hc_coord, axis=1)
    # df.to_csv(out_file)


def update_missing_values(df, out_file):
    index = df['hc_x'].index[df['hc_x'].apply(isnan)]
    print("Missing {} Values, updating now".format(len(index)))
    for ind in index:
        row = df.loc[df.index[ind]]
        hc_data = calc_hc_coord(row)

        df.loc[df.index[ind], ["hc_x", "hc_y", "hc_z"]] = [float(hc_data.x.value), hc_data.y.value, hc_data.z.value]
    df.to_csv(out_file)

if __name__ == "__main__":
    
    try:
        data = pd.read_csv("modified_data.csv")
    except:
        data = pd.read_csv("data.csv")
        data["hc_x"] = float("NaN")
        data["hc_y"] = float("NaN")
        data["hc_z"] = float("NaN")

    df_size = data.shape[0]
    print("Size of full dataset: {}".format(df_size))
    output_filename = "modified_data.csv"

    # update_full_df(data, output_filename)
    update_missing_values(data, output_filename)
    
    
    
