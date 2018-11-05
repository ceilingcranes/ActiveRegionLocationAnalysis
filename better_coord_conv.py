import pandas as pd
import sunpy.coordinates
from astropy.coordinates import SkyCoord
from math import isnan
import astropy.units as u 

def test_func(row):
    return "test"

def calc_hc_coord(row):
    lat, lon = row["LAT_FWT"], row["LON_FWT"]
    try:
        hgs = SkyCoord(lon*u.deg, lat*u.deg, obstime = row["T_REC"], frame = "heliographic_stonyhurst")
        ret = hgs.transform_to("heliocentric")
    except:
        ret = None
    return ret
    return pd.Series(ret.x.value, ret.y.value, ret.z.value)
    
if __name__ == "__main__":
    
#    try:
    data = pd.read_csv("modified_data.csv")
 #   except:
  #      data = pd.read_csv("data.csv")
   #     data["hc_x"] = float("NaN")
    #    data["hc_y"] = float("NaN")
     #   data["hc_z"] = float("NaN")

    batch_size = 10000
    df_size = data.shape[0]
    print("Size of full dataset: {}".format(df_size))
    batch_count = 1
    output_filename = "modified_data.csv"
    nan_count = 0

    while batch_size*batch_count < df_size:
        end = (batch_count*batch_size) - 1
        start = (batch_count*batch_size) - batch_size
        # print("Start: {}, End: {}".format(start, end))
        # If last row in batch hasn't been updated, update whole batch
        if isnan(data.loc[data.index[end],"hc_x"]):
            ind=start
            while ind < end:
                row = data.loc[data.index[ind]]
                if isnan(row["LON_FWT"]) or isnan(row["LAT_FWT"]):
                    nan_count += 1

                hc_data = calc_hc_coord(row)

                data.loc[data.index[ind], ["hc_x", "hc_y", "hc_z"]] = [float(hc_data.x.value), hc_data.y.value, hc_data.z.value]

                ind +=1
        batch_count += 1
        data.to_csv(output_filename)
        print("Shape of dataframe: {}".format(data.shape))
        print("Batch {} Complete, {} rows remaining".format(batch_count, df_size-end))

    last = (batch_count*batch_size) - batch_size  
    data.at[data.index[last]:, ["hc_x", "hc_y", "hc_z"]] = data.loc[data.index[last]:].apply(calc_hc_coord, axis=1)
    data.to_csv(output_filename)
    print("Total NaN coordinate values: ", nan_count)
    

    
    
    
