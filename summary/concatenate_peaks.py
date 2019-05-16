import pandas as pd
import glob



if __name__ == '__main__':
    paths = glob.glob('peaks/*')

    dfs = []
    for path in paths:
        df = pd.read_csv(path)
        if df.shape[0] > 1:  # TODO this is wrong. 
         dfs.append(df)

    df_con = pd.concat(dfs)
    df_con.to_csv("all_peaks.csv", index=False)
