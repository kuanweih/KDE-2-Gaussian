import os
import errno
import pandas as pd
from typing import List


def create_dir(dir_name: str):
    """ Create directory with a name 'dir_name' """
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def df_concat(paths: List[str]) -> pd.DataFrame:
    """ Concatenate multiple pandas dataframe.

    : paths : a list of paths for csv files
    : return : concatenated dataframe
    """
    dfs = [pd.read_csv(path) for path in paths]
    df_con = pd.concat(dfs)
    return  df_con
