import numpy as np
import pandas as pd

def hkl2dfhkl(hkl_path):
    def sci2dec(sci):
        dec = np.round(np.float64(sci), 1)
        return dec

    with open(hkl_path, "r") as f:
        lines = f.readlines()
    intensity_lines = []
    found_data_start = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("# H") and "|Fc|^2" in line:
            found_data_start = True
            continue
        if not found_data_start or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 6:
            h, k, l, mult, d, intensity_sci = parts[:6]
            intensity_dec = sci2dec(intensity_sci)
            intensity_lines.append("%3s %3s %3s %12s %8s" % (h, k, l, d, intensity_dec))
    rows = [line.split() for line in intensity_lines]

    # convert to DataFrame
    df = pd.DataFrame(rows, columns=['h', 'k', 'l', 'd', 'intensity'])

    # convert numeric columns to proper types
    df = df.astype({
        'h': int,
        'k': int,
        'l': int,
        'd': float,
        'intensity': float
    })
    print(f'dfhkl head: \n{df.head()}')
    return df
