import numpy as np
import pandas as pd

def dfhkl2dfhklaxes(df, min_intensity, factory, geometry, detector, sample, user):
    rows = []
    engines = factory.create_new_engine_list()
    engines.init(geometry, detector, sample)
    engines.get()
    engine_hkl = engines.engine_get_by_name("hkl")
    for refl in df.itertuples(index=False):
        h = refl.h
        k = refl.k
        l = refl.l
        d = refl.d
        inten = refl.intensity
        if inten > min_intensity:
            try:
                solutions = engine_hkl.pseudo_axis_values_set([h, k, l], user)
                for _, item in solutions.items():
                    read = item.geometry_get().axis_values_get(user)
                    if read is not None:
                        rows.append({
                            'h': h,
                            'k': k,
                            'l': l,
                            'd': d,
                            'intensity': inten,
                            'omega': read[0],
                            'chi': read[1],
                            'phi': read[2],
                            'tth': read[3]
                        })
            except Exception as e:
                print(e)
    new_df = pd.DataFrame(rows, columns=['h', 'k', 'l', 'omega', 'chi', 'phi', 'tth', 'd', 'intensity'])
    return new_df
