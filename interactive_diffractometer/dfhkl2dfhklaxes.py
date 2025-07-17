import numpy as np
import pandas as pd

def dfhkl2dfhklaxes(df, factory, geometry, detector, sample, user):
    new_df = pd.DataFrame(columns=['h', 'k', 'l', 'omega', 'chi', 'phi', 'tth', 'd', 'intensity'])    
    engines = factory.create_new_engine_list()
    engines.init(geometry, detector, sample)
    engines.get()
    engine_hkl = engines.engine_get_by_name("hkl")
    for idx, refl in df.iterrows():
        h = refl['h']
        k = refl['k']
        l = refl['l']
        d = refl['d']
        inten = refl['intensity']
        solutions = engine_hkl.pseudo_axis_values_set([h,k,l], user)
        # similar to apply_axes_solns in hkl.py
        for i, item in enumerate(solutions.items()):
            read = item.geometry_get().axis_values_get(user)
            if read is not None:
                new_row = pd.DataFrame([{'h':h, \
                                        'k':k, \
                                        'l':l, \
                                        'd':4.5, \
                                        'intensity':inten, \
                                        'omega':read[0], \
                                        'chi':read[1], \
                                        'phi':read[2], \
                                        'tth':read[3]}])
                new_df = pd.concat([new_df, new_row], ignore_index=True)    
    return new_df
