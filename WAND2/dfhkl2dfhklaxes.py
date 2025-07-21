import numpy as np
import pandas as pd
from tqdm import tqdm
import logging

logger = logging.getLogger(__name__)

def dfhkl2dfhklaxes(df, min_intensity, factory, geometry, detector, sample, user):
    logger.info("Starting dfhkl2dfhklaxes with %d reflections", len(df))
    rows = []
    new_df = pd.DataFrame(columns=['h', 'k', 'l', 'mu', 'omega', 'chi', 'phi', 'gamma','delta', 'd', 'intensity']) 
    engines = factory.create_new_engine_list()
    engines.init(geometry, detector, sample)
    engines.get()
    engine_hkl = engines.engine_get_by_name("hkl")
    engine_hkl.current_mode_set('lifting_detector_mu') # TODO CHECK THIS
    #axes = geometry.axis_names_get()
    #for axis in axes:
    #    tmp = geometry.axis_get(axis)
    #    #if (axis=='omega') or (axis=='chi') or (axis=='phi'):
    #    if (axis=='chi'):
    #        tmp.min_max_set(-0.01, 0.01, user)
    #        geometry.axis_set(axis, tmp)
    num_refl = len(df)
    print(f"Searching through {num_refl} reflections...")
    found = 0
    not_found = 0
    df = df[df['intensity'] > min_intensity]
    for idx, refl in tqdm(df.iterrows(), total=len(df), desc="Reflections"):
        h = refl['h']
        k = refl['k']
        l = refl['l']
        d = refl['d']
        inten = refl['intensity']
        try:
            solutions = engine_hkl.pseudo_axis_values_set([h,k,l], user)
            # similar to apply_axes_solns in hkl.py
            for i, item in enumerate(solutions.items()):
                read = item.geometry_get().axis_values_get(user)
                if read is not None:
                    rows.append({'h':h, \
                                 'k':k, \
                                 'l':l, \
                                 'd':d, \
                                 'intensity':inten, \
                                 'mu':read[0], \
                                 'omega':read[1], \
                                 'chi':read[2], \
                                 'phi':read[3], \
                                 'gamma':read[4], \
                                 'delta':read[5]})
                    found += 1
        except Exception as e:
            #logger.exception(f"Exception for hkl=({h},{k},{l}): {e}")
            logger.info(f"Exception for hkl=({h},{k},{l}): {e}")
            not_found += 1
            #print(e)
    new_df = pd.DataFrame(rows, columns=['h', 'k', 'l', 'mu', 'omega', 'chi', 'phi', 'gamma', 'delta', 'd', 'intensity'])
    foundrefl = num_refl-not_found
    print(f"found {found} motor positions in {foundrefl} reflections. Did not find positions for {not_found} reflections.")
    logger.info("Completed dfhkl2dfhklaxes. Output DataFrame has %d rows", len(new_df))
    if new_df is not None:
        return new_df
    else:
        print("empty dataframe, something went wrong")
        return
