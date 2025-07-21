import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def dfhkl2dfhklaxes(df, min_intensity, factory, geometry, detector, sample, user):
    logger.info("Starting dfhkl2dfhklaxes with %d reflections", len(df))
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
                n_solutions = solutions.n_items
                for i in range(n_solutions):
                    item = solutions.get(i)
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
                logger.exception(f"Exception for hkl=({h},{k},{l}): {e}")
    new_df = pd.DataFrame(rows, columns=['h', 'k', 'l', 'omega', 'chi', 'phi', 'tth', 'd', 'intensity'])
    logger.info("Completed dfhkl2dfhklaxes. Output DataFrame has %d rows", len(new_df))
    return new_df
