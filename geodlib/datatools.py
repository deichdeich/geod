import numpy as np
import os

def get_arr_from_bin(fname):
    outdat = np.fromfile(fname)
    outdat = outdat.reshape(int(len(outdat)/7), 7)
    return(outdat[1:])


def get_ics_from_data(dcty):
    datagen = data_generator(dcty)
    ics = []
    for _, data in datagen:
        ics.append(data[-1, 1:-1])
    return(np.array(ics))

def data_generator(dcty, ending = 'bin'):
    for i in os.listdir(dcty):
        if i.endswith(ending):
            data = get_arr_from_bin(dcty+'/'+i)
            if len(data) > 10:
                yield (i, data)
