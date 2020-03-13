import numpy as np
import matplotlib.pyplot as plt
import os
import pygeod
from scipy import interpolate

"""
Invariant point
"""
""" This gets the area of a section.
Note that this is not a true area, as it does not
take into account the non-concavity of some of the sections
but this is unimportant for the current task.
"""
def get_area_of_sec(data):
    inside_pt = np.array([np.mean(data[:,1]), np.mean(data[:,2])])
    centered = np.copy(data)
    centered[:,1] -= inside_pt[0]
    centered[:,2] -= inside_pt[1]
    sect = centered[:,1:3]
    pt0 = np.zeros((len(sect)+1,2))
    pt1 = np.zeros((len(sect)+1,2))
    pt0[:-1,:] = sect
    pt1[1:,:] = sect
    area = np.nansum(0.5 * (pt0[:,0] + pt1[:,0]) * (pt0[:,1] - pt1[:,1]))
    if (len(data) == 0) or area < 0:
        area = np.nan
    return(area)

def get_area_dict(dcty):
    print("getting area dictionary...")
    area_dict = {}
    for i in os.listdir(dcty):
        print(f' {i}          ',flush=True, end='\r')
        if i.endswith('bin'):
            data = pygeod.get_arr_from_bin(i)
            area_dict[i] = get_area_of_sect(i)
    return(area_dict)

def make_big_array(dcty):
    print("making huge array...")
    length = 0
    for i in os.listdir(dcty):
        print(f' {i}          ',flush=True, end='\r')
        if i.endswith('bin'):
            data = pygeod.get_arr_from_bin(f'{dcty}/'+i)
            length += len(data)

    big_arr = np.zeros((length, 3))

    start, stop = 0, 0
    for i in os.listdir(dcty):
        print(f' {i}          ', flush = True, end = '\r')
        if i.endswith('bin'):
            data = pygeod.get_arr_from_bin(f'{dcty}/'+i)
            stop += len(data)
            big_arr[start:stop, 0] = data[:, 1]
            big_arr[start:stop, 1] = data[:, 2]
            big_arr[start:stop, 2] = get_area_of_sec(data)
            start += len(data)

    return(big_arr)


def get_invariant_point(dcty):
    print("getting invariant point...")
    
    big_arr = make_big_array(dcty)[::2000]
    area_function = interpolate.interp2d(big_arr[:,0], big_arr[:,1], big_arr[:,2])
    x = np.linspace(np.nanmin(big_arr[:,0]),
                    np.nanmax(big_arr[:,0]),
                    int(1e3))

    y = np.linspace(np.nanmin(big_arr[:,1]),
                    np.nanmax(big_arr[:,1]),
                    int(1e3)) 
    
    area_map = area_function(x, y)

    minidx = np.argmin(area_map)
    #x_val = x[minidx[0]]
    #y_val = y[minidx[1]]

    return(area_function)

"""
Rotation curve
"""
def get_rotation_number(data, u0):
    sect = data[:,1:3]
    angles = np.zeros(len(sect))
    
    pt0 = np.zeros((len(sect)+1,2))
    pt1 = np.zeros((len(sect)+1,2))
    pt0[:-1,:] = sect
    pt1[1:,:] = sect
    v0 = pt0 - u0
    v1 = pt1 - u0
    angles = angle_between(v0,v1)[1:-1]
    angles = angles[np.nonzero(angles)[0]]
    return(angles.mean()/(2*np.pi))

def get_rot_curve(data, invpt, curve_list):
    try:
        r1, r2 = get_radius_of_section(data)
    except:
        r1, r2 = np.nan, np.nan
    rot_num = get_rotation_number(data, invpt)
    curve_list.append([r1, rot_num])
    curve_list.append([r2, rot_num])

def get_radius_of_section(section, tol = 0.0001):
    pr0_points = section[np.abs(section[:,2]) < tol]
    pr0_radii = pr0_points[:,1]
    bins = np.arange(np.min(pr0_radii), np.max(pr0_radii)+tol, tol)
    inds = np.digitize(pr0_radii, bins)
    rads = bins[inds]
    left, right = rad_finder(rads)
    return(np.array([left, right]))

def rad_finder(rads):
    mean = np.nanmean(rads)
    left = rads[rads < mean]
    right = rads[rads > mean]
    return(np.nanmean(left), np.nanmean(right))

"""
Misc.
"""
def data_generator(dcty, ending = 'bin'):
    for i in os.listdir(dcty):
        if i.endswith(ending):
            data = pygeod.get_arr_from_bin(dcty+'/'+i)
            if len(data) > 10:
                yield (i, data)

def unit_vector(vector):
    norm = np.linalg.norm(vector, axis=1)
    newvec = np.copy(vector)
    newvec[:,0] /= norm
    newvec[:,1] /= norm
    return newvec



def angle_between(v1, v2):
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        v_combo = np.zeros((len(v1_u), 2, 2))
        v_combo[:,0,:] = v1_u
        v_combo[:,1,:] = v2_u
        detwhere = np.where(np.linalg.det(v_combo) >= 0)
        notdetwhere = np.setdiff1d(np.arange(len(v_combo)), detwhere)
        angles = np.zeros(len(v1))
        angles[detwhere] = 2*np.pi-np.arccos(np.clip(np.sum(v1_u[detwhere]*v2_u[detwhere], axis=1), -1.0, 1.0))
        angles[notdetwhere] = np.arccos(np.clip(np.sum(v1_u[notdetwhere]*v2_u[notdetwhere], axis=1), -1.0, 1.0))
        return(angles)


"""
Plotting
"""
def plot_poin(filedata, ax, xlim, poin_ylim, s=0.0001, **kwargs):
    if not ax:
        fig, ax = plt.subplots(111)

    if not poin_ylim: poin_ylim = (-np.inf, np.inf)
    if not xlim: xlim = (-np.inf, np.inf)
    
    filedata = filedata[:,1:3]
    
    #x cut
    filedata = filedata[filedata[:,0] > xlim[0]]
    filedata = filedata[filedata[:,0] < xlim[1]]

    #ycut
    filedata = filedata[filedata[:,1] > poin_ylim[0]]
    filedata = filedata[filedata[:,1] < poin_ylim[1]]

    ax.scatter(filedata[:,0], filedata[:,1], s=s, **kwargs)
    return(ax)

def make_combo_plot(dcty,
                    invpt,
                    xlim = None,
                    poin_ylim = None,
                    rot_ylim = None,
                    fig = None,
                    poinplot = None,
                    rotplot = None,
                    s = 0.0001,
                    ex = 0,
                    **kwargs):
    
    if not fig:
        fig, (poinplot, rotplot) = plt.subplots(2,1,sharex=True)
        fig.subplots_adjust(hspace=0)
    
    datas = data_generator(dcty)
    
    curve_list = []
    
    for package in datas:
        filename, data = package
        print(f'\r {filename}: making poincare map...         ', flush = True, end='')
        plot_poin(data, poinplot, xlim, poin_ylim, s, **kwargs)
        print(f'\r {filename}: calculating rotation number...', flush = True, end='')
        get_rot_curve(data, invpt, curve_list)

    curve_list = np.array(curve_list)
    curve_list = curve_list[curve_list[:,0].argsort()]

    rotplot.plot(curve_list[:,0], curve_list[:,1]+ex, marker='o', ms=2, **kwargs)
    if xlim:
       poinplot.set_xlim(xlim[0], xlim[1])
    if poin_ylim:
       poinplot.set_ylim(poin_ylim[0], poin_ylim[1])
    if rot_ylim:
       rotplot.set_ylim(rot_ylim[0]+ex, rot_ylim[1]+ex)
    rotplot.set_ylabel(r'$\nu_\theta$',size=20)
    poinplot.set_ylabel(r'$P_r$', size=20)
    rotplot.set_xlabel(r'$r$', size=20)
    poinplot.set_title(make_title(dcty))
    return(fig, poinplot, rotplot)

def make_title(dcty):
    for i in os.listdir(dcty):
        if i.endswith('ini'):
            config_file = i
            break
    print(config_file)
    config = pygeod.get_config(dcty+'/'+config_file)
    configph = config['PHYSICAL PARAMETERS']
    title = r"$M = {m}, a = {a}, L = {jz}, E = {e}, \zeta = {z}$".format(m = configph['MASS'], a = configph['SPIN'], jz = configph['JZ'], e = configph['ENERGY'], z = configph['ZETA'])
    return(title)
