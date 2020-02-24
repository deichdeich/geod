import configparser
import os
import matplotlib.pyplot as plt
import numpy as np

homedir = os.path.expanduser('~')

"""
The main function
"""
def get_geodesic(filename,
                 init_state,
                 config_file,
                 return_data = False,
                 overwrite_header = False,
		 statelist = None,
                 printcmd = False):
    
    filename = filename.replace('~', homedir)
    config_file = config_file.replace('~', homedir)
    config = get_config(config_file)
    confint = config['INTEGRATION PARAMETERS']
    
    # make sure everything is set up nicely
    path = confint['out_dir']
    path = path.replace('~', homedir)
    
    if not header_check(path) or overwrite_header:
        print('Making header file from config file...', flush=True)
        make_header_file(path, config_file)
        print('Building code...\n', flush=True)
        os.system('make -C geod_src')
    else:
        print('Matching header file found...', flush=True)
    
    # make sure there's a place for the data to go
    make_output_dir(path, config_file, statelist)

    tol = confint['tolerance']
    h = confint['h_init']
    pointol = confint['poincare_tolerance']
    make_section = 0 if confint['make_section'] == 'False' else 1
    t0 = confint['t0']
    tmax = confint['tmax']
    
    
    # making the command
    syscmd = f"geod_src/geod {t0} {tmax} {h} {pointol} {make_section} "
    for thing in init_state:
        syscmd += f"{thing} "
    syscmd += f"{config['PHYSICAL PARAMETERS']['jz']} {path}/{filename}"
    if printcmd: print(syscmd)
    # run the command
    print("\nRunning integration...", flush = True)
    os.system(syscmd)
    
    data = 0
    if return_data:
        data = get_arr_from_bin(f'{path}/{filename}')
        
    return(data)

"""
data stuff
"""
def plot_traj(data, ax = False, **kwargs):
    from mpl_toolkits.mplot3d import Axes3D
    import mpl_toolkits.mplot3d.art3d as art3d 
    from matplotlib.patches import Circle
    r = data[:,1]
    ph = data[:,3]
    th = data[:,5]
    if not ax:
        ax = plt.axes(projection='3d')
    ax.scatter(r * np.sin(ph) * np.cos(th), r * np.sin(ph) * np.sin(th), r*np.cos(ph), **kwargs)
    ax.scatter(0,0,0, c='black', s=20)
    p = Circle((0,0), radius = 2, edgecolor='k', facecolor='none')
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p,z=0,zdir="z")
    bl = max(r)+max(r)/5
    ax.set_zlim(-bl,bl)
    ax.set_xlim(-bl,bl)
    ax.set_ylim(-bl,bl)
    return(ax)


def get_arr_from_bin(fname):
    outdat = np.fromfile(fname)
    outdat = outdat.reshape(int(len(outdat)/7), 7)
    return(outdat)


"""
Header file stuff
"""
def make_header_file(path, config_file = 'config.ini'):
    header_file = open('geod_src/definitions.h', 'w')
    config = get_config(config_file)
    confpp = config['PHYSICAL PARAMETERS']
    header_file.write(f'//{path}\n')
    for key in confpp:
        header_file.write('#define {confkey: <16} {value}\n'.format(confkey = key.upper(),
                                                                    value = confpp[key]))
    header_file.close()
    return 0    

def header_check(path):
    config_path = ''
    if 'definitions.h' in os.listdir('geod_src'):
        header_file = open('geod_src/definitions.h', 'r')
        header_text = header_file.readlines()
        if len(header_text) > 0:
            config_path = header_text[0][2:].strip('\n')
    
    retval = config_path==path

    return(retval)

"""
Misc
"""
def make_output_dir(path, config_file, statefile):
    path1 = ('/').join(path.split('/')[:-1])
    path2 = path.split('/')[-1]
    if path2 not in os.listdir(path1):
        print('Making output directory...', flush=True)
        os.mkdir(path)
        os.system(f'cp {config_file} {path1}/{path2}')
    if statefile is not None:
        os.system(f'cp {statefile} {path1}/{path2}')

def get_config(config_file = 'config.ini'):
    config = configparser.ConfigParser(delimiters=' ') 
    config.read(config_file)
    return(config)
    
def get_filestructure(outfile):
    file_structure = outfile.split('/')
    filename = file_structure[-1]
    path = ('/').join(file_structure[:-1])
    if path is '': path = '.'
    return(path, filename)

"""
if __name__ == "__main__":
    data = get_geodesic('out2.bin',
                        [91.08497395, 0., 1.04719755, 0., 0.],
                        return_data = True)
"""


