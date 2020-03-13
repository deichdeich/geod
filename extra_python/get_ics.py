import numpy as np
import pth_kso2
import pth_edgbo2
import pth_edgbo1SR

pth_dict = {'kerr_so2': pth_kso2.pth_sol,
            'edgb_o2': pth_edgbo2.pth_sol,
            'edgb_o1SR': pth_edgbo1SR.pth_sol}

def make_ic_range(theory, nstates, rmin, rmax, **kwargs):
    func = pth_dict[theory]
    rs = np.linspace(rmin,rmax,nstates)
    data = func(r = rs, **kwargs)
    ics = np.zeros((nstates,6))
    ics[:,0] = rs
    ics[:,1] = 0
    ics[:,2] = np.pi/2
    ics[:,3] = data
    ics[:,4] = 0
    ics[:,5] = kwargs['L']
    return(ics)

