import numpy as np

def get_bound_y0_range(th, Ee, Jz, a, M, Ny0, radius_bounds = None):
    if radius_bounds == None:
        radius_bounds = get_czv_radius_bounds(th, Ee, a, Jz, M, -1)
    rs = np.linspace(radius_bounds[0], radius_bounds[1], Ny0)
    y0s = np.zeros((Ny0, 6))
    for i in range(Ny0):
        y0s[i] = np.array([rs[i], 0, th, 0, 0, Jz], dtype = np.float64)
    return(y0s)

def get_bound_y0(th, Ee, Jz, a, M):
    radius_bounds = get_czv_radius_bounds(th, Ee, a, Jz, M, -1)
    r = np.random.random() * (radius_bounds[1] - radius_bounds[0]) + radius_bounds[0]
    y0 = np.array([r, 0, th, 0, 0, Jz], dtype=np.float64)
    return(y0)
    
def get_czv_radius_bounds(th, e, a, j, M, mu):
    isco = get_isco(a, M)
    x = np.arange(isco, isco*100, 0.001)
    bound_orbits = np.where(veff_sr_o2(x,th, e, a, j, M, mu)<0)[0]
    if len(bound_orbits) == 0:
        raise ValueError('This region of parameter space contains no bound orbit')
    return(x[bound_orbits[0]],x[bound_orbits[-1]])
    
def get_plunge_y0_range(th, Ee, Jz, a, M, Ny0):
    isco = get_isco(a, M)
    schw = 2 * M
    rs = np.linspace(schw, isco, Ny0)
    y0s = np.zeros((Ny0, 6))
    for i in range(Ny0):
        y0s[i] = np.array([rs[i], 0, th, 0, 0, Jz], dtype = np.float64)
    return(sch, isco)
    
def get_isco(a, M):
    Z1 = 1 + (1-a**2)**(1/3) * ((1+a)**(1/3) + (1-a)**(1/3))
    Z2 = np.sqrt((3 * a**2) + Z1**2)
    r_isco = M*(3 + Z2 - np.sqrt((3-Z1)*(3+Z1+(2*Z2))))
    return(r_isco)

def func(r, Ee, a, Jz, M, mu):
    th = cmath.atan(-cmath.sqrt(-((-(pow(a,4)*pow(Ee,2)) + pow(a,2)*pow(Jz,2) + 
          pow(a,4)*pow(mu,2) + 4*pow(a,2)*pow(Ee,2)*M*r - 4*a*Ee*Jz*M*r - 
          2*pow(a,2)*M*pow(mu,2)*r + 2*M*pow(mu,2)*pow(r,3) + 
          pow(Ee,2)*pow(r,4) - pow(mu,2)*pow(r,4) + 
          cmath.sqrt(-4*pow(a,2)*pow(Jz,2)*(pow(Ee,2) - pow(mu,2))*
             pow(pow(a,2) + r*(-2*M + r),2) + 
            pow(pow(a,4)*(pow(Ee,2) - pow(mu,2)) - 4*a*Ee*Jz*M*r + 
              pow(r,3)*(2*M*pow(mu,2) + (pow(Ee,2) - pow(mu,2))*r) + 
              pow(a,2)*(pow(Jz,2) + 
                 2*r*(M*pow(mu,2) + pow(Ee,2)*r - pow(mu,2)*r)),2)))/
        (pow(a,2)*(pow(Ee,2) - pow(mu,2))*(pow(a,2) + r*(-2*M + r))))),
   -cmath.sqrt((pow(a,4)*pow(Ee,2) + pow(a,2)*pow(Jz,2) - pow(a,4)*pow(mu,2) - 
        4*a*Ee*Jz*M*r + 2*pow(a,2)*M*pow(mu,2)*r + 
        2*pow(a,2)*pow(Ee,2)*pow(r,2) - 2*pow(a,2)*pow(mu,2)*pow(r,2) + 
        2*M*pow(mu,2)*pow(r,3) + pow(Ee,2)*pow(r,4) - pow(mu,2)*pow(r,4) + 
        cmath.sqrt(-4*pow(a,2)*pow(Jz,2)*(pow(Ee,2) - pow(mu,2))*
           pow(pow(a,2) + r*(-2*M + r),2) + 
          pow(pow(a,4)*(pow(Ee,2) - pow(mu,2)) - 4*a*Ee*Jz*M*r + 
            pow(r,3)*(2*M*pow(mu,2) + (pow(Ee,2) - pow(mu,2))*r) + 
            pow(a,2)*(pow(Jz,2) + 
               2*r*(M*pow(mu,2) + pow(Ee,2)*r - pow(mu,2)*r)),2)))/
      (pow(a,2)*(pow(Ee,2) - pow(mu,2))*(pow(a,2) + r*(-2*M + r)))))
    return(th)

def veff_exact(r, th, e, a, j, M, mu):
    retval = (pow(mu,2) + (pow(j,2)*(-1 + (2*M*r)/(pow(r,2) + pow(a,2)*pow(np.cos(th),2))) - 
        (4*a*e*j*M*r*pow(np.sin(th),2))/(pow(r,2) + pow(a,2)*pow(np.cos(th),2)) + 
        (pow(e,2)*pow(np.sin(th),2)*
           (pow(pow(a,2) + pow(r,2),2) - 
             pow(a,2)*(pow(a,2) + r*(-2*M + r))*pow(np.sin(th),2)))/
         (pow(r,2) + pow(a,2)*pow(np.cos(th),2)))/
      ((-4*pow(a,2)*pow(M,2)*pow(r,2)*pow(np.sin(th),4))/
         pow(pow(r,2) + pow(a,2)*pow(np.cos(th),2),2) + 
        ((-1 + (2*M*r)/(pow(r,2) + pow(a,2)*pow(np.cos(th),2)))*pow(np.sin(th),2)*
           (pow(pow(a,2) + pow(r,2),2) - 
             pow(a,2)*(pow(a,2) + r*(-2*M + r))*pow(np.sin(th),2)))/
         (pow(r,2) + pow(a,2)*pow(np.cos(th),2))))/2.
    return(retval)

def veff_sr_o2(r, th, e, a, j, M, mu):
    retval = (pow(mu,2) + (pow(j,2)*(-1 + (2*M)/r - 
           (2*pow(a,2)*M*pow(np.cos(th),2))/pow(r,3)) - 
        (4*a*e*j*M*pow(np.sin(th),2))/r + 
        pow(e,2)*(pow(r,2)*pow(np.sin(th),2) + 
           (pow(a,2)*(2*pow(M,2)*pow(np.sin(th),2) - 
                pow(M,2)*pow(np.cos(th),2)*pow(np.sin(th),2) - 
                pow(M,2)*pow(np.sin(th),4) + (2*pow(M,3)*pow(np.sin(th),4))/r))/
            pow(M,2)))/
      ((-4*pow(a,2)*pow(M,2)*pow(np.sin(th),4))/pow(r,2) + 
        (-1 + (2*M)/r - (2*pow(a,2)*M*pow(np.cos(th),2))/pow(r,3))*
         (pow(r,2)*pow(np.sin(th),2) + 
           (pow(a,2)*(2*pow(M,2)*pow(np.sin(th),2) - 
                pow(M,2)*pow(np.cos(th),2)*pow(np.sin(th),2) - 
                pow(M,2)*pow(np.sin(th),4) + (2*pow(M,3)*pow(np.sin(th),4))/r))/
            pow(M,2))))/2.
    return(retval)
