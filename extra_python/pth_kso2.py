import numpy as np

def csc(x):
    return(1/np.sin(x))

def cot(x):
    return(1/np.tan(x))

def pth_sol(M, a, Ee, L, Th, r):
    pth = (np.sqrt(1/(pow(a,2) + 2*pow(r,2) + pow(a,2)*np.cos(2*Th)))*np.sqrt(pow(a,2) + 2*pow(r,2) + pow(a,2)*np.cos(2*Th))*np.sqrt((pow(a,2) + 2*pow(r,2) + pow(a,2)*np.cos(2*Th))*(-16*a*Ee*L*M*pow(r,3) - 2*pow(a,4)*M*(M + r) + pow(a,2)*(-8*pow(L,2)*M*r + 6*pow(Ee,2)*M*pow(r,3) + 4*(-1 + pow(Ee,2))*pow(r,4)) + 4*pow(r,3)*(pow(L,2)*(4*M - 2*r) + pow(r,2)*(2*M + (-1 + pow(Ee,2))*r)) + (pow(a,4)*pow(M,2) + 16*a*Ee*L*M*pow(r,3) + 4*pow(r,5)*(-2*M + r - pow(Ee,2)*r) - 4*pow(a,2)*(2*pow(L,2)*M*r + pow(r,3)*(-r + pow(Ee,2)*(2*M + r))))*np.cos(2*Th) + 2*pow(a,2)*M*(pow(Ee,2)*pow(r,3) + pow(a,2)*(M + r))*np.cos(4*Th) - pow(a,4)*pow(M,2)*np.cos(6*Th))*pow(csc(Th),2)))/(2.*np.sqrt(4*pow(a,2)*pow(r,4) + 4*pow(r,5)*(-2*M + r) + 2*pow(a,4)*M*(M + 2*r) - 2*pow(a,4)*M*(-2*r*np.cos(2*Th) + M*np.cos(4*Th))))
    return(pth)
