import numpy as np
import matplotlib.pyplot as plt
pi = np.pi

def vi_Love(xs, t, f, c, x):
    xj = x[1]
    xsj = xs[1]
    R = np.sqrt((x[0]-xs[0])**2+(xj-xsj)**2)
    termA = (xj-xsj)/R**2
    termB = np.exp(-(pi*f*(t-R/c))**2)
    termC = -6*pi*pi*f*f*(t-R/c)
    termD = 4*(pi**4)*(f**4)*(t-R/c)**3
    v = termA*termB*(termC+termD)
    return v

def vj_Love(xs, t, f, c, x):
    xi = x[0]
    xsi = xs[0]
    R = np.sqrt((xi-xsi)**2+(x[1]-xs[1])**2)
    termA = (xsi-xi)/R**2
    termB = np.exp(-(pi*f*(t-R/c))**2)
    termC = -6*pi*pi*f*f*(t-R/c)
    termD = 4*(pi**4)*(f**4)*(t-R/c)**3
    v = termA*termB*(termC+termD)
    return v

def sigii_Love(xs, t, f, c, x):
    xi = x[0]
    xj = x[1]
    xsi = xs[0]
    xsj = xs[1]
    R = np.sqrt((xi-xsi)**2+(xj-xsj)**2)
    termA = (xi-xsi)*(xj-xsj)/(R**3)
    tau = t-R/c
    termB = np.exp(-(pi*f*tau)**2)
    termC = (6*pi*pi*f*f/c) + (12*pi*pi*f*f*tau/R) - (24*(pi*f)**4*(tau**2)/c) - (8*(pi*f)**4*(tau**3)/R) + (8*(pi*f)**6*(tau**4)/c)
    sig = termA*termB*termC
    return sig

def sigjj_Love(xs, t, f, c, x):
    xi = x[0]
    xj = x[1]
    xsi = xs[0]
    xsj = xs[1]
    R = np.sqrt((xi-xsi)**2+(xj-xsj)**2)
    termA = (xsi-xi)*(xj-xsj)/(R**3)
    tau = t-R/c
    termB = np.exp(-(pi*f*tau)**2)
    termC = (6*pi*pi*f*f/c) + (12*pi*pi*f*f*tau/R) - (24*(pi*f)**4*(tau**2)/c) - (8*(pi*f)**4*(tau**3)/R) + (8*(pi*f)**6*(tau**4)/c)
    sig = termA*termB*termC
    return sig

def vi_Rayleigh(xs, t, f, c, x):
    xi = x[0]
    xsi = xs[0]
    R = np.sqrt((xi-xsi)**2+(x[1]-xs[1])**2)
    termA = (xi-xsi)/R**2
    termB = np.exp(-(pi*f*(t-R/c))**2)
    termC = -6*pi*pi*f*f*(t-R/c)
    termD = 4*(pi**4)*(f**4)*(t-R/c)**3
    v = termA*termB*(termC+termD)
    return v

def vj_Rayleigh(xs, t, f, c, x):
    xj = x[1]
    xsj = xs[1]
    R = np.sqrt((x[0]-xs[0])**2+(xj-xsj)**2)
    termA = (xj-xsj)/R**2
    termB = np.exp(-(pi*f*(t-R/c))**2)
    termC = -6*pi*pi*f*f*(t-R/c)
    termD = 4*(pi**4)*(f**4)*(t-R/c)**3
    v = termA*termB*(termC+termD)
    return v

def sigii_Rayleigh(xs, t, f, c, x):
    xi = x[0]
    xsi = xs[0]
    R = np.sqrt((xi-xsi)**2+(x[1]-xs[1])**2)
    termA = (xi-xsi)**2/(R**3)
    tau = t-R/c
    termB = np.exp(-(pi*f*tau)**2)
    termC = (6*pi*pi*f*f/c) + (12*pi*pi*f*f*tau/R) - (24*(pi*f)**4*(tau**2)/c) - (8*(pi*f)**4*(tau**3)/R) + (8*(pi*f)**6*(tau**4)/c)
    sig = termA*termB*termC
    return sig

def sigjj_Rayleigh(xs, t, f, c, x):
    xj = x[1]
    xsj = xs[1]
    R = np.sqrt((x[0]-xs[0])**2+(xj-xsj)**2)
    termA = (xj-xsj)**2/(R**3)
    tau = t-R/c
    termB = np.exp(-(pi*f*tau)**2)
    termC = (6*pi*pi*f*f/c) + (12*pi*pi*f*f*tau/R) - (24*(pi*f)**4*(tau**2)/c) - (8*(pi*f)**4*(tau**3)/R) + (8*(pi*f)**6*(tau**4)/c)
    sig = termA*termB*termC
    return sig
    