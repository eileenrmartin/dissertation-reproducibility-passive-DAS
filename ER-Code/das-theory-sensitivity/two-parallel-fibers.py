from VelAndStrRate import *
import numpy as np
import matplotlib.pyplot as plt
import sys
fig = sys.argv[1]

# cross correlate two records 
def xcorrelation(record1, record2, nLagsPerSide):
	nsteps = record2.size
	record2Padded = np.zeros(nsteps+nLagsPerSide*2,dtype=np.float32)
	record2Padded[nLagsPerSide:-nLagsPerSide] = record2
	xcorr = np.zeros(1+nLagsPerSide*2,dtype=np.float32)
	for i in range(1+nLagsPerSide*2):
		xcorr[i] = np.dot(record1,record2Padded[i:nsteps+i])
	return xcorr

c = 2000.0 # velocity m/s
f = 30.0 # frequency Hz
minRad = 2000.0 # inner radius of source annulus m
maxRad = 3000.0 # outer radius of source annulus m
xAs = [(-425.0,-850.0),(-425.0,-637.5),(-425.0,-425.0),(-425.0,-212.5),(-425.0,0.0),(-425.0,212.5),(-425.0,425.0)] # coordinates of receiver A
xB = (425.0, 425.0) # coordinates of receiver B

# for stack over many cross-correlations
cL = 2400 # m/s Love wave velocity
cR = 2000 # m/s Rayleigh wave velocity
nSrc = 2000 # number of sources that will be recorded by A and B in separate records
dt = 0.005 # time step s
scalePlot = 20000
timesLong = np.arange(0.0,400000.0, dt) # times over which A and B record responses to random sources in one time series
ntlong = timesLong.size
nSrcLong = 10000
times = np.arange(0.0,2.0, dt) # times at which receiver A and B record responses to each source
nt = times.size # number of time steps recorded
maxLag = 1.0 # maximum lag for cross-correlation s
nLags = int(maxLag/dt) # number of time lag steps to do cross-correlation over
timesXcorrs = np.arange(-nLags*dt,(nLags+1)*dt,dt) # time lags at which cross correlation is calculated




for ixA,xA in enumerate(xAs):
    longRecAFib = np.zeros(ntlong,dtype=np.float32)
    longRecBFib = np.zeros(ntlong,dtype=np.float32)
    for i in range(nSrcLong):
        if i%4000 == 0:
            print(i)
        np.random.seed(i+2*nSrcLong*ixA)
        # set off a random Love wave source
        rSrcL = np.random.uniform(minRad,maxRad)
        phiSrcL = np.random.uniform(-pi/2,3*pi/2)
        xsrcL = (-np.cos(phiSrcL)*rSrcL,-np.sin(phiSrcL)*rSrcL)
        startIdxL = np.random.randint(0,ntlong-nt)
        longRecAFib[startIdxL:nt+startIdxL] = longRecAFib[startIdxL:nt+startIdxL] + sigjj_Love(xsrcL, times, f, cL, xA)
        longRecBFib[startIdxL:nt+startIdxL] = longRecBFib[startIdxL:nt+startIdxL] + sigjj_Love(xsrcL, times, f, cL, xB)
        # set off a random Rayleigh wave source
        rSrcR = np.random.uniform(minRad,maxRad)
        phiSrcR = np.random.uniform(-pi/2,3*pi/2)
        xsrcR = (-np.cos(phiSrcR)*rSrcR,-np.sin(phiSrcR)*rSrcR)
        startIdxR = np.random.randint(0,ntlong-nt)
        longRecAFib[startIdxR:nt+startIdxR] = longRecAFib[startIdxR:nt+startIdxR] + sigjj_Rayleigh(xsrcR, times, f, cR, xA)
        longRecBFib[startIdxR:nt+startIdxR] = longRecBFib[startIdxR:nt+startIdxR] + sigjj_Rayleigh(xsrcR, times, f, cR, xB)
   
    xcorrsLongFib = xcorrelation(longRecAFib, longRecBFib, nLags)
    colorFrac = float(ixA)/float(len(xAs))
    plt.plot(timesXcorrs,xA[1]+scalePlot*xcorrsLongFib,linewidth=0.5,c=(1.0-colorFrac,0,colorFrac))
    plt.ylim(-1200, 600)
    # put a marker at the time lags indicating true Rayleigh and Love wave velocities
    dist = np.sqrt((xA[1]-xB[1])**2+(xA[0]-xB[0])**2)
    RayleighTime = dist/cR
    LoveTime = dist/cL
    plt.scatter(RayleighTime,xA[1],c='y')
    plt.scatter(LoveTime,xA[1],c='k')
    plt.scatter(-RayleighTime,xA[1],c='y')
    plt.scatter(-LoveTime,xA[1],c='k')
    plt.title('single long cross-correlation, Rayleigh and Love sources')
plt.xlabel('time-lag (s)')
plt.ylabel('j coordinate of receiver')
plt.savefig(fig+'two-parallel-lines')