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
xA = (-600.0,0.0) # coordinates of receiver A
xB = (600.0, 0.0) # coordinates of receiver B

# for stack over many cross-correlations
nSrc = 5000 # number of sources that will be recorded by A and B in separate records
dt = 0.005 # time step s
times = np.arange(0.0,2.0, dt) # times at which receiver A and B record responses to each source
nt = times.size # number of time steps recorded
maxLag = 1.0 # maximum lag for cross-correlation s
nLags = int(maxLag/dt) # number of time lag steps to do cross-correlation over
timesXcorrs = np.arange(-nLags*dt,(nLags+1)*dt,dt) # time lags at which cross correlation is calculated

# for single cross-correlation
timesLong = np.arange(0.0,40000.0, dt) # times over which A and B record responses to random sources in one time series
ntlong = timesLong.size
nSrcLong = 1000 # number of sources randomly distributed throughout timesLong recorded by A and B

recordsAGeoR = np.zeros((nSrc,nt),dtype = np.float32)
recordsBGeoR = np.zeros((nSrc,nt),dtype = np.float32)
xcorrsGeoR = np.zeros((nSrc,1+2*nLags),dtype = np.float32)
recordsAFibR = np.zeros((nSrc,nt),dtype = np.float32)
recordsBFibR = np.zeros((nSrc,nt),dtype = np.float32)
xcorrsFibR = np.zeros((nSrc,1+2*nLags),dtype = np.float32)
phi = np.zeros(nSrc,dtype=np.float32)



plotFrac = 4 # only one of every plotFrac sources is plotted

# plot geometry
plt.clf()
for i in range(nSrc):
	if i % plotFrac == 1:
		# choose source location randomly
		np.random.seed(i)
		rSrc = np.random.uniform(minRad,maxRad)
		phiSrc = np.random.uniform(-pi/2,3*pi/2) 
		phi[i] = phiSrc
		xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)

		# plot source
		plt.scatter(xsrc[0],xsrc[1],marker='.',c='k',s=1)

plt.scatter(xA[0],xA[1],marker='_',c='b',s=40)
plt.scatter(xB[0],xB[1],marker='_',c='b',s=40)
plt.annotate('x1',(xA[0]-300,xA[1]+300))
plt.annotate('x2',(xB[0]-300,xB[1]+300))
plt.annotate(r'$\phi_{Src}$'+' \n '+r'$=$'+' \n '+r' $0$',(-maxRad-800,0))
plt.annotate(r'$\phi_{Src} = \pi/2$',(-500,-maxRad-300))
plt.annotate(r'$\phi_{Src} = 3\pi/2$',(-500,maxRad+300))
plt.annotate(r'$\phi_{Src}$'+' \n '+r'$=$'+' \n '+r' $\pi$',(maxRad+100,0))
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
plt.savefig(fig+'RRXcorrGeom')




scaleGeo = 0.6 # scale factor just to make plots look nice
scaleFib = 0.03 # scale factor just to make plots look nice
for i in range(nSrc):
	# choose source location randomly
	np.random.seed(i)
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2) 
	phi[i] = phiSrc
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)

	# geophone records and cross-correlation
	recordsAGeoR[i,:] = vi_Rayleigh(xsrc, times, f, c, xA)
	recordsAGeoR[i,:] = recordsAGeoR[i,:] / scaleGeo
	recordsBGeoR[i,:] = vi_Rayleigh(xsrc, times, f, c, xB)
	recordsBGeoR[i,:] = recordsBGeoR[i,:] / scaleGeo
	xcorrsGeoR[i,:] = xcorrelation(recordsAGeoR[i,:], recordsBGeoR[i,:], nLags)

	# fiber records and cross-correlation
	recordsAFibR[i,:] = sigii_Rayleigh(xsrc, times, f, c, xA)
	recordsAFibR[i,:] = recordsAFibR[i,:] / scaleFib
	recordsBFibR[i,:] = sigii_Rayleigh(xsrc, times, f, c, xB)
	recordsBFibR[i,:] = recordsBFibR[i,:] / scaleFib
	xcorrsFibR[i,:] = xcorrelation(recordsAFibR[i,:], recordsBFibR[i,:], nLags)
    
# plot the records (or some portion of them)
plt.clf()
plotFrac = 8 # only one of every plotFrac sources is plotted
plt.subplot(1, 2, 1)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsAGeoR[i,:],c='k',linewidth=0.2)
plt.title('records geophone 1')
plt.ylabel(r'$\phi_{Src}$',fontsize = 26)
plt.xlabel('t (s)')
plt.subplot(1, 2, 2)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsBGeoR[i,:],c='k',linewidth=0.2)
plt.title('records geophone 2')
plt.xlabel('t (s)')
plt.savefig(fig+'geoRecsRRRayleigh')

plt.clf()

# plot the records (or some portion of them)
plt.subplot(1, 2, 1)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsAFibR[i,:],c='r',linewidth=0.2)
plt.title('records fiber 1')
plt.ylabel(r'$\phi_{Src}$',fontsize=26)
plt.xlabel('t (s)')
plt.subplot(1, 2, 2)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsBFibR[i,:],c='r',linewidth=0.2)
plt.title('records B fiber')
plt.ylabel('phiSrc')
plt.xlabel('t (s)')
plt.savefig(fig+'fiberRecsRRRayelgih')


plt.clf()
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(timesXcorrs,phi[i]+xcorrsFibR[i,:],c='r',linewidth=0.3)
		plt.plot(timesXcorrs,phi[i]+xcorrsGeoR[i,:],c='k',linewidth=0.2)
plt.title('source-wise \n cross-correlations')
plt.xlim(-maxLag,maxLag)
plt.ylabel(r'$\phi_{Src}$',fontsize=26)
plt.xlabel('t (s)')
plt.savefig(fig+'XCorrsRRRayleigh')

plt.clf()
xcorrStackFibR = np.sum(xcorrsFibR,axis=0)
plt.plot(timesXcorrs,xcorrStackFibR,c='r',label="fiber")
xcorrStackGeoR = np.sum(xcorrsGeoR,axis=0)
plt.plot(timesXcorrs,xcorrStackGeoR,c='k',label="geophone")
plt.xlim(-maxLag,maxLag)
plt.title('avg. single-source cross-correlation')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'stackedXCorrsRRRayleigh')

# long record averages
longRecAGeoR = np.zeros(ntlong,dtype=np.float32)
longRecBGeoR = np.zeros(ntlong,dtype=np.float32)
longRecAFibR = np.zeros(ntlong,dtype=np.float32)
longRecBFibR = np.zeros(ntlong,dtype=np.float32)
for i in range(nSrcLong):
	if i%1000 == 0:
		print(i)
	np.random.seed(i)
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoR[startIdx:nt+startIdx] = longRecAGeoR[startIdx:nt+startIdx] + vi_Rayleigh(xsrc, times, f, c, xA)
	longRecBGeoR[startIdx:nt+startIdx] = longRecBGeoR[startIdx:nt+startIdx] + vi_Rayleigh(xsrc, times, f, c, xB)
	longRecAFibR[startIdx:nt+startIdx] = longRecAFibR[startIdx:nt+startIdx] + sigii_Rayleigh(xsrc, times, f, c, xA)
	longRecBFibR[startIdx:nt+startIdx] = longRecBFibR[startIdx:nt+startIdx] + sigii_Rayleigh(xsrc, times, f, c, xB)

plt.clf()
xcorrsLongGeoR = xcorrelation(longRecAGeoR/scaleGeo, longRecBGeoR/scaleGeo, nLags)
plt.plot(timesXcorrs,xcorrsLongGeoR,c='k',label="geophone")
xcorrsLongFibR = xcorrelation(longRecAFibR/scaleFib, longRecBFibR/scaleFib, nLags)
plt.plot(timesXcorrs,xcorrsLongFibR,c='r',label="fiber")
plt.xlim(-maxLag,maxLag)
plt.title('long RR cross-correlation (Rayleigh sources only)')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'longXcorrsRRRayleigh')




# long record averages rayleigh and love sources
longRecAGeoR = np.zeros(ntlong,dtype=np.float32)
longRecBGeoR = np.zeros(ntlong,dtype=np.float32)
longRecAFibR = np.zeros(ntlong,dtype=np.float32)
longRecBFibR = np.zeros(ntlong,dtype=np.float32)
for i in range(nSrcLong):
	if i%1000 == 0:
		print(i)
	np.random.seed(i)
	# Rayleigh wave source
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoR[startIdx:nt+startIdx] = longRecAGeoR[startIdx:nt+startIdx] + vi_Rayleigh(xsrc, times, f, c, xA)
	longRecBGeoR[startIdx:nt+startIdx] = longRecBGeoR[startIdx:nt+startIdx] + vi_Rayleigh(xsrc, times, f, c, xB)
	longRecAFibR[startIdx:nt+startIdx] = longRecAFibR[startIdx:nt+startIdx] + sigii_Rayleigh(xsrc, times, f, c, xA)
	longRecBFibR[startIdx:nt+startIdx] = longRecBFibR[startIdx:nt+startIdx] + sigii_Rayleigh(xsrc, times, f, c, xB)
	# Love wave source
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoR[startIdx:nt+startIdx] = longRecAGeoR[startIdx:nt+startIdx] + vi_Love(xsrc, times, f, c, xA)
	longRecBGeoR[startIdx:nt+startIdx] = longRecBGeoR[startIdx:nt+startIdx] + vi_Love(xsrc, times, f, c, xB)
	longRecAFibR[startIdx:nt+startIdx] = longRecAFibR[startIdx:nt+startIdx] + sigii_Love(xsrc, times, f, c, xA)
	longRecBFibR[startIdx:nt+startIdx] = longRecBFibR[startIdx:nt+startIdx] + sigii_Love(xsrc, times, f, c, xB)

plt.clf()
xcorrsLongGeoR = xcorrelation(longRecAGeoR/scaleGeo, longRecBGeoR/scaleGeo, nLags)
plt.plot(timesXcorrs,xcorrsLongGeoR,c='k',label="geophone")
xcorrsLongFibR = xcorrelation(longRecAFibR/scaleFib, longRecBFibR/scaleFib, nLags)
plt.plot(timesXcorrs,xcorrsLongFibR,c='r',label="fiber")
plt.xlim(-maxLag,maxLag)
plt.title('long RR cross-correlation (Love and Rayleigh sources)')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'longXcorrsRRBoth')