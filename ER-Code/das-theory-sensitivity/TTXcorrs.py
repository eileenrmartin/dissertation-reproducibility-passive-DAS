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

recordsAGeoT = np.zeros((nSrc,nt),dtype = np.float32)
recordsBGeoT = np.zeros((nSrc,nt),dtype = np.float32)
xcorrsGeoT = np.zeros((nSrc,1+2*nLags),dtype = np.float32)
recordsAFibT = np.zeros((nSrc,nt),dtype = np.float32)
recordsBFibT = np.zeros((nSrc,nt),dtype = np.float32)
xcorrsFibT = np.zeros((nSrc,1+2*nLags),dtype = np.float32)
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

plt.scatter(xA[0],xA[1],marker='|',c='m',s=40)
plt.scatter(xB[0],xB[1],marker='|',c='m',s=40)
plt.annotate('x1',(xA[0]-300,xA[1]+300))
plt.annotate('x2',(xB[0]-300,xB[1]+300))
plt.annotate(r'$\phi_{Src}$'+' \n '+r'$=$'+' \n '+r' $0$',(-maxRad-800,0))
plt.annotate(r'$\phi_{Src} = \pi/2$',(-500,-maxRad-300))
plt.annotate(r'$\phi_{Src} = 3\pi/2$',(-500,maxRad+300))
plt.annotate(r'$\phi_{Src}$'+' \n '+r'$=$'+' \n '+r' $\pi$',(maxRad+100,0))
plt.xlabel('x')
plt.ylabel('y')
plt.axis('scaled')
plt.savefig(fig+'TTXcorrGeom')


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
	recordsAGeoT[i,:] = vj_Love(xsrc, times, f, c, xA)
	recordsAGeoT[i,:] = recordsAGeoT[i,:] / scaleGeo
	recordsBGeoT[i,:] = vj_Love(xsrc, times, f, c, xB)
	recordsBGeoT[i,:] = recordsBGeoT[i,:] / scaleGeo
	xcorrsGeoT[i,:] = xcorrelation(recordsAGeoT[i,:], recordsBGeoT[i,:], nLags)

	# fiber records and cross-correlation
	recordsAFibT[i,:] = sigjj_Love(xsrc, times, f, c, xA)
	recordsAFibT[i,:] = recordsAFibT[i,:] / scaleFib
	recordsBFibT[i,:] = sigjj_Love(xsrc, times, f, c, xB)
	recordsBFibT[i,:] = recordsBFibT[i,:] / scaleFib
	xcorrsFibT[i,:] = xcorrelation(recordsAFibT[i,:], recordsBFibT[i,:], nLags)
    
# plot the records (or some portion of them)
plotFrac = 8
plt.subplot(1, 2, 1)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsAGeoT[i,:],c='k',linewidth=0.2)
plt.title('records geophone 1')
plt.ylabel(r'$\phi_{Src}$',fontsize = 26)
plt.xlabel('t (s)')
plt.subplot(1, 2, 2)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsBGeoT[i,:],c='k',linewidth=0.2)
plt.title('records geophone 2')
plt.xlabel('t (s)')
plt.savefig(fig+'geoRecsTTLove')

plt.clf()

# plot the records (or some portion of them)
plt.subplot(1, 2, 1)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsAFibT[i,:],c='r',linewidth=0.2)
plt.title('records fiber 1')
plt.ylabel(r'$\phi_{Src}$',fontsize=26)
plt.xlabel('t (s)')
plt.subplot(1, 2, 2)
for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(times,phi[i]+recordsBFibT[i,:],c='r',linewidth=0.2)
plt.title('records fiber 2')
plt.xlabel('t (s)')
plt.savefig(fig+'fiberRecsTTLove')


plt.clf()

for i in range(nSrc):
	if i%plotFrac == 1:
		plt.plot(timesXcorrs,phi[i]+xcorrsFibT[i,:],c='r',linewidth=0.3)
		plt.plot(timesXcorrs,phi[i]+xcorrsGeoT[i,:],c='k',linewidth=0.2)
plt.title('source-wise \n cross-correlations')
plt.xlim(-maxLag,maxLag)
plt.ylabel(r'$\phi_{Src}$',fontsize=26)
plt.xlabel('t (s)')
plt.savefig(fig+'XCorrsTTLove')

plt.clf()
xcorrStackFibT = np.sum(xcorrsFibT,axis=0)
plt.plot(timesXcorrs,xcorrStackFibT,c='r',label="fiber")
xcorrStackGeoT = np.sum(xcorrsGeoT,axis=0)
plt.plot(timesXcorrs,xcorrStackGeoT,c='k',label="geophone")
plt.xlim(-maxLag,maxLag)
#plt.legend(bbox_to_anchor=(.85, .5), loc=0, borderaxespad=0.5,prop={'size': 10})
plt.title('avg. single-source cross-correlation')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'stackedXCorrsTTLove')

# long record averages
longRecAGeoT = np.zeros(ntlong,dtype=np.float32)
longRecBGeoT = np.zeros(ntlong,dtype=np.float32)
longRecAFibT = np.zeros(ntlong,dtype=np.float32)
longRecBFibT = np.zeros(ntlong,dtype=np.float32)
for i in range(nSrcLong):
	if i%1000 == 0:
		print(i)
	np.random.seed(i)
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoT[startIdx:nt+startIdx] = longRecAGeoT[startIdx:nt+startIdx] + vj_Love(xsrc, times, f, c, xA)
	longRecBGeoT[startIdx:nt+startIdx] = longRecBGeoT[startIdx:nt+startIdx] + vj_Love(xsrc, times, f, c, xB)
	longRecAFibT[startIdx:nt+startIdx] = longRecAFibT[startIdx:nt+startIdx] + sigjj_Love(xsrc, times, f, c, xA)
	longRecBFibT[startIdx:nt+startIdx] = longRecBFibT[startIdx:nt+startIdx] + sigjj_Love(xsrc, times, f, c, xB)

plt.clf()
xcorrsLongGeoT = xcorrelation(longRecAGeoT/scaleGeo, longRecBGeoT/scaleGeo, nLags)
plt.plot(timesXcorrs,xcorrsLongGeoT,c='k',label="geophone")
xcorrsLongFibT = xcorrelation(longRecAFibT/scaleFib, longRecBFibT/scaleFib, nLags)
plt.plot(timesXcorrs,xcorrsLongFibT,c='r',label="fiber")
plt.xlim(-maxLag,maxLag)
#plt.legend(bbox_to_anchor=(.85, .5), loc=0, borderaxespad=0.5,prop={'size': 10})
plt.title('long TT cross-correlation (Love sources only)')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'longXcorrsTTLove')




# long record averages rayleigh and love sources
longRecAGeoT = np.zeros(ntlong,dtype=np.float32)
longRecBGeoT = np.zeros(ntlong,dtype=np.float32)
longRecAFibT = np.zeros(ntlong,dtype=np.float32)
longRecBFibT = np.zeros(ntlong,dtype=np.float32)
for i in range(nSrcLong):
	if i%1000 == 0:
		print(i)
	np.random.seed(i)
	# Rayleigh wave source
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoT[startIdx:nt+startIdx] = longRecAGeoT[startIdx:nt+startIdx] + vj_Rayleigh(xsrc, times, f, c, xA)
	longRecBGeoT[startIdx:nt+startIdx] = longRecBGeoT[startIdx:nt+startIdx] + vj_Rayleigh(xsrc, times, f, c, xB)
	longRecAFibT[startIdx:nt+startIdx] = longRecAFibT[startIdx:nt+startIdx] + sigjj_Rayleigh(xsrc, times, f, c, xA)
	longRecBFibT[startIdx:nt+startIdx] = longRecBFibT[startIdx:nt+startIdx] + sigjj_Rayleigh(xsrc, times, f, c, xB)
	# Love wave source
	rSrc = np.random.uniform(minRad,maxRad)
	phiSrc = np.random.uniform(-pi/2,3*pi/2)
	xsrc = (-np.cos(phiSrc)*rSrc,-np.sin(phiSrc)*rSrc)
	startIdx = np.random.randint(0,ntlong-nt)
	longRecAGeoT[startIdx:nt+startIdx] = longRecAGeoT[startIdx:nt+startIdx] + vj_Love(xsrc, times, f, c, xA)
	longRecBGeoT[startIdx:nt+startIdx] = longRecBGeoT[startIdx:nt+startIdx] + vj_Love(xsrc, times, f, c, xB)
	longRecAFibT[startIdx:nt+startIdx] = longRecAFibT[startIdx:nt+startIdx] + sigjj_Love(xsrc, times, f, c, xA)
	longRecBFibT[startIdx:nt+startIdx] = longRecBFibT[startIdx:nt+startIdx] + sigjj_Love(xsrc, times, f, c, xB)

plt.clf()
xcorrsLongGeoT = xcorrelation(longRecAGeoT/scaleGeo, longRecBGeoT/scaleGeo, nLags)
plt.plot(timesXcorrs,xcorrsLongGeoT,c='k',label="geophone")
xcorrsLongFibT = xcorrelation(longRecAFibT/scaleFib, longRecBFibT/scaleFib, nLags)
plt.plot(timesXcorrs,xcorrsLongFibT,c='r',label="fiber")
plt.xlim(-maxLag,maxLag)
#plt.legend(bbox_to_anchor=(.85, .5), loc=0, borderaxespad=0.5,prop={'size': 10})
plt.title('long TT cross-correlation (Love and Rayleigh sources)')
plt.xlabel('time-lag (s)')
plt.savefig(fig+'longXcorrsTTBoth')