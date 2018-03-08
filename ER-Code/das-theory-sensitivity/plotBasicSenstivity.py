import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
import sys
fig = sys.argv[1]

def LoveSensitivityGeo(c,lam,theta,phis):
    return 2*pi*c*np.absolute(-np.cos(theta)*np.sin(phis)+np.sin(theta)*np.cos(phis))/lam
def RayleighSensitivityGeo(c,lam,theta,phis):
    return 2*pi*c*np.absolute(np.cos(theta)*np.cos(phis)+np.sin(theta)*np.sin(phis))/lam
def PSensitivityGeo(c,lam,theta,phis1,phi2):
    return 2*pi*c*np.absolute(np.cos(phis1-theta)*np.cos(phi2))/lam
def SVSensitivityGeo(c,lam,theta,phis1,phi2):
    return 2*pi*c*np.absolute(np.cos(phis1-theta)*np.sin(phi2))/lam
def SHSensitivityGeo(c,lam,theta,phis1):
    return 2*pi*c*np.absolute(np.sin(phis1-theta))/lam

def LoveSensitivityStrain(c,lam,theta,phis):
    ccss = np.cos(theta)*np.cos(phis)+np.sin(theta)*np.sin(phis)
    sccs = np.sin(theta)*np.cos(phis)-np.cos(theta)*np.sin(phis)
    return np.absolute(c*ccss*sccs*(2*pi/lam)**2)
def RayleighSensitivityStrain(c,lam,theta,phis):
    ccss = (np.cos(theta)*np.cos(phis))+(np.sin(theta)*np.sin(phis))
    return np.absolute(c*(2*pi*ccss/lam)**2)  
def PSensitivityStrain(c,lam,theta,phis1,phi2):
    return np.absolute(c*(2*pi * np.cos(phis1-theta) * np.cos(phi2)/lam)**2)
def SVSensitivityStrain(c,lam,theta,phis1,phi2):
    return np.absolute(0.5*c*np.sin(2*phi2)*(2*pi*np.cos(phis1-theta)/lam)**2)
def SHSensitivityStrain(c,lam,theta,phis1):
    return np.absolute(c*0.5*(2*pi/lam)**2 * np.sin(2*(phis1-theta)) * np.cos(phi2))

def LoveSensitivityDAS(c,lam,theta,phis,g):
    sccs = np.sin(theta)*np.cos(phis)-np.cos(theta)*np.sin(phis)
    ccss = np.cos(theta)*np.cos(phis)+np.sin(theta)*np.sin(phis)
    return 4*pi*c*np.absolute( sccs * np.sin(pi*g*ccss/lam) )/(lam*g)
def RayleighSensitivityDAS(c,lam,theta,phis,g):
    ccss = np.cos(theta)*np.cos(phis)+np.sin(theta)*np.sin(phis)
    return 4*pi*c*np.absolute(ccss * np.sin(pi*g*ccss/lam) )/ (lam*g)
def PSensitivityDAS(c,lam,theta,phis1,phi2,g):
    return np.absolute(4*pi*c*np.cos(phis1-theta)*np.cos(phi2)*np.sin(0.5*g*2*pi*np.cos(phis1-theta)*np.cos(phi2)/lam)/(lam*g))
def SVSensitivityDAS(c,lam,theta,phis1,phi2,g):
    return np.absolute(c*(2*pi/lam)**2 * np.cos(phis1-theta)*np.sin(2*phi2)*np.sin(pi*g*np.cos(phis1-theta)*np.cos(phis)/lam)/(g*np.cos(phi2)))
def SHSensitivityDAS(c,lam,theta,phis1,g):
    return np.absolute(2*pi*c*np.sin(phis1-theta)*np.sin(pi*g*np.cos(phis1-theta)*np.cos(phi2)/lam)/(lam*g))






c = 400 # velocity
frequencies = np.arange(9,40,10.0)
theta = 0




# #################### Plots of particle velocity, point-wise strain rate and DAS 10 m gauge together #############
ff, axarr = plt.subplots(4, 3, subplot_kw=dict(projection='polar'))
for a in range(4):
        anglelabels = ['0$^\circ$','45$^\circ$','','135$^\circ$','180$^\circ$','225$^\circ$','','315$^\circ$']
        axarr[a,0].set_ylim(0,250)
        axarr[a,0].set_yticks([100,200])
        labels = ['100','200'] 
        axarr[a,0].set_yticklabels(labels,fontsize='x-small')
        axarr[a,0].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,0].annotate("$\lambda$="+str(int(c/frequencies[a])),xy=(-5*pi/6, 250),xytext=((-15*pi/16, 600)))
        axarr[a,1].set_ylim(0,160)
        axarr[a,1].set_yticks([60,120])
        labels = ['60','120'] 
        axarr[a,1].set_yticklabels(labels,fontsize='x-small')
        axarr[a,1].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,2].set_ylim(0,50)
        axarr[a,2].set_yticks([20,40])
        labels = ['20','40'] 
        axarr[a,2].set_yticklabels(labels,fontsize='x-small')
        axarr[a,2].set_xticklabels(anglelabels,fontsize='x-small')


#plt.figure()
for idx,f in enumerate(frequencies):
    nphis = 100+idx*30
    phis = np.arange(0,2*pi,2*pi/nphis)
    lam = c/f # wavelengths
    sensitivities = RayleighSensitivityGeo(c,lam,theta,phis)
    colorFrac = float(idx)/float(frequencies.size)
    axarr[idx, 0].plot(phis, sensitivities,c='g')
    #plt.polar(phis,sensitivities,c=[1-colorFrac*0.8,1,0],label = str(int(lam))+" m, Rayleigh")
    sensitivities = LoveSensitivityGeo(c,lam,theta,phis)
    axarr[idx, 0].plot(phis, sensitivities,c='r')
    #plt.polar(phis,sensitivities,c=[colorFrac,0,1-colorFrac],label = str(int(lam))+" m, Love")
axarr[0, 0].set_title('Particle Velocity')

#plt.clf()
#plt.figure()
for idx,f in enumerate(frequencies):
    nphis = 100+idx*30
    phis = np.arange(0,2*pi,2*pi/nphis)
    lam = c/f # wavelengths
    sensitivities = RayleighSensitivityStrain(c,lam,theta,phis)
    colorFrac = float(idx)/float(frequencies.size)
    axarr[idx, 1].plot(phis, sensitivities,c='g')
    sensitivities = LoveSensitivityStrain(c,lam,theta,phis)
    axarr[idx, 1].plot(phis, sensitivities,c='r')
axarr[0, 1].set_title('Point-wise Strain Rate')

g = 10.0 # 15 meter gauge length
for idx,f in enumerate(frequencies):
    nphis = 100+idx*30
    phis = np.arange(0,2*pi,2*pi/nphis)
    lam = c/f # wavelengths
    sensitivities = RayleighSensitivityDAS(c,lam,theta,phis,g)
    colorFrac = float(idx)/float(frequencies.size)
    axarr[idx, 2].plot(phis, sensitivities,c='g')
    sensitivities = LoveSensitivityDAS(c,lam,theta,phis,g)
    axarr[idx, 2].plot(phis, sensitivities,c='r')
axarr[0, 2].set_title('DAS, '+str(int(g))+' m gauge')
plt.savefig(fig+'geophone_ptStrRate_DAS_'+str(int(g))+'gauge_basic_sensitivity')










# #################### Plots of DAS 2 m, DAS 5 m, DAS 10 m and DAS 20 m gauge together #############
plt.clf()
gauges = [2.0, 5.0, 10.0, 20.0]
ff, axarr = plt.subplots(4, len(gauges), subplot_kw=dict(projection='polar'))

for a in range(4):
        anglelabels = ['0$^\circ$','45$^\circ$','','135$^\circ$','180$^\circ$','225$^\circ$','','315$^\circ$']
        axarr[a,0].set_ylim(0,160)
        ticks = [60,120]
        axarr[a,0].set_yticks(ticks)
        labels = [str(t) for t in ticks] 
        axarr[a,0].set_yticklabels(labels,fontsize='x-small')
        axarr[a,0].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,0].annotate("$\lambda$="+str(int(c/frequencies[a])),xy=(-5*pi/6, 160),xytext=((-15*pi/16, 360)))
        axarr[a,1].set_ylim(0,120)
        ticks = [50,100]
        axarr[a,1].set_yticks(ticks)
        labels = [str(t) for t in ticks] 
        axarr[a,1].set_yticklabels(labels,fontsize='x-small')
        axarr[a,1].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,2].set_ylim(0,50)
        ticks = [20,40]
        axarr[a,2].set_yticks(ticks)
        labels = [str(t) for t in ticks] 
        axarr[a,2].set_yticklabels(labels,fontsize='x-small')
        axarr[a,2].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,3].set_ylim(0,25)
        ticks = [10,20]
        axarr[a,3].set_yticks(ticks)
        labels = [str(t) for t in ticks] 
        axarr[a,3].set_yticklabels(labels,fontsize='x-small')
        axarr[a,3].set_xticklabels(anglelabels,fontsize='x-small')


for ig,g in enumerate(gauges):
    for idx,f in enumerate(frequencies):
        nphis = 100+idx*30
        phis = np.arange(0,2*pi,2*pi/nphis)
        lam = c/f # wavelengths
        sensitivities = RayleighSensitivityDAS(c,lam,theta,phis,g)
        colorFrac = float(idx)/float(frequencies.size)
        axarr[idx, ig].plot(phis, sensitivities,c='g')
        sensitivities = LoveSensitivityDAS(c,lam,theta,phis,g)
        axarr[idx, ig].plot(phis, sensitivities,c='r')
    axarr[0, ig].set_title(str(int(g))+' m gauge')
plt.savefig(fig+'multi_gauge_DAS_basic_sensitivity')






######### sensitivity plots for P, SV, SH waves ####################################



phi2List = [np.pi/8, np.pi/4, 3*np.pi/8]
# #################### Plots of particle velocity, point-wise strain rate and DAS 10 m gauge together #############
for phi2idx,phi2 in enumerate(phi2List):
    plt.clf()
    ff, axarr = plt.subplots(4, 3, subplot_kw=dict(projection='polar'))
    for a in range(4):
        anglelabels = ['0$^\circ$','45$^\circ$','','135$^\circ$','180$^\circ$','225$^\circ$','','315$^\circ$']
        axarr[a,0].set_ylim(0,250)
        axarr[a,0].set_yticks([100,200])
        labels = ['100','200'] 
        axarr[a,0].set_yticklabels(labels,fontsize='x-small')
        axarr[a,0].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,0].annotate("$\lambda$="+str(int(c/frequencies[a])),xy=(-5*pi/6, 250),xytext=((-15*pi/16, 600)))
        axarr[a,1].set_ylim(0,80)
        axarr[a,1].set_yticks([30,60])
        labels = ['30','60'] 
        axarr[a,1].set_yticklabels(labels,fontsize='x-small')
        axarr[a,1].set_xticklabels(anglelabels,fontsize='x-small')
        axarr[a,2].set_ylim(0,40)
        axarr[a,2].set_yticks([15,30])
        labels = ['15','30'] 
        axarr[a,2].set_yticklabels(labels,fontsize='x-small')
        axarr[a,2].set_xticklabels(anglelabels,fontsize='x-small')


    for idx,f in enumerate(frequencies):
        nphis = 100+idx*30
        phis = np.arange(0,2*pi,2*pi/nphis)
        lam = c/f # wavelengths
        sensitivities = PSensitivityGeo(c,lam,theta,phis,phi2)
        axarr[idx, 0].plot(phis, sensitivities,c='darkorange',linewidth=5)
        sensitivities = SVSensitivityGeo(c,lam,theta,phis,phi2)
        axarr[idx, 0].plot(phis, sensitivities,c='black',linewidth=3)
        sensitivities = SHSensitivityGeo(c,lam,theta,phis)
        axarr[idx, 0].plot(phis, sensitivities,c='blue',linewidth=2)
    axarr[0, 0].set_title('Particle Velocity')


    for idx,f in enumerate(frequencies):
        nphis = 100+idx*30
        phis = np.arange(0,2*pi,2*pi/nphis)
        lam = c/f # wavelengths
        sensitivities = PSensitivityStrain(c,lam,theta,phis,phi2)
        axarr[idx, 1].plot(phis, sensitivities,c='darkorange',linewidth=5)
        sensitivities = SVSensitivityStrain(c,lam,theta,phis,phi2)
        axarr[idx, 1].plot(phis, sensitivities,c='black',linewidth=3)
        sensitivities = SHSensitivityStrain(c,lam,theta,phis)
        axarr[idx, 1].plot(phis, sensitivities,c='blue',linewidth=2)
    axarr[0, 1].set_title('Point-wise Strain Rate')



    g = 10.0 # 15 meter gauge length
    for idx,f in enumerate(frequencies):
        nphis = 100+idx*30
        phis = np.arange(0,2*pi,2*pi/nphis)
        lam = c/f # wavelengths
        sensitivities = PSensitivityDAS(c,lam,theta,phis,phi2,g)
        axarr[idx, 2].plot(phis, sensitivities,c='darkorange',linewidth=5)
        sensitivities = SVSensitivityDAS(c,lam,theta,phis,phi2,g)
        axarr[idx, 2].plot(phis, sensitivities,c='black',linewidth=3)
        sensitivities = SHSensitivityDAS(c,lam,theta,phis,g)
        axarr[idx, 2].plot(phis, sensitivities,c='blue',linewidth=2)
    axarr[0, 2].set_title('DAS, '+str(int(g))+' m gauge')
    plt.savefig(fig+'geophone_ptStrRate_DAS_'+str(int(g))+'gauge_body_wave_sensitivity_phi2_'+str(phi2idx))
