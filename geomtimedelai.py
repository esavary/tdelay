from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *
axis_font = {'fontname':'Arial', 'size':'18'}


def rt(t):
    rt=10**(-5)+t/30.
    return rt
def flux(t):
    magnitude = 2 * 10 ** (-6) * (t - 20) ** 4 - 0.00022 * (t - 20) ** 3 + 0.0077 * (t - 20) ** 2 + 0.0025 * (t - 20) - 19.5
    f = 10 **(magnitude / (-2.5)) * 4000
    return f

def normint(gridstep,xmin,xmax,ymin,ymax,t):
    inttot=[]
    for x in np.arange(xmin,xmax,gridstep):
        for y in np.arange(ymin,ymax,gridstep):
            if (x < rt) and (y ** 2 < rt(t) ** 2 - x ** 2):
                delai=-np.sqrt(rt(t)**2-x**2-y**2)
                inttot.append(flux(t+delai)/(pi*2*rt(t+delai)**2))
    tot=np.array(inttot)
    return np.sum(tot)

def nonormmeandelaiint(gridstep,xmin,xmax,ymin,ymax,t):
    inttot=[]
    for x in np.arange(xmin,xmax,gridstep):
        for y in np.arange(ymin,ymax,gridstep):
            if (x<rt) and (y**2<rt(t)**2-x**2):
                delai=-np.sqrt(rt(t)**2-x**2-y**2)
                inttot.append(flux(t+delai)*delai/(pi*2*rt(t+delai)**2))
    tot=np.array(inttot)
    return np.sum(tot)

def meandelai(gridstep,xmin,xmax,ymin,ymax,t):

    meandelai=nonormmeandelaiint(gridstep,xmin,xmax,ymin,ymax,t)*1./normint(gridstep,xmin,xmax,ymin,ymax,t)
    return meandelai

####display
def plotmeandt(gridstep,xmin,xmax,ymin,ymax,timetab):
    mdt=[]
    for t in timetab:
        mdt.append(meandelai(gridstep,xmin,xmax,ymin,ymax,t))

    plt.plot(timetab,mdt)
    plt.tick_params(labelsize=16)
    plt.xlabel(' t [ld]',** axis_font)
    plt.ylabel('<$\delta t)$> [ld]',** axis_font)
    plt.show()
def drawcolormap(gridstep,xmin,xmax,ymin,ymax,t): #plot the map +return the array
    dimmapx=int((xmax-xmin)*1./gridstep)
    dimmapy = int((ymax - ymin) * 1. / gridstep)
    delaimap = np.zeros((dimmapx, dimmapy))
    for x in np.arange(xmin, xmax, gridstep):
        for y in np.arange(ymin, ymax, gridstep):
            if (x ** 2 + y ** 2 < rt(t) ** 2) and (x<rt(t)):
                deltat = -np.sqrt(rt(t) ** 2 - x ** 2 - y ** 2)
                delaimap[(x + (xmax-xmin)/2) * 1./gridstep, (y + (ymax-ymin)/2) * 1./gridstep] = deltat

    plt.figure()
    p = plt.pcolormesh(np.arange(xmin, xmax, gridstep), np.arange(ymin, ymax, gridstep), delaimap)
    plt.title('time delay at t='+str(t)+' ld [ld]')
    plt.colorbar(p)
    plt.show()
    return delaimap


if __name__ == '__main__':

    xmax=5
    ymax=5
    xmin=-5
    ymin=-5
    timetab=np.arange(0,40,1)
    gridstep=0.01
    plotmeandt(gridstep,xmin,xmax,ymin,ymax,timetab)
    delaimap=drawcolormap(gridstep,xmin,xmax,ymin,ymax,timetab[30])

