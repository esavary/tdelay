from sympy import *
import numpy as np
import matplotlib.pyplot as plt
axis_font = {'fontname':'Arial', 'size':'18'}

c=1.#unit light day/day
def rt(t): #radius of the supernova as a function of t
    rt=10**(-5)+t/30.
    return rt
def flux(t):
    magnitude = 2 * 10 ** (-6) * (t - 20) ** 4 - 0.00022 * (t - 20) ** 3 + 0.0077 * (t - 20) ** 2 + 0.0025 * (t - 20) - 19.5 #approximate fit of a supernova template
    f = 10 **(magnitude / (-2.5)) * 4000 # equation to transform magnitude into flux, 4000~= zero point for blue wavelength
    return f

#normalisation factor, the integral is computed by taking the sum on all pixels in the chosen grid
def normint(gridstep,xmin,xmax,ymin,ymax,t,z):
    inttot=[]
    for x in np.arange(xmin,xmax,gridstep):
        for y in np.arange(ymin,ymax,gridstep):
            if (x < rt) and (y ** 2 < rt(t) ** 2 - x ** 2):
                delai=-np.sqrt(rt(t)**2-x**2-y**2)
                inttot.append(flux(t+delai)/(pi*2*rt(t+delai)**2))
    tot=np.array(inttot)
    return np.sum(tot)

#non normalized mean delta t, the integral is computed by taking the sum on all pixels in the chosen grid
def nonormmeandelaiint(gridstep,xmin,xmax,ymin,ymax,t,z):
    inttot=[]
    for x in np.arange(xmin,xmax,gridstep):
        for y in np.arange(ymin,ymax,gridstep):
            if (x<rt) and (y**2<rt(t)**2-x**2):
                delai=-np.sqrt(rt(t)**2-x**2-y**2)
                inttot.append(flux(t+delai)*delai/(c*pi*2*rt(t+delai)**2))
    tot=np.array(inttot)
    return (1+z)*np.sum(tot)

def meandelai(gridstep,xmin,xmax,ymin,ymax,t,z):
    #mean delai normalized 
    meandelai=nonormmeandelaiint(gridstep,xmin,xmax,ymin,ymax,t,z)*1./normint(gridstep,xmin,xmax,ymin,ymax,t,z)
    return meandelai

####display function
def plotmeandt(gridstep,xmin,xmax,ymin,ymax,timetab,z):#plot the mean time delai for an array of time values
    mdt=[]
    for t in timetab:
        mdt.append(meandelai(gridstep,xmin,xmax,ymin,ymax,t,z))

    plt.plot(timetab,mdt)
    plt.tick_params(labelsize=16)
    plt.xlabel(' t [ld]',** axis_font)
    plt.ylabel('<$\delta t)$> [ld]',** axis_font)
    plt.show()
def drawcolormap(gridstep,xmin,xmax,ymin,ymax,t): 
    #draw the delay map at a given t +return the array containing the values.
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
    
    
    z=1.
    # min and max coordinates of the grid:
    xmax=5
    ymax=5
    xmin=-5
    ymin=-5
    timetab=np.arange(0,40,1)
    gridstep=0.01
    plotmeandt(gridstep,xmin,xmax,ymin,ymax,timetab,z)
    delaimap=drawcolormap(gridstep,xmin,xmax,ymin,ymax,timetab[30])

