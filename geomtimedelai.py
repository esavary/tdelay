from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle as pl


axis_font = {'fontname': 'Arial', 'size': '18'}

c = 1.  # unit light day/day
ldaystopixels=37.002954816# change this parameter for each map

def rt(t):  # radius of the supernova as a function of t
    rt = 10 ** (-5) + (t+30) / 30.

    return rt #returns lightdays in pixels


def flux(t):
    '''
    magnitude = 2 * 10 ** (-6) * (t - 20) ** 4 - 0.00022 * (t - 20) ** 3 + 0.0077 * (t - 20) ** 2 + 0.0025 * (
    t - 20) - 19.5  # approximate fit of a supernova template
    f = 10 ** (magnitude / (-2.5)) * 4000
    '''
    f=3.94651870932e-07*t**5-0.000113731703456*t**4+0.0111317577292*t**3-0.381176185494*t**2-0.112832155346*t+115.424606424
    # equation to transform magnitude into flux, 4000~= zero point for blue wavelength

    return f
def returncroppedmap(xcut,ycut,xmin, xmax, ymin, ymax,gridstep,nomagnification= False):
    repetitionfactor=1./gridstep
    if nomagnification==False:
        hdul = fits.open('map.fits')
        map=hdul[0].data
        croppedmap=map[xcut:xcut+(xmax-xmin),ycut:ycut+(ymax-ymin)]

    else:
        print 'nomagnification'
        croppedmap=np.ones((xmax-xmin,ymax-ymin))

    multmapx = np.repeat(np.array(croppedmap),repetitionfactor,axis=1)
    multmaptot = np.repeat(multmapx, repetitionfactor, axis=0)
    return multmaptot



# normalisation factor, the integral is computed by taking the sum on all pixels in the chosen grid
def normint(gridstep, xmin, xmax, ymin, ymax,xcut,ycut, t, z):
    magnificationmap = returncroppedmap(xcut, ycut, xmin, xmax, ymin, ymax, gridstep, nomagnification=False)
    inttot = []

    for x in np.arange(xmin, xmax, gridstep):
        for y in np.arange(ymin, ymax, gridstep):
            if (x < rt(t) * ldaystopixels) and (np.power(y, 2) < np.power(rt(t), 2) - np.power(x, 2)):
                delai = np.sqrt(np.power(rt(t), 2) - np.power(x / ldaystopixels, 2) - np.power(y / ldaystopixels, 2))#positive delai
                inttot.append(flux(t + delai)*magnificationmap[x,y] / (pi * np.power(rt(t+delai),2)))
    tot = np.array(inttot)
    return np.sum(tot)


# non normalized mean delta t, the integral is computed by taking the sum on all pixels in the chosen grid
def nonormmeandelaiint(gridstep, xmin, xmax, ymin, ymax,xcut,ycut, t, z):
    magnificationmap=returncroppedmap(xcut,ycut,xmin, xmax, ymin, ymax,gridstep,nomagnification= False)
    inttot = []
    for x in np.arange(xmin, xmax, gridstep):
        for y in np.arange(ymin, ymax, gridstep):
            if (x < rt(t)*ldaystopixels) and (np.power(y,2)< np.power(rt(t),2) - np.power(x,2)):
                delai = np.sqrt(np.power(rt(t),2) - np.power(x/ldaystopixels,2) - np.power(y/ldaystopixels,2))
                inttot.append(flux(t + delai) * delai *magnificationmap[x,y]/ (c * pi * 2 * np.power(rt(t+delai),2)))
    tot = np.array(inttot)
    return (1 + z) * np.sum(tot)


def meandelai(gridstep, xmin, xmax, ymin, ymax,xcut,ycut,t, z):
    # mean delai normalized
    meandelai = nonormmeandelaiint(gridstep, xmin, xmax, ymin, ymax,xcut,ycut, t, z) * 1. / normint(gridstep, xmin, xmax, ymin,
                                                                                          ymax,xcut,ycut, t, z)
    return meandelai
def returnmodifiedflux(gridstep, xmin, xmax, ymin, ymax, xcut,ycut,t, z,magnification=True):
    inttot = []
    if magnification==True:
        magnificationmap = returncroppedmap(xcut, ycut, xmin, xmax, ymin, ymax, gridstep, nomagnification=False)
        for x in np.arange(xmin, xmax, gridstep):
            for y in np.arange(ymin, ymax, gridstep):
                if (x < rt(t) * ldaystopixels) and (np.power(y, 2) < np.power(rt(t)* ldaystopixels, 2) - np.power(x, 2)):
                    delai = np.sqrt(np.power(rt(t), 2) - np.power(x / ldaystopixels, 2) - np.power(y / ldaystopixels, 2))
                    inttot.append(flux(t + delai)*magnificationmap[x,y] / (pi * np.power(rt(t+delai),2)))
        tot = np.array(inttot)*gridstep**2
    else:
        for x in np.arange(xmin, xmax, gridstep):
            for y in np.arange(ymin, ymax, gridstep):
                if (x < rt(t) * ldaystopixels) and (np.power(y, 2) < np.power(rt(t)* ldaystopixels, 2) - np.power(x, 2)):
                    delai = np.sqrt(np.power(rt(t), 2) - np.power(x / ldaystopixels, 2) - np.power(y / ldaystopixels, 2))
                    inttot.append(flux(t + delai) / (pi * np.power(rt(t+delai),2)))
        tot = np.array(inttot)*gridstep**2

    return np.sum(tot)


####display function
def plotmeandt(gridstep, xmin, xmax, ymin, ymax, xcut,ycut,timetab, z):  # plot the mean time delai for an array of time values
    mdt = []
    for t in timetab:
        mdt.append(meandelai(gridstep, xmin, xmax, ymin, ymax, xcut,ycut, t, z))
    fig=plt.figure()
    plt.plot(timetab, mdt)
    plt.tick_params(labelsize=16)
    plt.xlabel(' t [ld]', **axis_font)
    plt.ylabel('<$\delta t)$> [ld]', **axis_font)
    pl.dump(fig, file('meandtbigmagn.pickle', 'w'))
    plt.show()


def drawcolormap(gridstep, xmin, xmax, ymin, ymax, t):
    # draw the delay map at a given t +return the array containing the values.
    dimmapx = int((xmax - xmin) * 1. / gridstep)
    dimmapy = int((ymax - ymin) * 1. / gridstep)
    delaimap = np.zeros((dimmapx, dimmapy))
    for x in np.arange(xmin, xmax, gridstep):
        for y in np.arange(ymin, ymax, gridstep):
            if ((np.power(x,2) + np.power(y,2)) < np.power(rt(t)*ldaystopixels,2)) and (x < rt(t)*ldaystopixels):
                deltat = -np.sqrt(np.power(rt(t),2) - np.power(x/ldaystopixels,2) - np.power(y/ldaystopixels,2))
                delaimap[int((x + (xmax - xmin) / 2) * 1. / gridstep),int( (y + (ymax - ymin) / 2) * 1. / gridstep)] = deltat

    plt.figure()
    p = plt.pcolormesh(np.arange(xmin, xmax, gridstep), np.arange(ymin, ymax, gridstep), delaimap ,vmin=-1.1*rt(40),vmax=0)
    plt.title('time delay at t=' + str(t) + ' ld [ld]')
    plt.colorbar(p)
    plt.show()
    return delaimap


if __name__ == '__main__':
    z = 1.
    # min and max coordinates of the grid:
    minmax=200
    xmax = minmax
    ymax = minmax
    xmin = -minmax
    ymin = -minmax
    xcut=4100
    ycut=5700
    timetab = np.arange(-10, 20, 0.4)
    gridstep = 1
    #plotmeandt(gridstep, xmin, xmax, ymin, ymax,xcut,ycut, timetab,  z)
    #delaimap = drawcolormap(gridstep, xmin, xmax, ymin, ymax, timetab[10])
    essai=returncroppedmap(100,100, xmin, xmax, ymin, ymax,gridstep)

    fluxtrue=[]
    modified=[]

    for t in timetab:
        fluxtrue.append(flux(t))
        modified.append(returnmodifiedflux(gridstep, xmin, xmax, ymin, ymax, xcut,ycut,t, z,magnification=True))
    fig=plt.figure()
    plt.plot(timetab,fluxtrue,label='analytique')
    plt.plot(timetab, modified, label='modified with no magnification')
    plt.legend()

    pl.dump(fig, file('truefluxcompbigmagn.pickle', 'w'))
    plt.show()

    fluxmaxindice=np.argmax(fluxtrue)
    fluxmaxmodifiedindice=np.argmax(modified)
    print 'original max:',timetab[fluxmaxindice],' modified max:',timetab[fluxmaxmodifiedindice]
    np.save('truearray01magn.npy', fluxtrue, allow_pickle=True, fix_imports=True)
    np.save('modifiedarray01magn.npy', modified, allow_pickle=True, fix_imports=True)
    np.save('timerray015magn.npy',timetab, allow_pickle=True, fix_imports=True)
