def SED_loss(modf,obsf,obse):
    tobs = np.where(obsf > 0)[0]

    diffs = (((modf[tobs]-obsf[tobs])/obse[tobs])**2.)
    return np.sum(diffs)

def linear_ffact(x,a):
    return x*a

def h(theta,x):
    return np.sum((theta*x.transpose()).transpose(),axis=0)

def cost(theta,x,y,yerr):
    m = float(len(y))
    return np.sum(((h(theta,x)-y)/yerr)**2.)/(2.*m)

def grad(theta,x,y,yerr):
    grd = np.zeros(len(theta))
    for ii in range(len(theta)):
        tm      = (h(theta,x)-y)/yerr
        grd[ii] = np.sum(tm*x[ii])

    return grd
