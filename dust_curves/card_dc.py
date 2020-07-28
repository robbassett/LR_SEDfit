import numpy as np

def IRc(x):
    y=(x**(1.61))
    ax=(0.574)*y
    bx=(-0.527)*y

    return [ax,bx]

def OPc(x):

    y = (x-1.82)
    ax = (1.)+(.17699*y)-(0.50447*(y**2.))-(0.02427*(y**3.))+(0.72085*(y**4))+(0.01979*(y**5))-(0.7753*(y**6.))+(0.32999*(y**7))
    bx = (1.41338*y)+(2.28305*(y**2.))+(1.07233*(y**3.))-(5.38434*(y**4.))-(0.62251*(y**5.))+(5.30260*(y**6.))-(2.09002*(y**7.))

    return [ax,bx]

def Fa(x):

    ya = (x-5.9)
    ay = (-.04473)*(ya**2.)
    by = (-0.009779)*(ya**3.)

    return ay+by

def Fb(x):

    yb = (x-5.9)
    at = 0.213*(yb**2.)
    bt = 0.1207*(yb**3.)
    return at+bt

def UVc(x):

    at = (x-4.67)**2.
    bt = (x-4.62)**2.
    
    at+=0.341
    bt+=0.263

    F1=0.
    F2=0.
    if x > 5.9:
        F1=Fa(x)
        F2=Fb(x)

    ax = 1.752-(0.316*x)-(1.04/at)+F1
    bx = (-3.09)+(1.825*x)+(1.206/bt)+F2

    return [ax,bx]

def av(tmp,Rv):

    return tmp[0]+(tmp[1]/Rv)

def card_dc(ll,Rv):

    si = len(ll)
    dc = np.zeros(si)

    breaks = [3.3,100.,1.1,0.3]
    for i in range(si):

        x=1./(ll[i])

        if x > breaks[3] and x < breaks[2]: tmp=IRc(x)

        if x > breaks[2] and x < breaks[0]: tmp=OPc(x)

        if x > breaks[0]: tmp=UVc(x)

        dc[i]=av(tmp,Rv)

    #dc*=(Rv*ebv)

    return dc

        
