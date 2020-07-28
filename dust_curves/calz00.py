import numpy as np

def calz00(ll):
    si = len(ll)
    dc = np.zeros(si)
    for i in range(si):
        xt=1./(ll[i])
        fc=0.
        if ll[i] >= .63 and ll[i] <= 2.2:
            fc=1.04*xt
            fc-=1.857

        #if ll[i] >= .12 and ll[i] <= 0.63:
        if ll[i] <= 0.63:
            fc1=1.509*xt
            fc2=0.198*(xt**2.)
            fc3=0.011*(xt**3.)
            fc = fc1-fc2+fc3-2.156

        fc*=1.17
        if fc != 0: fc+=1.78

        dc[i]=fc

    return dc
