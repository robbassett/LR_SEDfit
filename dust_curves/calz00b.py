import numpy as np

def calz00b(ll):
    si = len(ll)
    dc = np.zeros(si)

    t = np.where((ll >= 0.63)&(ll <= 2.2))[0]
    dc[t] = (1.04/ll[t])-1.857
    t1 = np.where((ll >= 0.15)&(ll < .63))[0]
    dc[t1] = (1.509/ll[t1])-(0.198/(ll[t1]*ll[t1]))+(0.011/(ll[t1]*ll[t1]*ll[t1]))-2.156
    t2 = np.where(ll < 0.15)[0]
    dc[t2] = 4.126+(0.931/ll[t2])

    fc = dc[t1[0]]/dc[t2[-1]]

    dc[t2]*=fc
    
    
    return dc
