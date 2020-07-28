import numpy as np

def reddy16(ll):

    dc = np.zeros(len(ll))

    x1 = 1./ll
    x2 = 1./(ll*ll)
    x3 = 1./(ll*ll*ll)

    t1 = np.where((ll >= 0.15)&(ll < 0.6))[0]
    dc[t1] = (-5.726)+(4.004*x1[t1])-(0.525*x2[t1])+(0.029*x3[t1])+2.505
    t2 = np.where((ll >= 0.6)&(ll < 2.85))[0]
    dc[t2] = (-2.672)-(0.010*x1[t2])+(1.532*x2[t2])-(0.412*x3[t2])+2.505
    t3 = np.where(ll < 0.15)[0]
    dc[t3] = (2.191)+(0.974*x1[t3])

    return dc
