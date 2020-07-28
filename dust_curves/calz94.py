import numpy as np

def calz94(ll):
    Qe = -2.156 + (1.509/ll) - (0.198/(ll**2.)) + (0.011/(ll**3.))
    return Qe
