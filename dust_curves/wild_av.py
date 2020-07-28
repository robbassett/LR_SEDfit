import numpy as np

def Qlam(lam,p,ba,sSFR,n):

    l1 = .2175
    l2 = 0.30
    l3 = 0.8

    ba-=.6
    sSFR+=9.5

    s1 = p[0]+p[1]*ba+p[2]*sSFR
    s2 = p[3]+p[4]*ba+p[5]*sSFR
    s3 = p[6]+p[7]*ba+p[8]*sSFR
    s4 = 1.6

    lb1 = l1
    lb2 = ((lb1**s2)/(l2**(s2-s3)))**(1./s3)
    lb3 = ((lb2**s3)/(l3**(s3-s4)))**(1./s4)

    Ql  = ( ((lam/lb1)**(n*s1)) + ((lam/lb1)**(n*s2)) + ((lam/lb2)**(n*s3)) + ((lam/lb3)**(n*s4)) )**(-1./n)

    return Ql

def wild_av(lam):

    avh=1.
    P_high = [.2,.9,.1,1.1,.3,-.2,1.3,.6,.1]

    Qt = Qlam(lam,P_high,.6,-9.5,20.)

    return Qt
    
