from matplotlib import pyplot as plt
from card_dc import *
from calz94  import *
from calz00  import *
from wild_av import *
from SMCdc   import *
from calz00b import *
from reddy16 import *

def dcnorm(tv,tb,dc):

  nv = dc[tb]-dc[tv]
  si=len(dc)
  for i in range(si): dc[i]/=nv

  return dc

ll = np.arange(0.05,2.,.0001)
tv = np.where(abs(ll-.551) == min(abs(ll-.551)))
tb = np.where(abs(ll-.445) == min(abs(ll-.445)))
tuv = np.where(abs(ll-.15) == min(abs(ll-.15)))

rv=3.1
ebv=.5

dc1 = calz00(ll)
dc2=reddy16(ll)
dc3=calz00b(ll)

dc1 = dcnorm(tv[0][0],tb[0][0],dc1)
dc2= dcnorm(tv[0][0],tb[0][0],dc2)
dc3= dcnorm(tv[0][0],tb[0][0],dc3)

#print dc1[tb[0][0]]-dc1[tv[0][0]]
#print dc2[tb[0][0]]-dc2[tv[0][0]]
#print dc3[tb[0][0]]-dc3[tv[0][0]]

ll*=10000.

F = plt.figure()
ax = F.add_subplot(111)
ax.plot(ll,dc1,'k-')
ax.plot(ll,dc2,'r-')
ax.plot(ll,dc3,'b-')
ax.axvline(ll[tv[0][0]])
ax.set_xlim(500,5000)
plt.show()
