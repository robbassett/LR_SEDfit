import numpy as np
from flux_converters import *

class Photometry(object):

    def __init__(self,input_file,label_file):
        phot = np.loadtxt(input_file,float)
        labs = np.loadtxt(label_file,str)
        flx_ver = open(label_file,'r')
        for i in range(2): line=flx_ver.readline()
        flx_ver = line.split(' ')[-2]

        self.data = {}
        for i,ID in enumerate(phot[:,0].astype(int)):
            for j,filt in enumerate(labs[:,0]):
                if j == 0:
                    self.data[str(ID)] = {'redshift':0.0,'phot':{}}
                elif j == 1:
                    self.data[str(ID)]['redshift'] = phot[i,int(labs[j,1])]
                else:
                    in1,in2 = phot[i,int(labs[j,1])],phot[i,int(labs[j,2])]
                    if flx_ver == 'AB':
                        tflx = AB_mag(in1)
                        terr = tflx-AB_mag(in1+in2)
                    if flx_ver == 'Zfourge':
                        tflx,terr = ZFOURGE_flx(in1),ZFOURGE_flx(in2)
                    if flx_ver == 'Jy':
                        tflx,terr = Jy(in1),Jy(in2)
                    else:
                        tflx,terr = in1,in2
                    self.data[str(ID)]['phot'][filt] = [tflx,terr]

    def get_fit_filters(self,ID,filters):

        tmf = self.data[str(ID)]['phot']
        filters_match = {}
        for i,k in enumerate(tmf.keys()):
            if tmf[k][0] > 0:
                filters_match[k] = {'wav':filters[k]['wav'],'T':filters[k]['T']}

        return list(filters_match.keys())
        
