import numpy as np

class Photometry(object):

    def __init__(self,input_file,label_file):
        phot = np.loadtxt(input_file,float)
        labs = np.loadtxt(label_file,str)

        self.phot = {}
        for i,ID in enumerate(phot[:,0].astype(int)):
            self.phot[str(ID)] = {}
            for j,filt in enumerate(labs[:,0]):
                self.phot[str(ID)][filt] = [phot[i,int(labs[j,1])],phot[i,int(labs[j,2])]]

    def ID_2_arrays(self,ID,filters):

        tmf = self.phot[str(ID)]
        wav = np.zeros(len(list(tmf.keys())))
        flx = np.zeros(len(list(tmf.keys())))
        err = np.zeros(len(list(tmf.keys())))
        for i,k in enumerate(tmf.keys()):
            wav[i] = (np.average(filters[k]['wav'],weights=filters[k]['T']))
            flx[i] = tmf[k][0]
            err[i] = tmf[k][1]

        return wav,flx,err
        
