import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np
import glob
import sys

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

import BPASS_tools as bpt
import LR_mod as lrm
import dust_curves as dm
from load_photometry import Photometry as pht
from LRsed import *


dust_dic = {'SMC':'SMC',
            'Reddy (2016)':'reddy16',
            'Wild (2011)':'wild',
            'Calzetti (2000)':'calz00b'}

class startup_gui(QMainWindow):

    def __init__(self):
        super().__init__()

        self.filters = None
        self.photometry = None
        self.BPASS_folder = None
        self.dust_mod = None
        self.initUI()

    def initUI(self):

        exitAct = QAction(QIcon('assets/exit.png'), 'Exit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.triggered.connect(qApp.quit)

        loadfilt = QAction(QIcon('assets/filt.png'), 'Load Filters', self)
        loadfilt.triggered.connect(self.load_filters)

        loadphot = QAction(QIcon('assets/tele.png'), 'Load Photometry', self)
        loadphot.triggered.connect(self.load_photometry)

        bpfold = QAction(QIcon('assets/bp.png'),'Select BPASS Folder',self)
        bpfold.triggered.connect(self.get_BPfold)

        gobutt = QAction(QIcon('assets/go.png'),'Start SED Fitting!',self)
        gobutt.triggered.connect(self.go)

        self.toolbar = QToolBar(self)
        self.addToolBar(self.toolbar)
        self.toolbar.addAction(loadfilt)
        self.toolbar.addAction(loadphot)
        self.toolbar.addAction(bpfold)
        self.toolbar.addAction(gobutt)
        self.toolbar.addAction(exitAct)
        
        self.setGeometry(100, 100, 320, 185)
        self.setWindowTitle('Startup Widget')

        dlab = QLabel(self)
        dlab.setText('Dust Curve:')
        dlab.setGeometry(QRect(15,40,130,50))
        self.combo_dust = QComboBox(self)
        self.combo_dust.setGeometry(QRect(85,40,130,50))
        self.combo_dust.setObjectName('dust_mod')
        for k in dust_dic.keys(): self.combo_dust.addItem(k)

        rvlab = QLabel(self)
        rvlab.setText('Rv')
        rvlab.setGeometry(QRect(220,40,30,50))
        self.Rvbox = QLineEdit(self)
        self.Rvbox.setGeometry(QRect(245,55,60,20))
        self.Rvbox.setText('2.74')

        filt_lab = QLabel(self)
        filt_lab.setText('             Filters?')
        filt_lab.setGeometry(QRect(15,75,120,50))
        filt_lab.setFont(QFont("Times", 18,QFont.Bold))
        self.ch1 = QLabel(self)
        self.ch1.setText('         ')
        self.ch1.setPixmap(QPixmap('assets/x.png'))
        self.ch1.setGeometry(QRect(150,88,20,20))

        phot_lab = QLabel(self)
        phot_lab.setText('    Photometry?')
        phot_lab.setGeometry(QRect(15,105,150,50))
        phot_lab.setFont(QFont("Times", 18,QFont.Bold))
        self.ch2 = QLabel(self)
        self.ch2.setText('         ')
        self.ch2.setPixmap(QPixmap('assets/x.png'))
        self.ch2.setGeometry(QRect(150,118,20,20))

        fold_lab = QLabel(self)
        fold_lab.setText('BPASS Folder?')
        fold_lab.setGeometry(QRect(15,135,150,50))
        fold_lab.setFont(QFont("Times", 18,QFont.Bold))
        self.ch3 = QLabel(self)
        self.ch3.setText('         ')
        self.ch3.setPixmap(QPixmap('assets/x.png'))
        self.ch3.setGeometry(QRect(150,148,20,20))
        
        self.show()

    def load_filters(self):
        msg = QMessageBox()
        msg.setWindowTitle('Load Filters')
        msg.setText('Select the filter definition file (.npy)')
        msg.exec_()
        
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Select Filter File", "","Filter Files (*.npy)", options=options)
        if fileName:
            self.filters = np.load(fileName,allow_pickle=True).item()
            self.ch1.setPixmap(QPixmap('assets/check.png'))

    def load_photometry(self):
        msg = QMessageBox()
        msg.setWindowTitle('Load Photometry')
        msg.setText('Select the photometry text file')
        msg.exec_()
        
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        photfile, _ = QFileDialog.getOpenFileName(self,"Select Photometric Data File", "","All Files (*.*)", options=options)
        if photfile:
            msg = QMessageBox()
            msg.setWindowTitle('Load Photometry')
            msg.setText('Select the photometry definition file')
            msg.exec_()
            
            options = QFileDialog.Options()
            options |= QFileDialog.DontUseNativeDialog
            photdef, _ = QFileDialog.getOpenFileName(self,"Select Photometry Definition File", "","All Files (*.*)", options=options)

            if photdef:
                self.photometry = pht(photfile,photdef)
                self.ch2.setPixmap(QPixmap('assets/check.png'))

    def get_BPfold(self):
        msg = QMessageBox()
        msg.setWindowTitle('Select BPASS Folder')
        msg.setText('Select folder containing your BPASS spectra')
        msg.exec_()

        options = QFileDialog.getExistingDirectory(self)
        self.BPASS_folder = options
        self.ch3.setPixmap(QPixmap('assets/check.png'))

    def go(self):
        self.Rv = float(self.Rvbox.text())
        self.bdc = dm.DustCurve(dust_dic[self.combo_dust.currentText()])
        fit_one_SED(14759,
                    self.bdc,
                    self.filters,
                    self.BPASS_folder,
                    self.photometry,
                    Rv=self.Rv
        )

                
def main():
    app = QApplication(sys.argv)
    ex = startup_gui()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
