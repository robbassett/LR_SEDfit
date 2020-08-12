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
from plotting_scripts import *


dust_dic = {'Reddy (2016)':'reddy16',
            'SMC':'SMC',
            'Wild (2011)':'wild',
            'Calzetti (2000)':'calz00b'}

class sed_window(QMainWindow):
    def __init__(self, outputs, parent=None):
        super(sed_window, self).__init__(parent)
        self.outputs = outputs

        self.initUI()

    def initUI(self):
        self.setGeometry(0, 250, 600, 75)
        self.setWindowTitle('SED Fitting Window')

        idlab = QLabel(self)
        idlab.setText('Galaxy ID:')
        idlab.setGeometry(QRect(15,10,130,50))
        self.combo_id = QComboBox(self)
        self.combo_id.setGeometry(QRect(85,10,130,50))
        self.combo_id.setObjectName('tmp_ID')
        for k in self.outputs.keys(): self.combo_id.addItem(k)

        spbutt = QPushButton('Show SED',self)
        spbutt.setToolTip('Display Spectrum')
        spbutt.setGeometry(QRect(230,20,100,25))
        spbutt.clicked.connect(self.disp_spec)

        shbutt = QPushButton('Show SF-hist.',self)
        shbutt.setToolTip('Display Spectrum')
        shbutt.setGeometry(QRect(340,20,120,25))
        shbutt.clicked.connect(self.disp_sfh)

        svbutt = QPushButton('Save Output',self)
        svbutt.setToolTip('Display Spectrum')
        svbutt.setGeometry(QRect(470,20,100,25))
        svbutt.clicked.connect(self.save_out)

    def disp_spec(self):
        SED_plot(self.outputs,self.combo_id.currentText())

    def disp_sfh(self):
        SFhist_plot(self.outputs,self.combo_id.currentText())

    def save_out(self):
        text, okPressed = QInputDialog.getText(self, "Save Outputs","Filename (no ext.):", QLineEdit.Normal, "")
        if okPressed and text != '':
            np.save(f'./{text}.npy',self.outputs)

class phd_window(QMainWindow):
    def __init__(self, filters, parent=None):
        super(phd_window, self).__init__(parent)
        self.filters=filters
        self.initUI()

    def initUI(self):
        self.setGeometry(400, 0, 320, 185)
        self.setWindowTitle('Photometry Definition Widget')

        fnlab = QLabel(self)
        fnlab.setText('Filename:')
        fnlab.setGeometry(QRect(20,5,75,50))
        self.Fnbox = QLineEdit(self)
        self.Fnbox.setGeometry(QRect(85,20,200,20))
        self.Fnbox.setText('my_labels.dat')

        mlab = QLabel(self)
        mlab.setText('Mag. or Flux?:')
        mlab.setGeometry(QRect(20,40,75,50))
        self.combo_sys = QComboBox(self)
        self.combo_sys.setGeometry(QRect(105,40,100,50))
        self.combo_sys.setObjectName('morf')
        for k in ['mag','fnu']: self.combo_sys.addItem(k)

        ulab = QLabel(self)
        ulab.setText('Units?:')
        ulab.setGeometry(QRect(20,75,75,50))
        self.combo_unit = QComboBox(self)
        self.combo_unit.setGeometry(QRect(105,75,150,50))
        self.combo_unit.setObjectName('unit')
        for k in ['AB (mag)','Zfourge (fnu)','Jy (fnu)','muJy (fnu)']:
            self.combo_unit.addItem(k)

        gobutt = QPushButton('Generate!',self)
        gobutt.setToolTip('Generate the blank photometry definition file')
        gobutt.setGeometry(QRect(20,118,280,50))
        gobutt.clicked.connect(self.generate)

    def generate(self):
        lf_out = open(f'{self.Fnbox.text()}','w')
        lf_out.write('# LRsed photometry definition file\n')
        lf_out.write(f'# photometric units: {self.combo_unit.currentText()}\n')
        lf_out.write('#\n')
        lf_out.write('# Instructions:\n')
        lf_out.write('# You should create a text file containing photometry and associated errors for each band you\'d\n')
        lf_out.write('# like to use in your fitting. Each row corresponds to a single object, simply separate values with spaces.\n')
        lf_out.write('# You will also need to include a column with an object ID and a redshift. Once this is done you can then\n')
        lf_out.write(f'# replace the X\'s below with columns corresponding to {self.combo_sys.currentText()} and error for each filter\n')
        lf_out.write('# with indexing beginning at 0 (it\'s python people!).\n')
        lf_out.write('# To ignore filters simply remove the corresponding line from this file, it will not be considered in the SED fit.\n')
        lf_out.write('# (note, error column not used for ID or redshift, however a placeholder is required)\n')
        lf_out.write('# FILTER_NAME  FLUX_COLUMN  ERROR_COLUMN\n')
        lf_out.write('#\n')
        lf_out.write('ID         X   X  \n')
        lf_out.write('redshift   X   X  \n')
        for k in self.filters.keys():
            lf_out.write(f'{k}       X   X\n')
        lf_out.close()

        msg = QMessageBox()
        msg.setWindowTitle('Photo Def Generated')
        msg.setText(f'Photometric definition file {self.Fnbox.text()} generated. Open this file to find instructions on setting it up for use with LRsed')
        msg.exec_()
        
        self.close()

class startup_gui(QMainWindow):
    def __init__(self, parent=None):
        super(startup_gui, self).__init__(parent)

        self.filters = False
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
        
        self.setGeometry(0, 0, 320, 185)
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
        query = QMessageBox()
        query.setWindowTitle('New Photometry Definition?')
        query.setText('Do you need to generate a photometry definition file?')
        query.setIcon(QMessageBox.Question)
        query.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
        res = query.exec_()

        if res == QMessageBox.Yes:
            if not self.filters:
                self.load_filters()
            self.phdw = phd_window(self.filters)
            self.phdw.show()
        else:
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

        self.SEDoutputs = {}
        for k in self.photometry.data.keys():
            tmf = fit_one_SED(k,
                    self.bdc,
                    self.filters,
                    self.BPASS_folder,
                    self.photometry,
                    Rv=self.Rv
            )
            self.SEDoutputs[k] = tmf
            
        self.sedw = sed_window(self.SEDoutputs)
        self.sedw.show()
        
        
def main():
    app = QApplication(sys.argv)
    main = startup_gui()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
