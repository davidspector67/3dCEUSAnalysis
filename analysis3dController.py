from analysis3dGUI import Contrast3dAnalysisGUI
import os
import utils as ut
import lognormalFunctions as lf
from itertools import chain
import platform
from ticEditor import TICEditorGUI
import scipy.interpolate as interpolate
from PyQt5.QtCore import Qt, QLine
from PyQt5.QtGui import QPixmap, QImage, QPainter, QBitmap, QColor
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
from scipy.spatial import ConvexHull
import pyvista as pv
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QFileDialog
import shutil


# When running for testing, use the below code to silence invalid value in log warning.
# NOTE: vtkDelaunay3D warning is dealt with in the voi3dInterpolation function
# The +[CATransaction synchronize] error is a pyqt fault but has not shown to increase 
# the risk of user-end failures

system = platform.system()

class Contrast3dAnalysisController(Contrast3dAnalysisGUI):
    def __init__(self):
        super().__init__()

        self.chooseFileButton.clicked.connect(self.getTextInput) #get the text input of nifti path
        self.clearFileButton.clicked.connect(self.clearInputFilePath) #clear nifti file path

        self.inputXmlFileLocation = ""
        self.outputNiftiFileLocation = ""
        self.chooseInputFolderButton.clicked.connect(self.getInputFolder)
        self.chooseOutputFolderButton.clicked.connect(self.getOutputFolder)
        self.clearInputFolderButton.clicked.connect(self.clearInputBrowsing)
        self.clearOutputFolderButton.clicked.connect(self.clearOutputBrowsing)

        self.convertXmlReady = False
        self.calculateParamapButton.setHidden(True)


    # FIRST STEP: get the input text of the nifti input line edit
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
    # if file path is entered, check that it exists and get it, otherwise, let user
    # browse themselves for a valid nifti file; also makes directories to store drawings and masks made
    def getTextInput(self) :
        fileName, _ = QFileDialog.getOpenFileName(None, 'Open File', filter = '*.nii *.nii.gz')
        if fileName != '':
            self.niftiLineEdit.setText(fileName)
            self.inputTextPath = self.niftiLineEdit.text()
            self.feedbackText.setText("Valid file chosen - Press 'Open Image'")
            self.openIMGButton.clicked.connect(self.openInitialImageSlices) #opens initial slices

    def getOutputFolder(self): 
        fname = QFileDialog.getExistingDirectory(None, 'Select Directory')
        if fname == '':
            return
        else:
            self.niftiDestinationLineEdit.setText(fname)
            self.outputNiftiFileLocation = fname
            self.feedbackText.setText("Valid destination folder chosen")

        if self.inputXmlFileLocation:
            self.convertXmlButton.clicked.connect(self.convertXmltoNifti)
            self.feedbackText.setText("Valid folders chosen - Ready to convert to Nifti")
    
    def getInputFolder(self):
        fname = QFileDialog.getExistingDirectory(None, 'Select Directory')
        if fname != '':
            self.xmlLineEdit.setText(fname)
            self.inputXmlFileLocation = fname
            self.feedbackText.setText("Valid Folder Chosen")
            if self.outputNiftiFileLocation:
                self.convertXmlButton.clicked.connect(self.convertXmltoNifti)
                self.convertXmlReady = True

    def clearInputBrowsing(self):
        self.xmlLineEdit.setText("")
        self.inputXmlFileLocation = ""
        if self.convertXmlReady:
            self.convertXmlButton.clicked.disconnect()
            self.convertXmlReady = False
    
    def clearOutputBrowsing(self):
        self.niftiDestinationLineEdit.setText("")
        self.outputNiftiFileLocation = ""
        if self.convertXmlReady:
            self.convertXmlButton.clicked.disconnect()
            self.convertXmlReady = False


    def clearInputFilePath(self):
        self.niftiLineEdit.clear()
        if self.convertXmlReady:
            self.openIMGButton.clicked.disconnect()
            self.convertXmlReady = False



# SECOND STEP: does a double check for valid file type and then displays the initial slices into QLabels
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
    # displays all the initial 2D slices for axial, sag and coronal planes
    # sets the max of scroll bars and gets dimensions of each slice
    def openInitialImageSlices(self):
        if os.path.exists("niftiROIs"):
            shutil.rmtree("niftiROIs") #will also remove all the drawings made
        os.mkdir("niftiROIs") #for drawings

        self.feedbackText.setText("Use 'Draw ROI' to start drawing ROIs to construct a VOI")
        
        self.nibImg = nib.load(self.inputTextPath, mmap=False)
        self.dataNibImg = self.nibImg.get_fdata()
        self.dataNibImg = self.dataNibImg.astype(np.uint8)
        self.windowsComputed = False

        self.OGData4dImg = self.dataNibImg.copy()

        self.data4dImg = self.dataNibImg
        self.x, self.y, self.z, self.numSlices = self.data4dImg.shape
        self.maskCoverImg = np.zeros([self.x, self.y, self.z,4])
        self.slicesChanger.setMaximum(self.numSlices-1)
        self.curSlices.setText(str(self.curSlice+1))
        self.totalSlices.setText(str(self.numSlices))
        self.slicesChanger.valueChanged.connect(self.sliceValueChanged)
        self.slicesChanger.setDisabled(False)

        self.x -= 1
        self.y -= 1
        self.z -= 1

        self.sliceArray = np.array(list(range(self.numSlices)))

        self.totalFramesAx.setText(str(self.z+1))
        self.totalFramesSag.setText(str(self.x+1))
        self.totalFramesCor.setText(str(self.y+1))

        self.currentFrameAx.setText("1")
        self.currentFrameSag.setText("1")
        self.currentFrameCor.setText("1")

        tempAx = self.maskCoverImg[:,:,0,:] #2D data for axial
        tempAx = np.flipud(tempAx) #flipud
        tempAx = np.rot90(tempAx,3) #rotate ccw 270
        tempAx = np.require(tempAx,np.uint8, 'C')

        tempSag = self.maskCoverImg[0,:,:,:] #2D data for sagittal
        tempSag = np.flipud(tempSag) #flipud
        tempSag = np.rot90(tempSag,2) #rotate ccw 180
        tempSag = np.fliplr(tempSag)
        tempSag = np.require(tempSag,np.uint8,'C')

        tempCor = self.maskCoverImg[:,0,:,:] #2D data for coronal
        tempCor = np.rot90(tempCor,1) #rotate ccw 90
        tempCor = np.flipud(tempCor) #flipud
        tempCor = np.require(tempCor,np.uint8,'C')

        self.maskAxH, self.maskAxW = tempAx[:,:,0].shape #getting height and width for each plane
        self.maskSagH, self.maskSagW = tempSag[:,:,0].shape
        self.maskCorH, self.maskCorW = tempCor[:,:,0].shape

        self.maskBytesLineAx, _ = tempAx[:,:,0].strides #in order to create proper QImage, need to know bytes/line
        self.maskBytesLineSag, _ = tempSag[:,:,0].strides
        self.maskBytesLineCor, _ = tempCor[:,:,0].strides

        self.curMaskAxIm = QImage(tempAx, self.maskAxW, self.maskAxH, self.maskBytesLineAx, QImage.Format_ARGB32) #creating QImage
        self.curMaskSagIm = QImage(tempSag, self.maskSagW, self.maskSagH, self.maskBytesLineSag, QImage.Format_ARGB32)
        self.curMaskCorIm = QImage(tempCor, self.maskCorW, self.maskCorH, self.maskBytesLineCor, QImage.Format_ARGB32)

        self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(331,311)) #displaying QPixmap in the QLabels
        self.maskLayerSag.setPixmap(QPixmap.fromImage(self.curMaskSagIm).scaled(331,311))
        self.maskLayerCor.setPixmap(QPixmap.fromImage(self.curMaskCorIm).scaled(331,311))
        self.maskLayerAx.setMouseTracking(True)
        self.maskLayerSag.setMouseTracking(True)
        self.maskLayerCor.setMouseTracking(True)


        #create a list for frames in combobox for easy navigation
        self.listAxial = []
        for n in range(self.z):
            self.listAxial.append("Slice " + str(n+1))
        self.listSag = []
        for m in range(self.x):
            self.listSag.append("Slice " + str(m+1))
        self.listCor = []
        for l in range(self.y):
            self.listCor.append("Slice " + str(l+1))

        self.drawPolygonButton.setCheckable(True)

        #getting initial image data for axial, sag, coronal slices
        self.data2dAx = self.data4dImg[:,:,0, self.curSlice] #2D data for axial
        self.data2dAx = np.flipud(self.data2dAx) #flipud
        self.data2dAx = np.rot90(self.data2dAx,3) #rotate ccw 270
        self.data2dAx = np.require(self.data2dAx,np.uint8, 'C')

        self.data2dSag = self.data4dImg[0,:,:, self.curSlice] #2D data for sagittal
        self.data2dSag = np.flipud(self.data2dSag) #flipud
        self.data2dSag = np.rot90(self.data2dSag,2) #rotate ccw 180
        self.data2dSag = np.fliplr(self.data2dSag)
        self.data2dSag = np.require(self.data2dSag,np.uint8,'C')

        self.data2dCor = self.data4dImg[:,0,:, self.curSlice] #2D data for coronal
        self.data2dCor = np.rot90(self.data2dCor,1) #rotate ccw 90
        self.data2dCor = np.flipud(self.data2dCor) #flipud
        self.data2dCor = np.require(self.data2dCor,np.uint8,'C')

        self.heightAx, self.widthAx = self.data2dAx.shape #getting height and width for each plane
        self.heightSag, self.widthSag = self.data2dSag.shape
        self.heightCor, self.widthCor = self.data2dCor.shape

        self.bytesLineAx, _ = self.data2dAx.strides #in order to create proper QImage, need to know bytes/line
        self.bytesLineSag, _ = self.data2dSag.strides
        self.bytesLineCor, _ = self.data2dCor.strides

        self.qImgAx = QImage(self.data2dAx, self.widthAx, self.heightAx, self.bytesLineAx, QImage.Format_Grayscale8) #creating QImage
        self.qImgSag = QImage(self.data2dSag, self.widthSag, self.heightSag, self.bytesLineSag, QImage.Format_Grayscale8)
        self.qImgCor = QImage(self.data2dCor, self.widthCor, self.heightCor, self.bytesLineCor, QImage.Format_Grayscale8)

        self.pixmapAx = QPixmap.fromImage(self.qImgAx).scaled(331,311) #creating QPixmap from QImage
        self.pixmapSag = QPixmap.fromImage(self.qImgSag).scaled(331,311)
        self.pixmapCor = QPixmap.fromImage(self.qImgCor).scaled(331,311)

        self.axialPlane.setPixmap(self.pixmapAx) #displaying QPixmap in the QLabels
        self.sagPlane.setPixmap(self.pixmapSag)
        self.corPlane.setPixmap(self.pixmapCor)

        self.scrolling = True
        self.axCoverLabel.setCursor(Qt.BlankCursor)
        self.sagCoverLabel.setCursor(Qt.BlankCursor)
        self.corCoverLabel.setCursor(Qt.BlankCursor)
        self.newXVal = 0
        self.newYVal = 0
        self.newZVal = 0

        self.curAlpha.setDisabled(False)
        self.curAlpha.valueChanged.connect(self.alphaValueChanged)
        self.curAlpha.setValue(255)

        self.expandAxialButton.clicked.connect(self.enlargeAxImg) #expands certain slices
        self.expandSagButton.clicked.connect(self.enlargeSagImg)
        self.expandCorButton.clicked.connect(self.enlargeCorImg)        
        self.closeExpandButton.clicked.connect(self.closeExpandedImg)

        self.acceptPolygonButton.clicked.connect(self.acceptPolygon) #called to exit the paint function
        self.undoLastPtButton.clicked.connect(self.undoLastPoint) #deletes last drawn rectangle if on sag or cor slices

        self.undoLastROIButton.clicked.connect(self.undoLastROI)
        self.drawPolygonButton.clicked.connect(self.startROIDraw)

        self.interpolateVOIButton.clicked.connect(lambda:  self.feedbackText.setText("Must have have at least 1 ROI per plane to generate VOI"))


# THIRD STEP: if incorrect file path, make easier for user to clear input
#----------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------
    #to clear the file path to input a new one as well as clear all the directories created for saving drawings
    #of ROIs and masks of that particular image

    def sliceValueChanged(self):
        self.curSlice = int(self.sliceArray[self.slicesChanger.value()])
        self.curSlices.setText(str(self.curSlice+1))
        self.changeAxialSlices()
        self.changeSagSlices()
        self.changeCorSlices()

    def alphaValueChanged(self):
        self.alphaTracker.setValue(int(self.curAlpha.value()))
        for i in range(len(self.pointsPlotted)):
            self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2],3] = int(self.curAlpha.value())
        self.changeAxialSlices()
        self.changeSagSlices()
        self.changeCorSlices()


# FOURTH STEP: if correct file path and type, be able to scroll through all slices
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def changeAxialSlices(self):

        #if image is not expanded
        if self.newZVal > 98 and self.axialTextLabel.isHidden() is False:
            self.currentFrameAx.move(762,340)

        #if image is expanded
        if self.newZVal > 98 and self.axialTextLabel.isHidden() : 
            self.currentFrameAx.move(1040, 680) # David --> was 1032

        self.currentFrameAx.setText(str(self.newZVal+1))

        self.data2dAx = self.data4dImg[:,:,self.newZVal, self.curSlice]#, self.curSlice] #defining 2D data for axial
        self.data2dAx = np.flipud(self.data2dAx) #flipud
        self.data2dAx = np.rot90(self.data2dAx,3) #rotate
        self.data2dAx = np.require(self.data2dAx,np.uint8,'C')

        self.bytesLineAx, _ = self.data2dAx.strides
        self.qImgAx = QImage(self.data2dAx,self.widthAx, self.heightAx, self.bytesLineAx, QImage.Format_Grayscale8)

        tempAx = self.maskCoverImg[:,:,self.newZVal,:] #2D data for axial
        tempAx = np.flipud(tempAx) #flipud
        tempAx = np.rot90(tempAx,3) #rotate ccw 270
        tempAx = np.require(tempAx,np.uint8, 'C')

        self.curMaskAxIm = QImage(tempAx, self.maskAxW, self.maskAxH, self.maskBytesLineAx, QImage.Format_ARGB32) #creating QImage

        if self.id == 0:
            self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(331,311)) #displaying QPixmap in the QLabels
            self.axialPlane.setPixmap(QPixmap.fromImage(self.qImgAx).scaled(331,311)) #otherwise, would just display the normal unmodified q_img
        elif self.id == 1:
            self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(680, 638))
            self.axialPlane.setPixmap(QPixmap.fromImage(self.qImgAx).scaled(680,638)) #otherwise, would just display the normal unmodified q_img


    def changeSagSlices(self):

        #if image is not expanded
        if self.newXVal > 98 and self.axialTextLabel.isHidden() is False:
            self.currentFrameSag.move(1092,340)

        #if image is expanded
        if self.newXVal > 98 and self.axialTextLabel.isHidden() :
            self.currentFrameSag.move(1032, 680)

        self.currentFrameSag.setText(str(self.newXVal+1))

        self.data2dSag = self.data4dImg[self.newXVal,:,:, self.curSlice]#, self.curSlice]
        self.data2dSag = np.flipud(self.data2dSag) #flipud
        self.data2dSag = np.rot90(self.data2dSag,2) #rotate
        self.data2dSag = np.fliplr(self.data2dSag)
        self.data2dSag = np.require(self.data2dSag,np.uint8,'C')

        self.bytesLineSag, _ = self.data2dSag.strides
        self.qImgSag = QImage(self.data2dSag,self.widthSag, self.heightSag, self.bytesLineSag, QImage.Format_Grayscale8)

        tempSag = self.maskCoverImg[self.newXVal,:,:,:] #2D data for sagittal
        tempSag = np.flipud(tempSag) #flipud
        tempSag = np.rot90(tempSag,2) #rotate ccw 180
        tempSag = np.fliplr(tempSag)
        tempSag = np.require(tempSag,np.uint8,'C')
        
        self.curMaskSagIm = QImage(tempSag, self.maskSagW, self.maskSagH, self.maskBytesLineSag, QImage.Format_ARGB32)

        if self.id == 0:
            self.maskLayerSag.setPixmap(QPixmap.fromImage(self.curMaskSagIm).scaled(331,311))
            self.sagPlane.setPixmap(QPixmap.fromImage(self.qImgSag).scaled(331,311))
        elif self.id == 2:
            self.maskLayerSag.setPixmap(QPixmap.fromImage(self.curMaskSagIm).scaled(680, 638))
            self.sagPlane.setPixmap(QPixmap.fromImage(self.qImgSag).scaled(680,638))


    def changeCorSlices(self):

        #if image is not expanded
        if self.newYVal > 98 and self.axialTextLabel.isHidden() is False:
            self.currentFrameCor.move(1090,700)

        #if image is expanded
        if self.newYVal > 98 and self.axialTextLabel.isHidden() :
            self.currentFrameCor.move(1022, 680)

        self.currentFrameCor.setText(str(self.newYVal+1))

        self.data2dCor = self.data4dImg[:,self.newYVal,:, self.curSlice]#, self.curSlice]
        self.data2dCor = np.rot90(self.data2dCor,1) #rotate
        self.data2dCor = np.flipud(self.data2dCor) #flipud
        self.data2dCor = np.require(self.data2dCor, np.uint8,'C')

        self.bytesLineCor, _ = self.data2dCor.strides
        self.qImgCor = QImage(self.data2dCor,self.widthCor,self.heightCor, self.bytesLineCor, QImage.Format_Grayscale8)

        tempCor = self.maskCoverImg[:,self.newYVal,:,:] #2D data for coronal
        tempCor = np.rot90(tempCor,1) #rotate ccw 90
        tempCor = np.flipud(tempCor) #flipud
        tempCor = np.require(tempCor,np.uint8,'C')

        self.curMaskCorIm = QImage(tempCor, self.maskCorW, self.maskCorH, self.maskBytesLineCor, QImage.Format_ARGB32)

        if self.id == 0:
            self.maskLayerCor.setPixmap(QPixmap.fromImage(self.curMaskCorIm).scaled(331,311))
            self.corPlane.setPixmap(QPixmap.fromImage(self.qImgCor).scaled(331,311))

        elif self.id == 3:
            self.maskLayerCor.setPixmap(QPixmap.fromImage(self.curMaskCorIm).scaled(680, 638))
            self.corPlane.setPixmap(QPixmap.fromImage(self.qImgCor).scaled(680,638))


# FIFTH STEP: if want to enlarge images of certain cut, can do so
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
    def enlargeAxImg(self):
        if self.id != 0:
            self.closeExpandedImg()
        if self.id != 1:
            self.id = 1 #to help painter identify which plane was enlarged
            self.axialPlane.setHidden(False)
            self.ofTextAx.setHidden(False)
            self.currentFrameAx.setHidden(False)
            self.totalFramesAx.setHidden(False)
            # self.save_seg_Axial.setHidden(False)
            self.alphaTracker.setHidden(True)
            self.curAlpha.setHidden(True)
            self.curSlices.setHidden(True)
            self.totalSlices.setHidden(True)
            self.slicesChanger.setHidden(True)
            self.legend.setHidden(True)
            self.slicesLabel.setHidden(True)
            self.alphaLabel.setHidden(True)
            self.slicesOf.setHidden(True)
            self.axialPlane.move(470,30)
            self.axialPlane.resize(680,638)
            self.axCoverLabel.move(470, 30)
            self.axCoverLabel.resize(680, 638)
            self.maskLayerAx.move(470, 30)
            self.maskLayerAx.resize(680, 638)
            tempAx = self.maskCoverImg[:,:,self.newZVal,:] #2D data for axial
            tempAx = np.flipud(tempAx) #flipud
            tempAx = np.rot90(tempAx,3) #rotate ccw 270
            tempAx = np.require(tempAx,np.uint8, 'C')
            self.curMaskAxIm = QImage(tempAx, self.maskAxW, self.maskAxH, self.maskBytesLineAx, QImage.Format_ARGB32) #creating QImage
            self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(680,638)) #displaying QPixmap in the QLabels
            self.maskLayerSag.setHidden(True)
            self.maskLayerCor.setHidden(True)
            self.axCoverPixmap = QPixmap(680, 638)
            self.axCoverPixmap.fill(Qt.transparent)
            self.axCoverLabel.setPixmap(self.axCoverPixmap)
            self.sagCoverLabel.setHidden(True)
            self.corCoverLabel.setHidden(True)
            self.compressLabel.setHidden(True)
            self.compressValue.setHidden(True)
            self.windowHeightLabel.setHidden(True)
            self.windowDepthLabel.setHidden(True)
            self.windowWidthLabel.setHidden(True)
            self.windowHeightValue.setHidden(True)
            self.windowDepthValue.setHidden(True)
            self.windowWidthValue.setHidden(True)
            self.legend.setHidden(True)

            self.sagPlane.setHidden(True) #hide other labels which would get in way of displaying expanded image
            self.corPlane.setHidden(True)

            self.axialTextLabel.setHidden(True)
            self.sagTextLabel.setHidden(True)
            self.corTextLabel.setHidden(True)

            self.ofTextSag.setHidden(True)
            self.ofTextCor.setHidden(True)
            self.currentFrameSag.setHidden(True)
            self.currentFrameCor.setHidden(True)
            self.totalFramesSag.setHidden(True)
            self.totalFramesCor.setHidden(True)

            self.currentFrameAx.move(1040, 680)
            self.ofTextAx.move(1060, 680)
            self.totalFramesAx.move(1075, 680)

            self.axialPlane.setPixmap(QPixmap.fromImage(self.qImgAx).scaled(680,638))

            tempAx = self.maskCoverImg[:,:,self.newZVal,:] #2D data for axial
            tempAx = np.flipud(tempAx) #flipud
            tempAx = np.rot90(tempAx,3) #rotate ccw 270  
            tempAx = np.require(tempAx,np.uint8, 'C')

            self.curMaskAxIm = QImage(tempAx, self.maskAxW, self.maskAxH, self.maskBytesLineAx, QImage.Format_ARGB32) #creating QImage
            self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(680,638)) #displaying QPixmap in the QLabels
            self.axCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.axCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            axVertLine = QLine(int(self.newXVal/self.x*680), 0, int(self.newXVal/self.x*680), 638)
            axLatLine = QLine(0, int(self.newYVal/self.y*638), 680, int(self.newYVal/self.y*638))
            painter.drawLines([axVertLine, axLatLine])
            painter.end()

    def enlargeSagImg(self):
        if self.id != 0:
            self.closeExpandedImg()
        if self.id != 2:
            self.sagCoverLabel.move(470, 30)
            self.sagCoverLabel.resize(680, 638)
            self.sagCoverPixmap = QPixmap(680, 638)
            self.sagCoverPixmap.fill(Qt.transparent)
            self.sagCoverLabel.setPixmap(self.sagCoverPixmap)
            self.axCoverLabel.setHidden(True)
            self.corCoverLabel.setHidden(True)
            self.alphaTracker.setHidden(True)
            self.curAlpha.setHidden(True)
            self.curSlices.setHidden(True)
            self.totalSlices.setHidden(True)
            self.slicesChanger.setHidden(True)
            self.slicesLabel.setHidden(True)
            self.alphaLabel.setHidden(True)
            self.slicesOf.setHidden(True)
            self.legend.setHidden(True)
            self.id = 2
            self.sagPlane.setHidden(False)
            self.ofTextSag.setHidden(False)
            self.currentFrameSag.setHidden(False)
            self.totalFramesSag.setHidden(False)
            self.sagPlane.move(470,30)
            self.sagPlane.resize(680,638)
            self.maskLayerSag.move(470, 30)
            self.maskLayerSag.resize(680, 638)
            self.maskLayerAx.setHidden(True)
            self.maskLayerCor.setHidden(True)
            self.compressLabel.setHidden(True)
            self.compressValue.setHidden(True)
            self.windowHeightLabel.setHidden(True)
            self.windowDepthLabel.setHidden(True)
            self.windowWidthLabel.setHidden(True)
            self.windowHeightValue.setHidden(True)
            self.windowDepthValue.setHidden(True)
            self.windowWidthValue.setHidden(True)
            self.legend.setHidden(True)

            self.axialPlane.setHidden(True)
            self.corPlane.setHidden(True)

            self.axialTextLabel.setHidden(True)
            self.sagTextLabel.setHidden(True)
            self.corTextLabel.setHidden(True)

            self.ofTextAx.setHidden(True)
            self.ofTextCor.setHidden(True)
            self.currentFrameAx.setHidden(True)
            self.currentFrameCor.setHidden(True)
            self.totalFramesAx.setHidden(True)
            self.totalFramesCor.setHidden(True)

            self.currentFrameSag.move(1040, 680)
            self.ofTextSag.move(1060,680)
            self.totalFramesSag.move(1075,680)

            self.sagPlane.setPixmap(QPixmap.fromImage(self.qImgSag).scaled(680,638)) #otherwise, would just display the normal unmodified q_img

            self.sagCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.sagCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            sagVertLine = QLine(int(self.newZVal/self.z*680), 0, int(self.newZVal/self.z*680), 638)
            sagLatLine = QLine(0, int(self.newYVal/self.y*638), 680, int(self.newYVal/self.y*638))
            painter.drawLines([sagVertLine, sagLatLine])
            painter.end()
            
            tempSag = self.maskCoverImg[self.newXVal,:,:,:] #2D data for sagittal
            tempSag = np.flipud(tempSag) #flipud
            tempSag = np.rot90(tempSag,2) #rotate ccw 180
            tempSag = np.fliplr(tempSag)
            tempSag = np.require(tempSag,np.uint8,'C')

            self.curMaskSagIm = QImage(tempSag, self.widthSag, self.heightSag, self.maskBytesLineSag, QImage.Format_ARGB32)
            self.maskLayerSag.setPixmap(QPixmap.fromImage(self.curMaskSagIm).scaled(680, 638))
       

    def enlargeCorImg(self):
        if self.id != 0:
            self.closeExpandedImg()
        if self.id != 3:
            self.corCoverLabel.move(470, 30)
            self.corCoverLabel.resize(680, 638)
            self.corCoverPixmap = QPixmap(680, 638)
            self.corCoverPixmap.fill(Qt.transparent)
            self.corCoverLabel.setPixmap(self.corCoverPixmap)
            self.axCoverLabel.setHidden(True)
            self.sagCoverLabel.setHidden(True)
            self.alphaTracker.setHidden(True)
            self.curAlpha.setHidden(True)
            self.curSlices.setHidden(True)
            self.totalSlices.setHidden(True)
            self.slicesChanger.setHidden(True)
            self.slicesLabel.setHidden(True)
            self.alphaLabel.setHidden(True)
            self.slicesOf.setHidden(True)
            self.legend.setHidden(True)
            self.id = 3
            self.corPlane.setHidden(False)
            self.ofTextCor.setHidden(False)
            self.currentFrameCor.setHidden(False)
            self.totalFramesCor.setHidden(False)
            self.corPlane.move(470,30)
            self.corPlane.resize(680,638)
            self.maskLayerCor.move(470, 30)
            self.maskLayerCor.resize(680, 638)
            self.maskLayerAx.setHidden(True)
            self.maskLayerSag.setHidden(True)
            self.compressLabel.setHidden(True)
            self.compressValue.setHidden(True)
            self.windowHeightLabel.setHidden(True)
            self.windowDepthLabel.setHidden(True)
            self.windowWidthLabel.setHidden(True)
            self.windowHeightValue.setHidden(True)
            self.windowDepthValue.setHidden(True)
            self.windowWidthValue.setHidden(True)
            self.legend.setHidden(True)

            self.axialPlane.setHidden(True)
            self.sagPlane.setHidden(True)

            self.axialTextLabel.setHidden(True)
            self.sagTextLabel.setHidden(True)
            self.corTextLabel.setHidden(True)

            self.ofTextAx.setHidden(True)
            self.ofTextSag.setHidden(True)
            self.currentFrameAx.setHidden(True)
            self.currentFrameSag.setHidden(True)
            self.totalFramesAx.setHidden(True)
            self.totalFramesSag.setHidden(True)

            self.currentFrameCor.move(1040,680)
            self.ofTextCor.move(1060,680)
            self.totalFramesCor.move(1075,680)

            self.corPlane.setPixmap(QPixmap.fromImage(self.qImgCor).scaled(680,638)) #otherwise, would just display the normal unmodified q_img

            tempCor = self.maskCoverImg[:,self.newYVal,:,:] #2D data for coronal
            tempCor = np.rot90(tempCor,1) #rotate ccw 90
            tempCor = np.flipud(tempCor) #flipud
            tempCor = np.require(tempCor,np.uint8,'C')

            self.curMaskCorIm = QImage(tempCor, self.maskCorW, self.maskCorH, self.maskBytesLineCor, QImage.Format_ARGB32)
            self.maskLayerCor.setPixmap(QPixmap.fromImage(self.curMaskCorIm).scaled(680,638))

            self.corCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.corCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            corVertLine = QLine(int(self.newXVal/self.x*680), 0, int(self.newXVal/self.x*680), 638)
            corLatLine = QLine(0, int(self.newZVal/self.z*638), 680, int(self.newZVal/self.z*638))
            painter.drawLines([corVertLine, corLatLine])
            painter.end()      


# SIXTH STEP: when close expanded image, return to normal layout and functions
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
    def closeExpandedImg(self):
        if self.id != 0:
            if self.id == 1:
                self.axCoverLabel.move(470, 30)
                self.axCoverLabel.resize(331, 311)
                self.axCoverPixmap = QPixmap(331, 311)
                self.axCoverPixmap.fill(Qt.transparent)
                self.axCoverLabel.setPixmap(self.axCoverPixmap)
                self.xCur = int(self.newXVal*331/(self.widthAx-1)) + 479
                self.yCur = int(self.newYVal*311/(self.heightAx-1)) + 30
                self.maskLayerSag.setHidden(False)
                self.maskLayerCor.setHidden(False)
                self.maskLayerAx.move(470, 30)
                self.maskLayerAx.resize(331, 311)
                tempAx = self.maskCoverImg[:,:,self.newZVal,:] #2D data for axial
                tempAx = np.flipud(tempAx) #flipud
                tempAx = np.rot90(tempAx,3) #rotate ccw 270  
                tempAx = np.require(tempAx,np.uint8, 'C')
                self.curMaskAxIm = QImage(tempAx, self.maskAxW, self.maskAxH, self.maskBytesLineAx, QImage.Format_ARGB32) #creating QImage
                self.maskLayerAx.setPixmap(QPixmap.fromImage(self.curMaskAxIm).scaled(331,311)) #displaying QPixmap in the QLabels
            elif self.id == 2:
                self.sagCoverLabel.move(820, 30)
                self.sagCoverLabel.resize(331, 311)
                self.sagCoverPixmap = QPixmap(331, 311)
                self.sagCoverPixmap.fill(Qt.transparent)
                self.sagCoverLabel.setPixmap(self.sagCoverPixmap)
                self.xCur = int(self.newZVal*331/(self.widthSag-1)) + 820
                self.yCur = int(self.newYVal*311/(self.heightSag-1)) + 30
                self.maskLayerAx.setHidden(False)
                self.maskLayerCor.setHidden(False)
                self.maskLayerSag.move(820, 30)
                self.maskLayerSag.resize(331, 311)
                tempSag = self.maskCoverImg[self.newXVal,:,:,:] #2D data for sagittal
                tempSag = np.flipud(tempSag) #flipud
                tempSag = np.rot90(tempSag,2) #rotate ccw 180
                tempSag = np.fliplr(tempSag)
                tempSag = np.require(tempSag,np.uint8,'C')
                self.curMaskSagIm = QImage(tempSag, self.widthSag, self.heightSag, self.maskBytesLineSag, QImage.Format_ARGB32)
                self.maskLayerSag.setPixmap(QPixmap.fromImage(self.curMaskSagIm).scaled(331, 311))
            elif self.id == 3:
                self.corCoverLabel.move(820, 390)
                self.corCoverLabel.resize(331, 311)
                self.corCoverPixmap = QPixmap(331, 311)
                self.corCoverPixmap.fill(Qt.transparent)
                self.corCoverLabel.setPixmap(self.corCoverPixmap)
                self.maskLayerAx.setHidden(False)
                self.maskLayerSag.setHidden(False)
                self.maskLayerCor.move(820, 390)
                self.maskLayerCor.resize(331, 311)
                tempCor = self.maskCoverImg[:,self.newYVal,:,:] #2D data for coronal
                tempCor = np.rot90(tempCor,1) #rotate ccw 90
                tempCor = np.flipud(tempCor) #flipud
                tempCor = np.require(tempCor,np.uint8,'C')
                self.curMaskCorIm = QImage(tempCor, self.maskCorW, self.maskCorH, self.maskBytesLineCor, QImage.Format_ARGB32)
                self.maskLayerCor.setPixmap(QPixmap.fromImage(self.curMaskCorIm).scaled(331,311))
            self.id = 0

            self.axCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.axCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            axVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            axLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([axVertLine, axLatLine])
            painter.end()
            self.sagCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.sagCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            sagVertLine = QLine(int(self.newZVal/self.z*331), 0, int(self.newZVal/self.z*331), 311)
            sagLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([sagVertLine, sagLatLine])
            painter.end()
            self.corCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.corCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            corVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            corLatLine = QLine(0, int(self.newZVal/self.z*311), 331, int(self.newZVal/self.z*311))
            painter.drawLines([corVertLine, corLatLine])
            painter.end()
            self.update()

            if self.aucParamapButton.isChecked() or self.peParamapButton.isChecked() or self.tpParamapButton.isChecked() or self.mttParamapButton.isChecked():
                self.legend.setHidden(False)
            self.compressLabel.setHidden(False)
            self.compressValue.setHidden(False)
            self.windowHeightLabel.setHidden(False)
            self.windowDepthLabel.setHidden(False)
            self.windowWidthLabel.setHidden(False)
            self.windowHeightValue.setHidden(False)
            self.windowDepthValue.setHidden(False)
            self.windowWidthValue.setHidden(False)
            if self.ticComputed:
                self.legend.setHidden(False)
            self.axCoverLabel.setHidden(False)
            self.sagCoverLabel.setHidden(False)
            self.corCoverLabel.setHidden(False)
            self.axialPlane.move(470,30)
            self.axialPlane.resize(331,311)
            self.sagPlane.move(820,30)
            self.sagPlane.resize(331,311)
            self.corPlane.move(820,390)
            self.corPlane.resize(331,311)
            self.alphaTracker.setHidden(False)
            self.curAlpha.setHidden(False)
            self.curSlices.setHidden(False)
            self.totalSlices.setHidden(False)
            self.slicesChanger.setHidden(False)
            self.slicesLabel.setHidden(False)
            self.alphaLabel.setHidden(False)
            self.slicesOf.setHidden(False)

            self.ofTextAx.move(790,340)
            self.ofTextSag.move(1120,340)
            self.ofTextCor.move(1120,700)

            self.totalFramesAx.move(805,340)
            self.totalFramesAx.resize(31,16)

            self.totalFramesSag.move(1135,340)
            self.totalFramesSag.resize(31,16)

            self.totalFramesCor.move(1135,700)
            self.totalFramesCor.resize(31,16)

            self.currentFrameAx.move(770,340)
            self.currentFrameAx.resize(31,16)

            self.currentFrameSag.move(1100,340)
            self.currentFrameSag.resize(31,16)

            self.currentFrameCor.move(1090,700)
            self.currentFrameCor.resize(31,16)

            self.totalFramesAx.setText(str(self.z+1))
            self.totalFramesSag.setText(str(self.x+1))
            self.totalFramesCor.setText(str(self.y+1))

            self.axialPlane.setHidden(False)
            self.sagPlane.setHidden(False)
            self.corPlane.setHidden(False)

            self.axialTextLabel.setHidden(False)
            self.sagTextLabel.setHidden(False)
            self.corTextLabel.setHidden(False)

            self.ofTextAx.setHidden(False)
            self.ofTextSag.setHidden(False)
            self.ofTextCor.setHidden(False)
            self.currentFrameAx.setHidden(False)
            self.currentFrameSag.setHidden(False)
            self.currentFrameCor.setHidden(False)
            self.totalFramesAx.setHidden(False)
            self.totalFramesSag.setHidden(False)
            self.totalFramesCor.setHidden(False)

            # self.save_seg_Axial.setHidden(False)
            # self.save_seg_Sag.setHidden(False)
            # self.save_seg_Cor.setHidden(False)

            self.axialPlane.setPixmap(QPixmap.fromImage(self.qImgAx).scaled(331,311)) #otherwise, would just display the normal unmodified q_img
            self.sagPlane.setPixmap(QPixmap.fromImage(self.qImgSag).scaled(331,311)) #otherwise, would just display the normal unmodified q_img
            self.corPlane.setPixmap(QPixmap.fromImage(self.qImgCor).scaled(331,311)) #otherwise, would just display the normal unmodified q_img



# SEVENTH STEP: user must draw a region of interest in each plane in order to generate a volume of interest
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
    def paintEvent(self,event):
        scrolling = "none"
        if self.id ==0 and self.scrolling:
            if self.xCur < 811 and self.xCur > 478 and self.yCur < 342 and self.yCur > 29 and (self.painted == "none" or self.painted == "ax"):
                self.actualX = int((self.xCur - 479)*(self.widthAx-1)/331)
                self.actualY = int((self.yCur - 30)*(self.heightAx-1)/311)
                scrolling = "ax"
                self.axCoverLabel.pixmap().fill(Qt.transparent)
                painter = QPainter(self.axCoverLabel.pixmap())
                painter.setPen(Qt.yellow)
                axVertLine = QLine(self.xCur - 479, 0, self.xCur - 479, 311)
                axLatLine = QLine(0, self.yCur - 30, 331, self.yCur - 30)
                painter.drawLines([axVertLine, axLatLine])
                painter.end()
            elif self.xCur < 1152 and self.xCur > 819 and self.yCur < 342 and self.yCur > 29 and (self.painted == "none" or self.painted == "sag"):
                self.actualX = int((self.xCur-820)*(self.widthSag-1)/331)
                self.actualY = int((self.yCur-30)*(self.heightSag-1)/311)
                scrolling = "sag"
                self.sagCoverLabel.pixmap().fill(Qt.transparent)
                painter = QPainter(self.sagCoverLabel.pixmap())
                painter.setPen(Qt.yellow)
                sagVertLine = QLine(self.xCur - 820, 0, self.xCur - 820, 311)
                sagLatLine = QLine(0, self.yCur - 30, 331, self.yCur - 30)
                painter.drawLines([sagVertLine, sagLatLine])
                painter.end()
            elif self.xCur < 1152 and self.xCur > 819 and self.yCur < 702 and self.yCur > 389 and (self.painted == "none" or self.painted == "cor"):
                self.actualX = int((self.xCur-820)*(self.widthCor-1)/331)
                self.actualY = int((self.yCur-390)*(self.heightCor-1)/311)
                scrolling = "cor"
                self.corCoverLabel.pixmap().fill(Qt.transparent)
                painter = QPainter(self.corCoverLabel.pixmap())
                painter.setPen(Qt.yellow)
                corVertLine = QLine(self.xCur - 820, 0, self.xCur - 820, 311)
                corLatLine = QLine(0, self.yCur - 390, 331, self.yCur-390)
                painter.drawLines([corVertLine, corLatLine])
                painter.end()

        elif self.id != 0 and self.scrolling:
            if self.xCur < 1151 and self.xCur > 469 and self.yCur < 669 and self.yCur > 29:
                if self.id == 1:
                    self.actualX = int((self.widthAx-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightAx-1)*(self.yCur-30)/638)
                    scrolling = "ax"
                    self.axCoverLabel.pixmap().fill(Qt.transparent)
                    painter = QPainter(self.axCoverLabel.pixmap())
                    painter.setPen(Qt.yellow)
                    axVertLine = QLine(self.xCur - 470, 0, self.xCur - 470, 638)
                    axLatLine = QLine(0, self.yCur - 30, 680, self.yCur - 30)
                    painter.drawLines([axVertLine, axLatLine])
                    painter.end()
                elif self.id == 2:
                    self.actualX = int((self.widthSag-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightSag-1)*(self.yCur-30)/638)
                    scrolling = "sag"
                    self.sagCoverLabel.pixmap().fill(Qt.transparent)
                    painter = QPainter(self.sagCoverLabel.pixmap())
                    painter.setPen(Qt.yellow)
                    axVertLine = QLine(self.xCur - 470, 0, self.xCur - 470, 638)
                    axLatLine = QLine(0, self.yCur - 30, 680, self.yCur - 30)
                    painter.drawLines([axVertLine, axLatLine])
                    painter.end()
                elif self.id == 3:
                    self.actualX = int((self.widthCor-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightCor-1)*(self.yCur-30)/638)
                    scrolling = "cor"
                    self.corCoverLabel.pixmap().fill(Qt.transparent)
                    painter = QPainter(self.corCoverLabel.pixmap())
                    painter.setPen(Qt.yellow)
                    axVertLine = QLine(self.xCur - 470, 0, self.xCur - 470, 638)
                    axLatLine = QLine(0, self.yCur - 30, 680, self.yCur - 30)
                    painter.drawLines([axVertLine, axLatLine])
                    painter.end()

        if scrolling == "ax":
            self.newXVal = self.actualX
            self.newYVal = self.actualY
            self.changeSagSlices()
            self.changeCorSlices()
            self.sagCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.sagCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            sagVertLine = QLine(int(self.newZVal/self.z*331), 0, int(self.newZVal/self.z*331), 311)
            sagLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([sagVertLine, sagLatLine])
            painter.end()
            
            self.corCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.corCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            corVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            corLatLine = QLine(0, int(self.newZVal/self.z*311), 331, int(self.newZVal/self.z*311))
            painter.drawLines([corVertLine, corLatLine])
            painter.end()
            self.update()

        elif scrolling == "sag":
            self.newZVal = self.actualX
            self.newYVal = self.actualY
            self.changeAxialSlices()
            self.changeCorSlices()
            self.axCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.axCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            axVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            axLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([axVertLine, axLatLine])
            painter.end()
            
            self.corCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.corCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            corVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            corLatLine = QLine(0, int(self.newZVal/self.z*311), 331, int(self.newZVal/self.z*311))
            painter.drawLines([corVertLine, corLatLine])
            painter.end()
            self.update()

        elif scrolling == "cor":
            self.newXVal = self.actualX
            self.newZVal = self.actualY
            self.changeAxialSlices()
            self.changeSagSlices()
            self.axCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.axCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            axVertLine = QLine(int(self.newXVal/self.x*331), 0, int(self.newXVal/self.x*331), 311)
            axLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([axVertLine, axLatLine])
            painter.end()

            self.sagCoverLabel.pixmap().fill(Qt.transparent)
            painter = QPainter(self.sagCoverLabel.pixmap())
            painter.setPen(Qt.yellow)
            sagVertLine = QLine(int(self.newZVal/self.z*331), 0, int(self.newZVal/self.z*331), 311)
            sagLatLine = QLine(0, int(self.newYVal/self.y*311), 331, int(self.newYVal/self.y*311))
            painter.drawLines([sagVertLine, sagLatLine])
            painter.end()
            self.update()


#----------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------
# defines mouse events here so that paint function works
    def mousePressEvent(self,event):
        self.xCur = event.x()
        self.yCur = event.y()
        self.newPointPlotted = False
        if self.drawPolygonButton.isChecked() and self.scrolling:
            # Plot ROI points
            if self.id == 0:
                if (self.xCur < 811 and self.xCur > 478 and self.yCur < 342 and self.yCur > 29) and (self.painted == "none" or self.painted == "ax"):
                    self.actualX = int((self.xCur - 479)*(self.widthAx-1)/331)
                    self.actualY = int((self.yCur - 30)*(self.heightAx-1)/311)
                    self.maskCoverImg[self.actualX, self.actualY, self.newZVal] = [0, 0, 255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "ax"
                    self.sagCoverLabel.setCursor(Qt.ArrowCursor)
                    self.corCoverLabel.setCursor(Qt.ArrowCursor)
                    self.curROIDrawn = False
                elif (event.x() < 1152 and event.x() > 819 and event.y() < 342 and event.y() > 29) and (self.painted == "none" or self.painted == "sag"):
                    self.actualX = int((self.xCur-820)*(self.widthSag-1)/331)
                    self.actualY = int((self.yCur-30)*(self.heightSag-1)/311)
                    self.maskCoverImg[self.newXVal, self.actualY, self.actualX] = [0,0,255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "sag"
                    self.axCoverLabel.setCursor(Qt.ArrowCursor)
                    self.corCoverLabel.setCursor(Qt.ArrowCursor)
                    self.curROIDrawn = False
                elif (event.x() < 1152 and event.x() > 819 and event.y() < 702 and event.y() > 389) and (self.painted == "none" or self.painted == "cor"):
                    self.actualX = int((self.xCur-820)*(self.widthCor-1)/331)
                    self.actualY = int((self.yCur-390)*(self.heightCor-1)/311)
                    self.maskCoverImg[self.actualX, self.newYVal, self.actualY] = [0,0,255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "cor"
                    self.axCoverLabel.setCursor(Qt.ArrowCursor)
                    self.sagCoverLabel.setCursor(Qt.ArrowCursor)
                    self.curROIDrawn = False
            elif self.id != 0 and self.xCur < 1151 and self.xCur > 469 and self.yCur < 669 and self.yCur > 29:
                if self.id == 1 and self.axialPlane.isHidden() is False and self.sagPlane.isHidden()  and self.corPlane.isHidden() :
                    self.actualX = int((self.widthAx-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightAx-1)*(self.yCur-30)/638)
                    self.maskCoverImg[self.actualX, self.actualY, self.newZVal] = [0,0,255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "ax"
                elif self.id == 2 and self.sagPlane.isHidden() is False and self.corPlane.isHidden()  and self.axialPlane.isHidden() :
                    self.actualX = int((self.widthSag-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightSag-1)*(self.yCur-30)/638)
                    self.maskCoverImg[self.newXVal, self.actualY, self.actualX] = [0,0,255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "sag"
                elif self.id == 3 and self.corPlane.isHidden() is False and self.axialPlane.isHidden()  and self.sagPlane.isHidden() :
                    self.actualX = int((self.widthCor-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightCor-1)*(self.yCur-30)/638)
                    self.maskCoverImg[self.actualX, self.newYVal, self.actualY] = [0,0,255,int(self.curAlpha.value())]
                    self.curPointsPlottedX.append(self.actualX)
                    self.curPointsPlottedY.append(self.actualY)
                    self.newPointPlotted = True
                    self.painted = "cor"
            self.changeSagSlices()
            self.changeCorSlices()
            self.changeAxialSlices()


    def mouseDoubleClickEvent(self, event):
        if self.scrolling:
            # Stop cross-hair pointer updates
            self.scrolling = False
            self.axCoverLabel.setCursor(Qt.ArrowCursor)
            self.sagCoverLabel.setCursor(Qt.ArrowCursor)
            self.corCoverLabel.setCursor(Qt.ArrowCursor)
        else:
            # Re-start cross-hair pointer updates
            self.scrolling = True
            if self.painted == "ax":
                self.axCoverLabel.setCursor(Qt.BlankCursor)
            elif self.painted == "sag":
                self.sagCoverLabel.setCursor(Qt.BlankCursor)
            elif self.painted == "cor":
                self.corCoverLabel.setCursor(Qt.BlankCursor)
            else:
                self.axCoverLabel.setCursor(Qt.BlankCursor)
                self.sagCoverLabel.setCursor(Qt.BlankCursor)
                self.corCoverLabel.setCursor(Qt.BlankCursor)
            self.paintEvent(event)
        if self.drawPolygonButton.isChecked():
            # Disregard previously plotted point
            if self.id == 0:
                if self.xCur < 811 and self.xCur > 478 and self.yCur < 342 and self.yCur > 29 and (self.painted == "none" or self.painted == "ax"):
                    self.actualX = int((self.xCur - 479)*(self.widthAx-1)/331)
                    self.actualY = int((self.yCur - 30)*(self.heightAx-1)/311)
                elif self.xCur < 1152 and self.xCur > 819 and self.yCur < 342 and self.yCur > 29 and (self.painted == "none" or self.painted == "sag"):
                    self.actualX = int((self.xCur-820)*(self.widthSag-1)/331)
                    self.actualY = int((self.yCur-30)*(self.heightSag-1)/311)
                elif self.xCur < 1152 and self.xCur > 819 and self.yCur < 702 and self.yCur > 389 and (self.painted == "none" or self.painted == "cor"):
                    self.actualX = int((self.xCur-820)*(self.widthCor-1)/331)
                    self.actualY = int((self.yCur-390)*(self.heightCor-1)/311)
            elif self.id != 0 and self.xCur < 1151 and self.xCur > 469 and self.yCur < 669 and self.yCur > 29:
                if self.id == 1 and self.axialPlane.isHidden() is False and self.sagPlane.isHidden()  and self.corPlane.isHidden() :
                    self.actualX = int((self.widthAx-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightAx-1)*(self.yCur-30)/638)
                elif self.id == 2 and self.sagPlane.isHidden() is False and self.corPlane.isHidden()  and self.axialPlane.isHidden() :
                    self.actualX = int((self.widthSag-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightSag-1)*(self.yCur-30)/638)
                elif self.id == 3 and self.corPlane.isHidden() is False and self.axialPlane.isHidden()  and self.sagPlane.isHidden() :
                    self.actualX = int((self.widthCor-1)*(self.xCur-470)/680)
                    self.actualY = int((self.heightCor-1)*(self.yCur-30)/638)
            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()
            if self.newPointPlotted:
                self.undoLastPoint()
            self.newPointPlotted = False
        if len(self.curPointsPlottedX) == 0:
            self.painted = "none"


    def mouseMoveEvent(self, event):
        self.xCur = event.x()
        self.yCur = event.y()


# EIGHTH STEP: when accept polygon is clicked, the draw polygon button
# returns to normal - is unchecked
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
    def acceptPolygon(self):
        # 2d interpolation
        if len(self.curPointsPlottedX) and self.voiComputed == False:
            self.drawPolygonButton.setChecked(False)
            self.feedbackScrollBar.setHidden(True)
            self.feedbackText.setText("VOI Created Successfully")
            self.curPointsPlottedX.append(self.curPointsPlottedX[0])
            self.curPointsPlottedY.append(self.curPointsPlottedY[0])
            self.maskCoverImg.fill(0)
            x, y = calculateSpline(self.curPointsPlottedX, self.curPointsPlottedY)
            newROI = []
            for i in range(len(x)):
                if self.painted == "ax":
                    if len(newROI) == 0 or newROI[-1] != (int(x[i]), int(y[i]), self.newZVal):
                        newROI.append([int(x[i]), int(y[i]), self.newZVal])
                elif self.painted == "sag":
                    if len(newROI) == 0 or newROI[-1] != (self.newXVal, int(y[i]), int (x[i])):
                        newROI.append([self.newXVal, int(y[i]), int(x[i])])
                elif self.painted == "cor":
                    if len(newROI) == 0 or newROI[-1] != (int(x[i]), self.newYVal, int(y[i])):
                        newROI.append([int(x[i]), self.newYVal, int(y[i])])
            self.pointsPlotted.append(newROI)
            for i in range(len(self.pointsPlotted)):
                for j in range(len(self.pointsPlotted[i])):
                    self.maskCoverImg[self.pointsPlotted[i][j][0], self.pointsPlotted[i][j][1], self.pointsPlotted[i][j][2]] = [0,0,255,int(self.curAlpha.value())]
            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()
            self.curPointsPlottedX = []
            self.curPointsPlottedY = []
            self.planesDrawn.append(self.painted)
            self.painted = "none"
            self.curROIDrawn = True
            self.undoLastROIButton.setHidden(False)
            self.acceptPolygonButton.setHidden(True)
            if len(self.pointsPlotted) >= 3 and ((self.planesDrawn[0]!=self.planesDrawn[1]) and (self.planesDrawn[1]!=self.planesDrawn[2]) and (self.planesDrawn[2]!=self.planesDrawn[0])):
                self.interpolateVOIButton.clicked.disconnect()
                self.interpolateVOIButton.clicked.connect(self.voi3dInterpolation)


# NINTH STEP: undo region of interest drawings
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

    def undoLastPoint(self):
        if len(self.curPointsPlottedX) != 0:
            self.maskCoverImg[self.curPointsPlottedX[-1]]
            self.curPointsPlottedX.pop()
            self.curPointsPlottedY.pop()
            self.maskCoverImg.fill(0)
            for i in range(len(self.pointsPlotted)):
                for j in range(len(self.pointsPlotted[i])):
                    self.maskCoverImg[self.pointsPlotted[i][j][0], self.pointsPlotted[i][j][1], self.pointsPlotted[i][j][2]] = [0,0,255, int(self.curAlpha.value())]
            for i in range(len(self.curPointsPlottedX)):
                if self.painted == "ax":
                    self.maskCoverImg[int(self.curPointsPlottedX[i]), int(self.curPointsPlottedY[i]), self.newZVal] = [0,0,255,int(self.curAlpha.value())]
                elif self.painted == "sag":
                    self.maskCoverImg[self.newXVal, int(self.curPointsPlottedY[i]), int(self.curPointsPlottedX[i])] = [0,0,255,int(self.curAlpha.value())]
                elif self.painted == "cor":
                    self.maskCoverImg[int(self.curPointsPlottedX[i]), self.newYVal, int(self.curPointsPlottedY[i])] = [0,0,255,int(self.curAlpha.value())]

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()
        if len(self.curPointsPlottedX) == 0:
            self.painted == "none"


    def startROIDraw(self):
        if self.drawPolygonButton.isChecked():
            self.acceptPolygonButton.setHidden(False)
            self.undoLastROIButton.setHidden(True)
        else:
            self.acceptPolygonButton.setHidden(True)
            self.undoLastROIButton.setHidden(False)

    def undoLastROI(self):
        if not self.voiComputed and len(self.planesDrawn):

            if len(self.pointsPlotted) == 0:
                self.feedbackText.setText("Unable to remove ROI")
            else:
                self.pointsPlotted.pop()
                self.planesDrawn.pop()
                self.maskCoverImg.fill(0)
                for i in range(len(self.pointsPlotted)):
                    for j in range(len(self.pointsPlotted[i])):
                        self.maskCoverImg[self.pointsPlotted[i][j][0], self.pointsPlotted[i][j][1], self.pointsPlotted[i][j][2]] = [0,0,255,int(self.curAlpha.value())]
                self.changeAxialSlices()
                self.changeSagSlices()
                self.changeCorSlices()
            self.interpolateVOIButton.clicked.disconnect()
            self.interpolateVOIButton.clicked.connect(lambda:  self.feedbackText.setText("Must have at least 1 ROI per plane to generate VOI"))
            self.update()
        
        elif not len(self.planesDrawn):
            self.feedbackText.setText("No ROIs to remove")

    def voi3dInterpolation(self):
        if self.voiComputed == False:
            points = calculateSpline3D(list(chain.from_iterable(self.pointsPlotted)))

            self.pointsPlotted = []
            self.maskCoverImg.fill(0)
            
            for point in points:
                if max(self.data4dImg[tuple(point)]) != 0:
                    self.maskCoverImg[tuple(point)] = [0,0,255,int(self.curAlpha.value())]
                    self.pointsPlotted.append(tuple(point))
            if len(self.pointsPlotted) == 0:
                self.feedbackText.setText("VOI not in US image.\nDraw new VOI over US image")
                self.interpolateVOIButton.clicked.disconnect()
                self.interpolateVOIButton.clicked.connect(lambda:  self.feedbackText.setText("Must have at least 1 ROI per plane to generate VOI"))
                self.maskCoverImg.fill(0)
                self.changeAxialSlices()
                self.changeSagSlices()
                self.changeCorSlices()
                return
            
            mask = np.zeros((self.maskCoverImg.shape[0], self.maskCoverImg.shape[1], self.maskCoverImg.shape[2]))

            for point in self.pointsPlotted:
                mask[point] = 1
            from scipy.ndimage import binary_fill_holes
            for i in range(mask.shape[2]):
                figure = plt.figure()
                ax = figure.add_subplot()
                border = np.where(mask[:,:,i] == 1)
                if (not len(border[0])) or (max(border[0]) == min(border[0])) or (max(border[1]) == min(border[1])):
                    continue
                border = np.array(border).T
                hull = ConvexHull(border)
                vertices = border[hull.vertices]
                shape = vertices.shape
                vertices = np.reshape(np.append(vertices, vertices[0]), (shape[0]+1, shape[1]))

                # Linear interpolation of 2d convex hull
                tck, u_ = interpolate.splprep(vertices.T, s=0.0, k=1)
                splineX, splineY = np.array(interpolate.splev(np.linspace(0, 1, 1000), tck))

                mask[:,:,i] = np.zeros((mask.shape[0], mask.shape[1]))
                for j in range(len(splineX)):
                    mask[int(splineX[j]), int(splineY[j]), i] = 1
                ax.imshow(mask[:,:,i])
                filledMask = binary_fill_holes(mask[:,:,i])
                ax.imshow(filledMask)
                maskPoints = np.array(np.where(filledMask == True))
                for j in range(len(maskPoints[0])):
                    self.maskCoverImg[maskPoints[0][j], maskPoints[1][j], i] = [0,0,255,int(self.curAlpha.value())]
                    self.pointsPlotted.append((maskPoints[0][j], maskPoints[1][j], i))
            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()
            self.voiComputed = True
            self.drawPolygonButton.setCheckable(False)
            self.computeTICButton.clicked.connect(self.showTic)
            self.interpolateVOIButton.clicked.disconnect()

    def acceptTIC(self):
        ax = self.fig.add_subplot(111)
        self.ticEditor.ticY -= min(self.ticEditor.ticY)
        ax.plot(self.ticEditor.ticX[:,0], self.ticEditor.ticY)

        self.sliceArray = self.ticEditor.ticX[:,1]
        if self.slicesChanger.value() >= len(self.sliceArray):
            self.slicesChanger.setValue(len(self.sliceArray)-1)
            self.sliceValueChanged()
        self.slicesChanger.setMaximum(len(self.sliceArray)-1)

        if system == 'Windows':
            ax.set_xlabel("Time (s)", fontsize=8, labelpad=0.5)
            ax.set_ylabel("Signal Amplitude", fontsize=8, labelpad=0.5)
            ax.set_title("Time Intensity Curve (TIC)", fontsize=10, pad=1.5)
            ax.tick_params('both', pad=0.3, labelsize=7.2)
            plt.xticks(fontsize=6)
            plt.yticks(fontsize=6)
        else:
            ax.set_xlabel("Time (s)", fontsize=4, labelpad=0.5)
            ax.set_ylabel("Signal Amplitude", fontsize=4, labelpad=0.5)
            ax.set_title("Time Intensity Curve (TIC)", fontsize=5, pad=1.5)
            ax.tick_params('both', pad=0.3, labelsize=3.6)
            plt.xticks(fontsize=3)
            plt.yticks(fontsize=3)
        tmppv = np.max(self.ticEditor.ticY)
        self.ticEditor.ticY = self.ticEditor.ticY/tmppv;

        # Bunch of checks
        if np.isnan(np.sum(self.ticEditor.ticY)):
            print('STOPPED:NaNs in the VOI')
            return;
        if np.isinf(np.sum(self.ticEditor.ticY)):
            print('STOPPED:InFs in the VOI')
            return;

        # Do the fitting
        try:
            params, popt, wholecurve = lf.data_fit([self.ticEditor.ticX[:,0], self.ticEditor.ticY], tmppv);
            ax.plot(self.ticEditor.ticX[:,0], wholecurve)
        except RuntimeError:
            print('RunTimeError')
            params = np.array([np.max(self.ticEditor.ticY)*tmppv, np.trapz(self.ticEditor.ticY*tmppv, x=self.ticEditor.ticX[:,0]), self.ticEditor.ticX[:,0], [np.argmax(self.ticEditor.ticY),0], np.max(self.ticEditor.ticX[:,0])*2, 0]);
        self.ticComputed = True
        self.fig.subplots_adjust(left=0.1, right=0.97, top=0.85, bottom=0.25)
        self.canvas.draw()
        self.feedbackLabel.setHidden(True)
        self.feedbackText.setHidden(True)
        self.ticAucLabel.setHidden(False)
        self.ticVarLabel.setHidden(False)
        self.ticPeLabel.setHidden(False)
        self.ticTpLabel.setHidden(False)
        self.ticMttLabel.setHidden(False)
        self.ticTmppvLabel.setHidden(False)
        self.ticAucVal.setHidden(False)
        self.ticPeVal.setHidden(False)
        self.ticTpVal.setHidden(False)
        self.ticMttVal.setHidden(False)
        self.ticTmppvVal.setHidden(False)
        self.ticAucVal.setText(str(int(popt[0]*1000)/1000))
        self.ticPeVal.setText(str(int(params[0]*1000)/1000))
        self.ticTpVal.setText(str(int(params[2]*100)/100))
        self.ticMttVal.setText(str(int(params[3]*100)/100))
        self.ticTmppvVal.setText(str(int(tmppv*100)/100))
        self.aucParamapButton.setCheckable(True)
        self.peParamapButton.setCheckable(True)
        self.tpParamapButton.setCheckable(True)
        self.mttParamapButton.setCheckable(True)
        xlist = []
        ylist = []
        zlist = []
        for i in self.pointsPlotted:
            xlist.append(i[0])
            ylist.append(i[1])
            zlist.append(i[2])
        # self.masterParamap = ut.paramap(self.OGData4dImg, xlist, ylist, zlist, self.header[1:4], self.header[4], 'BolusLognormal', self.compressValue.value(), int(self.windowHeightValue.value()*self.header[1]), int(self.windowWidthValue.value()*self.header[2]), int(self.windowDepthValue.value()*self.header[3]))
        # self.maxAuc = 0
        # self.minAuc = 9999
        # self.maxPe = 0
        # self.minPe = 9999
        # self.maxTp = 0
        # self.minTp = 9999
        # self.maxMtt = 0
        # self.minMtt = 9999
        # for i in range(len(self.pointsPlotted)):
        #     if self.masterParamap[self.pointsPlotted[i][0], self.pointsPlotted[i][1],self.pointsPlotted[i][2]][0] > self.maxAuc:
        #         self.maxAuc = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][0]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][0] < self.minAuc:
        #         self.minAuc = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][0]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][1] > self.maxPe:
        #         self.maxPe = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][1]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][1] < self.minPe:
        #         self.minPe = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][1] 
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][2] > self.maxTp:
        #         self.maxTp = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][2]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][2] < self.minTp:
        #         self.minTp = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][2]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][3] > self.maxMtt:
        #         self.maxMtt = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][3]
        #     if self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][3] < self.minMtt:
        #         self.minMtt = self.masterParamap[self.pointsPlotted[i][0],self.pointsPlotted[i][1],self.pointsPlotted[i][2]][3]
        self.windowsComputed = True
        self.ticEditor.close()
        # if self.curSlice >= numUsedSlices:
        #     self.slicesChanger.setValue(numUsedSlices)
        #     self.sliceValueChanged()

    def showTic(self):
        if not self.windowsComputed and self.voiComputed:
            self.header = self.nibImg.header['pixdim'] # [dims, voxel dims (3 vals), timeconst, 0, 0, 0]
            times = [i*self.header[4] for i in range(1, self.OGData4dImg.shape[3]+1)]
            self.voxelScale = self.header[1]*self.header[2]*self.header[3] # mm^3
            self.voxelScale /= len(self.pointsPlotted)
            simplifiedMask = self.maskCoverImg[:,:,:,2]
            TIC = ut.generate_TIC(self.OGData4dImg, simplifiedMask, times, self.compressValue.value(),  self.voxelScale)

            # Bunch of checks
            if np.isnan(np.sum(TIC[:,1])):
                print('STOPPED:NaNs in the VOI')
                return;
            if np.isinf(np.sum(TIC[:,1])):
                print('STOPPED:InFs in the VOI')
                return;

            self.ticEditor = TICEditorGUI()
            self.ticEditor.show()
            self.ticEditor.accept.clicked.connect(self.acceptTIC)
            self.ticX = np.array([[TIC[i,0],i] for i in range(len(TIC[:,0]))])
            self.ticY = TIC[:,1]
            self.ticEditor.graph(self.ticX, self.ticY)
            self.ticEditor.initT0()
            self.ticEditor.t0Scroll.setValue(int(min(self.ticEditor.ticX[:, 0])))
            self.ticEditor.t0ScrollValueChanged()

    def showAuc(self):
        if self.aucParamapButton.isChecked():
            if self.id == 0:
                self.legend.setHidden(False)
            self.peParamapButton.setChecked(False)
            self.tpParamapButton.setChecked(False)
            self.mttParamapButton.setChecked(False)
            cmapStruct = plt.get_cmap('viridis')
            cmap = cmapStruct.colors
            self.maskCoverImg.fill(0)

            for i in range(len(self.pointsPlotted)):
                color = cmap[int((255/(self.maxAuc-self.minAuc))*(self.masterParamap[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]][0]-self.minAuc))]
                color = [color[i]*255 for i in range(len(color))]
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [color[2], color[1], color[0], int(self.curAlpha.value())]

            self.figLeg.clear()
            a = np.array([[0,1]])
            self.figLeg = plt.figure()
            img = plt.imshow(a, cmap='viridis')
            plt.gca().set_visible(False)
            if system == 'Windows':
                cax = plt.axes([0, 0.1, 0.35, 0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3, 0.45, "AUC", rotation=270, size=9) 
                plt.tick_params('y', labelsize=8, pad=0.5)
            else:
                cax = plt.axes([0, 0.1, 0.25, 0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3.1, 0.45, "AUC", rotation=270, size=4.5) 
                plt.tick_params('y', labelsize=4, pad=0.5)
            cax.set_yticks([0,0.25, 0.5, 0.75, 1])
            cax.set_yticklabels([int(self.minAuc*10)/10,int((((self.maxAuc-self.minAuc)/4)+self.minAuc)*10)/10,int((((self.maxAuc - self.minAuc)/2)+self.minAuc)*10)/10,int(((3*(self.maxAuc-self.minAuc)/4)+self.minAuc)*10)/10,int(self.maxAuc*10)/10])
            self.horizLayoutLeg.removeWidget(self.canvasLeg)
            self.canvasLeg = FigureCanvas(self.figLeg)
            self.horizLayoutLeg.addWidget(self.canvasLeg)
            self.canvasLeg.draw()

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

        elif self.aucParamapButton.isCheckable():
            self.maskCoverImg.fill(0)
            for i in range(len(self.pointsPlotted)):
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [0,0,255,int(self.curAlpha.value())]
            self.legend.setHidden(True)

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

            
    def showPe(self):
        if self.peParamapButton.isChecked():
            if self.id == 0:
                self.legend.setHidden(False)
            self.aucParamapButton.setChecked(False)
            self.tpParamapButton.setChecked(False)
            self.mttParamapButton.setChecked(False)
            cmapStruct = plt.get_cmap('magma')
            cmap = cmapStruct.colors
            self.maskCoverImg.fill(0)

            for i in range(len(self.pointsPlotted)):
                color = cmap[int((255/(self.maxPe-self.minPe))*(self.masterParamap[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]][1]-self.minPe))]
                color = [color[i]*255 for i in range(len(color))]
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [color[2], color[1], color[0], int(self.curAlpha.value())]

            self.figLeg.clear()
            a = np.array([[0,1]])
            self.figLeg = plt.figure()
            img = plt.imshow(a, cmap='magma')
            plt.gca().set_visible(False)
            if system == 'Windows':
                cax = plt.axes([0,0.1,0.35,0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3,0.45,"PE", rotation=270, size=9)
                plt.tick_params('y', labelsize=8, pad=0.5)
            else:
                cax = plt.axes([0, 0.1, 0.25, 0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3.1, 0.45, "PE", rotation=270, size=4.5) 
                plt.tick_params('y', labelsize=4, pad=0.5)
            cax.set_yticks([0,0.25, 0.5, 0.75, 1])
            cax.set_yticklabels([int(self.minPe*10)/10,int((((self.maxPe-self.minPe)/4)+self.minPe)*10)/10,int((((self.maxPe - self.minPe)/2)+self.minPe)*10)/10,int(((3*(self.maxPe-self.minPe)/4)+self.minPe)*10)/10,int(self.maxPe*10)/10])
            self.horizLayoutLeg.removeWidget(self.canvasLeg)
            self.canvasLeg = FigureCanvas(self.figLeg)
            self.horizLayoutLeg.addWidget(self.canvasLeg)
            self.canvasLeg.draw()

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

        elif self.peParamapButton.isCheckable():
            self.maskCoverImg.fill(0)
            self.figLeg = plt.figure()
            for i in range(len(self.pointsPlotted)):
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [0,0,255,int(self.curAlpha.value())]
            self.legend.setHidden(True)

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

    def showTp(self):
        if self.tpParamapButton.isChecked():
            if self.id == 0:
                self.legend.setHidden(False)
            self.aucParamapButton.setChecked(False)
            self.peParamapButton.setChecked(False)
            self.mttParamapButton.setChecked(False)
            cmapStruct = plt.get_cmap('plasma')
            cmap = cmapStruct.colors
            self.maskCoverImg.fill(0)

            for i in range(len(self.pointsPlotted)):
                color = cmap[int((255/(self.maxTp-self.minTp))*(self.masterParamap[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]][2]-self.minTp))]
                color = [color[i]*255 for i in range(len(color))]
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [color[2], color[1], color[0], int(self.curAlpha.value())]

            self.figLeg.cs()
            a = np.array([[0,1]])
            self.figLeg = plt.figure()
            img = plt.imshow(a, cmap='plasma')
            plt.gca().set_visible(False)
            if system == 'Windows':
                cax = plt.axes([0,0.1,0.35,0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3,0.45,"TP",rotation=270,size=9)
                plt.tick_params('y', labelsize=8, pad=0.5)
            else:
                cax = plt.axes([0, 0.1, 0.25, 0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3, 0.45, "TP", rotation=270, size=4.5) 
                plt.tick_params('y', labelsize=4, pad=0.5)
            cax.set_yticks([0,0.25, 0.5, 0.75, 1])
            cax.set_yticklabels([int(self.minTp*10)/10,int((((self.maxTp-self.minTp)/4)+self.minTp)*10)/10,int((((self.maxTp - self.minTp)/2)+self.minTp)*10)/10,int(((3*(self.maxTp-self.minTp)/4)+self.minTp)*10)/10,int(self.maxTp*10)/10])
            self.horizLayoutLeg.removeWidget(self.canvasLeg)
            self.canvasLeg = FigureCanvas(self.figLeg)
            self.horizLayoutLeg.addWidget(self.canvasLeg)
            self.canvasLeg.draw()

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

        elif self.tpParamapButton.isCheckable():
            self.maskCoverImg.fill(0)
            for i in range(len(self.pointsPlotted)):
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [0,0,255,int(self.curAlpha.value())]
            self.legend.setHidden(True)

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

    def showMtt(self):
        if self.mttParamapButton.isChecked():
            if self.id == 0:
                self.legend.setHidden(False)
            self.aucParamapButton.setChecked(False)
            self.peParamapButton.setChecked(False)
            self.tpParamapButton.setChecked(False)
            cmapStruct = plt.get_cmap('cividis')
            cmap = cmapStruct.colors
            self.maskCoverImg.fill(0)

            for i in range(len(self.pointsPlotted)):
                color = cmap[int((255/(self.maxMtt-self.minMtt))*(self.masterParamap[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]][3]-self.minMtt))]
                color = [color[i]*255 for i in range(len(color))]
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [color[2], color[1], color[0], int(self.curAlpha.value())]

            self.figLeg.clear()
            a = np.array([[0,1]])
            self.figLeg = plt.figure()
            img = plt.imshow(a, cmap='cividis')
            plt.gca().set_visible(False)
            if system == 'Windows':
                cax = plt.axes([0,0.1,0.35,0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3,0.43,"MTT",rotation=270,size=9)
                plt.tick_params('y',labelsize=8,pad=0.5)
            else:
                cax = plt.axes([0, 0.1, 0.25, 0.8])
                cbar = plt.colorbar(orientation='vertical', cax=cax)
                plt.text(3, 0.43, "MTT", rotation=270, size=4.5) 
                plt.tick_params('y', labelsize=4, pad=0.5)
            cax.set_yticks([0,0.25, 0.5, 0.75, 1])
            cax.set_yticklabels([int(self.minMtt*10)/10,int((((self.maxMtt-self.minMtt)/4)+self.minMtt)*10)/10,int((((self.maxMtt - self.minMtt)/2)+self.minMtt)*10)/10,int(((3*(self.maxMtt-self.minMtt)/4)+self.minMtt)*10)/10,int(self.maxMtt*10)/10])
            self.horizLayoutLeg.removeWidget(self.canvasLeg)
            self.canvasLeg = FigureCanvas(self.figLeg)
            self.horizLayoutLeg.addWidget(self.canvasLeg)
            self.canvasLeg.draw()

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

        elif self.mttParamapButton.isCheckable():
            self.maskCoverImg.fill(0)
            for i in range(len(self.pointsPlotted)):
                self.maskCoverImg[self.pointsPlotted[i][0], self.pointsPlotted[i][1], self.pointsPlotted[i][2]] = [0,0,255,int(self.curAlpha.value())]
            self.legend.setHidden(True)

            self.changeAxialSlices()
            self.changeSagSlices()
            self.changeCorSlices()

    

    def convertXmltoNifti(self):
        self.inputTextPath = ut.xml2nifti(self.inputXmlFileLocation, self.outputNiftiFileLocation)

        if self.inputTextPath:
            self.openInitialImageSlices()


def calculateSpline(xpts, ypts): # 2D spline interpolation
    cv = []
    for i in range(len(xpts)):
        cv.append([xpts[i], ypts[i]])
    cv = np.array(cv)
    if len(xpts) == 2:
        tck, u_ = interpolate.splprep(cv.T, s=0.0, k=1)
    elif len(xpts) == 3:
        tck, u_ = interpolate.splprep(cv.T, s=0.0, k=2)
    else:
        tck, u_ = interpolate.splprep(cv.T, s=0.0, k=3)
    x,y = np.array(interpolate.splev(np.linspace(0, 1, 1000), tck))
    return x, y

def ellipsoidFitLS(pos):

    # centre coordinates on origin
    pos = pos - np.mean(pos, axis=0)

    # build our regression matrix
    A = pos**2

    # vector of ones
    O = np.ones(len(A))

    # least squares solver
    B, resids, rank, s = np.linalg.lstsq(A, O, rcond=None)

    # solving for a, b, c
    a_ls = np.sqrt(1.0/B[0])
    b_ls = np.sqrt(1.0/B[1])
    c_ls = np.sqrt(1.0/B[2])

    return (a_ls, b_ls, c_ls)

def calculateSpline3D(points):
    # Calculate ellipsoid of best fit
    # points = np.array(points)
    # a,b,c = ellipsoidFitLS(points)
    # output = set()


    # u = np.linspace(0., np.pi*2., 1000)
    # v = np.linspace(0., np.pi, 1000)
    # u, v = np.meshgrid(u,v)

    # x = a*np.cos(u)*np.sin(v)
    # y = b*np.sin(u)*np.sin(v)
    # z = c*np.cos(v)

    # # turn this data into 1d arrays
    # x = x.flatten()
    # y = y.flatten()
    # z = z.flatten()
    # x += np.mean(points, axis=0)[0]
    # y += np.mean(points, axis=0)[1]
    # z += np.mean(points, axis=0)[2]

    # for i in range(len(x)):
    #     output.add((int(x[i]), int(y[i]), int(z[i])))
    # return output

    cloud = pv.PolyData(points, force_float=False)
    volume = cloud.delaunay_3d(alpha=100.)
    shell = volume.extract_geometry()
    final = shell.triangulate()
    final.smooth(n_iter=1000)
    faces = final.faces.reshape((-1, 4))
    faces = faces[:, 1:]
    arr = final.points[faces]

    arr = np.array(arr)

    output = set()
    for tri in arr:
        slope_2 = (tri[2]-tri[1])
        start_2 = tri[1]
        slope_3 = (tri[0]-tri[1])
        start_3 = tri[1]
        for i in range(100, -1, -1):
            bound_one = start_2 + ((i/100)*slope_2)
            bound_two = start_3 + ((i/100)*slope_3)
            cur_slope = bound_one-bound_two
            cur_start = bound_two
            for j in range(100, -1, -1):
                cur_pos = cur_start + ((j/100)*cur_slope)
                output.add((int(cur_pos[0]), int(cur_pos[1]), int(cur_pos[2])))
    
    return output