from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import platform
import numpy as np
from PyQt5.QtWidgets import QMainWindow, QLabel, QFrame, QSizePolicy, \
    QHBoxLayout, QPushButton, QApplication, QSlider
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, QRect, QMetaObject, QCoreApplication

from matplotlib.widgets import RectangleSelector

import lognormalFunctions as lf

system = platform.system()


selectedPoints = []
removedPointsX = []
removedPointsY = []
ticX = []
ticY = []
pickT0 = True
t0Point = [-1,-1,-1] #[xVal, yVal, indexVal]

global normalizer

class TICEditorGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setObjectName("TIC Editor")
        self.resize(1202, 731)
        self.removePointsLabel = QLabel(self)
        self.removePointsLabel.setGeometry(QRect(410, 0, 381, 81))
        font = QFont()
        font.setPointSize(25)
        self.removePointsLabel.setFont(font)
        self.removePointsLabel.setLayoutDirection(Qt.LeftToRight)
        self.removePointsLabel.setFrameShape(QFrame.NoFrame)
        self.removePointsLabel.setFrameShadow(QFrame.Plain)
        self.removePointsLabel.setAlignment(Qt.AlignCenter)
        self.removePointsLabel.setObjectName("removePointsLabel")
        self.removePointsLabel.setHidden(True)
        self.selectT0Label = QLabel(self)
        self.selectT0Label.setGeometry(QRect(370, 0, 461, 81))
        self.selectT0Label.setFont(font)
        self.selectT0Label.setLayoutDirection(Qt.LeftToRight)
        self.selectT0Label.setFrameShape(QFrame.NoFrame)
        self.selectT0Label.setFrameShadow(QFrame.Plain)
        self.selectT0Label.setAlignment(Qt.AlignCenter)
        self.selectT0Label.setObjectName("selectT0Label")
        self.selectT0Label.setHidden(True)
        self.chooseMethodLabel = QLabel(self)
        self.chooseMethodLabel.setGeometry(QRect(370, 0, 461, 81))
        self.chooseMethodLabel.setFont(font)
        self.chooseMethodLabel.setFrameShape(QFrame.NoFrame)
        self.chooseMethodLabel.setAlignment(Qt.AlignCenter)
        self.chooseMethodLabel.setText("Select t0 selection Method")
        self.chooseMethodLabel.setHidden(True)
        self.frame = QFrame(self)
        self.frame.setGeometry(QRect(50, 90, 1111, 551))
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.frame.setObjectName("frame")
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame.sizePolicy().hasHeightForWidth())
        self.frame.setSizePolicy(sizePolicy)
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.horizLayout = QHBoxLayout(self.frame)
        self.fig = plt.figure()
        self.canvas = FigureCanvas(self.fig)
        self.horizLayout.addWidget(self.canvas)
        self.canvas.draw()
        self.deselect = QPushButton(self)
        self.deselect.setGeometry(QRect(50, 660, 200, 50))
        self.deselect.setText("De-Select Last Point")
        self.deselect.setFocusPolicy(Qt.NoFocus)
        self.deselect.setHidden(True)
        self.remove = QPushButton(self)
        self.remove.setGeometry(QRect(350, 660, 200, 50))
        self.remove.setText("Remove Selected Points")
        self.remove.setFocusPolicy(Qt.NoFocus)
        self.remove.setHidden(True)
        self.restore = QPushButton(self)
        self.restore.setGeometry(QRect(650,660,200,50))
        self.restore.setText("Restore Last Points")
        self.restore.setFocusPolicy(Qt.NoFocus)
        self.restore.setHidden(True)
        self.accept = QPushButton(self)
        self.accept.setGeometry(QRect(950, 660, 200, 50))
        self.accept.setText("Accept Curve")
        self.accept.setFocusPolicy(Qt.NoFocus)
        self.accept.setHidden(True)

        self.acceptT0Button = QPushButton(self)
        self.acceptT0Button.setGeometry(QRect(660, 660, 250, 50))
        self.acceptT0Button.setFont(font)
        self.acceptT0Button.setText("Accept Start Time")
        self.acceptT0Button.setFocusPolicy(Qt.NoFocus)
        self.acceptT0Button.clicked.connect(self.acceptT0)
        self.acceptT0Button.setHidden(True)

        self.selectT0Button = QPushButton(self)
        self.selectT0Button.setGeometry(QRect(480, 660, 250, 50))
        self.selectT0Button.setFont(font)
        self.selectT0Button.setText("Select T0")
        self.selectT0Button.setFocusPolicy(Qt.NoFocus)

        self.t0index = -1
        self.frontPointsX = []
        self.frontPointsY = []

        self.deselect.clicked.connect(self.deselectLast)
        self.remove.clicked.connect(self.removeSelectedPoints)
        self.restore.clicked.connect(self.restoreLastPoints)
        self.selectT0Button.clicked.connect(self.initT0)

        self.retranslateUi(self)
        QMetaObject.connectSlotsByName(self)

    def t0ScrollValueChanged(self):
        self.prevLine.remove()
        self.prevLine = self.ax.axvline(x = self.t0Scroll.value(), color = 'green', label = 'axvline - full height')
        self.canvas.draw()

    def displayTicEditButtons(self):
        self.accept.setHidden(False)
        self.restore.setHidden(False)
        self.remove.setHidden(False)
        self.deselect.setHidden(False)
        self.removePointsLabel.setHidden(False)

    def acceptT0(self):
        global pickT0, selectedPoints, removedPointsX, removedPointsY
        if self.t0index == -1:
            for i in range(len(self.ticX[:,0])):
                if self.ticX[:,0][i] > self.t0Scroll.value():
                    break
            self.t0index = i
        self.displayTicEditButtons()
        self.acceptT0Button.setHidden(True)
        self.t0Scroll.setHidden(True)
        self.selectT0Label.setHidden(True)
        selectedPoints = list(range(self.t0index))
        if len(selectedPoints) > 0:
            self.removeSelectedPoints()
            self.frontPointsX = removedPointsX[-1]
            self.frontPointsY = removedPointsY[-1]
            removedPointsX.pop()
            removedPointsY.pop()
        self.ticX[:,0] -= (min(self.ticX[:,0]) - 1)
        self.fig.clear()
        self.graph(self.ticX, self.ticY)

    def rect_highlight(self, event1, event2):
        global selectedPoints
        self.mask |= self.inside(event1, event2)
        x = self.ticX[:,0][self.mask]
        y = self.ticY[self.mask]
        addedIndices = np.array(list(range(len(self.ticY))))[self.mask]
        for index in addedIndices:
            selectedPoints.append(index) 
        self.ax.scatter(x, y, color='orange')
        self.canvas.draw()

    def inside(self, event1, event2):
        # Returns a boolean mask of the points inside the rectangle defined by
        # event1 and event2
        x0, x1 = sorted([event1.xdata, event2.xdata])
        y0, y1 = sorted([event1.ydata, event2.ydata])
        mask = ((self.ticX[:,0] > x0) & (self.ticX[:,0] < x1) &
                (self.ticY > y0) & (self.ticY < y1))
        return mask

    def graph(self,x,y):
        global ticX, ticY
        y -= min(y)
        self.ticX = x
        self.ticY = y
        ticX = self.ticX
        ticY = self.ticY
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(x[:,0],y)
        self.ax.scatter(x[:,0],y,color='r')
        if system == 'Windows':
            self.ax.set_xlabel("Time (s)", fontsize=8, labelpad=0.5)
            self.ax.set_ylabel("Signal Amplitude", fontsize=8, labelpad=0.5)
            self.ax.set_title("Time Intensity Curve (TIC)", fontsize=10, pad=1.5)
            self.ax.tick_params('both', pad=0.3, labelsize=7.2)
            plt.xticks(fontsize=6)
            plt.yticks(fontsize=6)
        else:
            self.ax.set_xlabel("Time (s)", fontsize=11, labelpad=1)
            self.ax.set_ylabel("Signal Amplitude", fontsize=11, labelpad=1)
            self.ax.set_title("Time Intensity Curve (TIC)", fontsize=14, pad=1.5)
            self.ax.tick_params('both', pad=0.3, labelsize=7.2)
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
        plt.xticks(np.arange(0, int(max(self.ticX[:,0]))+10, 10))
        self.fig.subplots_adjust(left=0.1, right=0.97, top=0.9, bottom=0.1)

        if self.t0index > -1:
            self.mask = np.zeros(self.ticX[:,0].shape, dtype=bool)
            self.selector = RectangleSelector(self.ax, self.rect_highlight, useblit=True, props = dict(facecolor='cyan', alpha=0.2))

        self.canvas.draw()

    def initT0(self):
        self.acceptT0Button.setHidden(False)
        self.selectT0Button.setHidden(True)

        self.t0Scroll = QSlider(self)
        self.t0Scroll.setGeometry(QRect(400, 660, 200, 50))
        self.t0Scroll.setOrientation(Qt.Horizontal)
        self.t0Scroll.setHidden(False)

        self.t0Scroll.setMinimum(int(min(self.ticX[:,0])))
        self.t0Scroll.setMaximum(int(max(self.ticX[:,0])))
        self.t0Scroll.valueChanged.connect(self.t0ScrollValueChanged)
        self.t0Scroll.setValue(int(min(self.ticX[:, 0])))
        self.prevLine = self.ax.axvline(x = self.t0Scroll.value(), color = 'green', label = 'axvline - full height')
        self.canvas.draw()

        # def manualChecked(self):
    #     self.hidePickingButtons()
    #     self.acceptT0Button.setHidden(False)
        # self.t0Scroll.setMinimum(int(min(self.ticX[:,0])))
        # self.t0Scroll.setMaximum(int(max(self.ticX[:,0])))
        # self.t0Scroll.valueChanged.connect(self.t0ScrollValueChanged)
        # self.t0Scroll.setValue(int(min(self.ticX[:,0])))
        # self.prevLine = self.ax.axvline(x = self.t0Scroll.value(), color = 'green', label = 'axvline - full height')
        # self.canvas.draw()
        # self.t0Scroll.setHidden(False)
        # self.selectT0Label.setHidden(False)

    def deselectLast(self):
        global selectedPoints
        if len(selectedPoints) > 0:
            lastPt = selectedPoints[-1]
            selectedPoints.pop()
            self.ax.scatter(self.ticX[lastPt][0],self.ticY[lastPt],color='red')
            self.canvas.draw()

    def removeSelectedPoints(self):
        global selectedPoints, removedPointsX, removedPointsY, ticX, ticY
        if len(selectedPoints) > 0:
            selectedPoints.sort()
            j = 0
            curRemovedX = []
            curRemovedY = []
            for i in range(len(self.ticX)):
                if i == selectedPoints[j]:
                    curRemovedX.append(ticX[i])
                    curRemovedY.append(ticY[i])
                    j += 1
                    if j == len(selectedPoints):
                        break
            self.ticX = np.delete(self.ticX, selectedPoints, axis=0)
            self.ticY = np.delete(self.ticY, selectedPoints)
            self.fig.clear()
            self.graph(self.ticX,self.ticY)
            ticX = self.ticX
            ticY = self.ticY
            removedPointsX.append(curRemovedX)
            removedPointsY.append(curRemovedY)
            selectedPoints = []

    def restoreLastPoints(self):
        global selectedPoints, removedPointsX, removedPointsY, ticX, ticY
        if len(removedPointsX) > 0:
            for i in range(len(selectedPoints)):
                self.deselectLast()
            selectedPoints = []
            j = 0
            i = 0
            max = self.ticX.shape[0] + len(removedPointsX[-1])
            while i < self.ticX.shape[0]-1:
                if self.ticX[i][0] < removedPointsX[-1][j][0] and removedPointsX[-1][j][0] < self.ticX[i+1][0]:
                    self.ticX = np.insert(self.ticX, i+1, removedPointsX[-1][j], axis=0)
                    self.ticY = np.insert(self.ticY, i+1, removedPointsY[-1][j])
                    j += 1
                    if j == len(removedPointsX[-1]):
                        break
                i += 1
            if i < max and j < len(removedPointsX[-1]):
                while j < len(removedPointsX[-1]):
                    self.ticX = np.insert(self.ticX, i+1, removedPointsX[-1][j], axis=0)
                    self.ticY = np.append(self.ticY, removedPointsY[-1][j])
                    j += 1
                    i += 1
            removedPointsX.pop()
            removedPointsY.pop()
            self.fig.clear()
            ticX = self.ticX
            ticY = self.ticY
            self.graph(ticX, ticY)


    def retranslateUi(self, roughGUI):
        _translate = QCoreApplication.translate
        roughGUI.setWindowTitle(_translate("TIC Editor", "TIC Editor"))
        self.removePointsLabel.setText(_translate("TIC Editor", "Select Points to Remove from TIC:"))
        self.selectT0Label.setText(_translate("TIC Editor", "Select Time to start Analysis"))



if __name__ == "__main__":
    import sys
    from numpy import genfromtxt
    # my_data = genfromtxt('/Users/davidspector/Home/Stanford/USImgAnalysisGui_v2/Data/newest_test_tic.csv', delimiter=',')[1:]
    my_data = genfromtxt('/Users/davidspector/Home/Stanford/USImgAnalysisGui_v2/Data/C3P13_original_tic.csv', delimiter=',')[1:]
    # my_data = genfromtxt('/Users/davidspector/Home/Stanford/USImgAnalysisGui_v2/Data/C3P13_original_tic.csv', delimiter=',')

    test_ticX = np.array([[my_data[i,0],i] for i in range(len(my_data[:,0]))])
    # test_ticY = my_data[:,1]
    test_ticY = my_data[:,1] - min(my_data[:,1])

    normalizer = max(test_ticY)

    print(np.average(my_data[:,1]))


    app = QApplication(sys.argv)
    ui = TICEditorGUI()
    ui.show()
    ui.graph(test_ticX, test_ticY)
    sys.exit(app.exec_())
