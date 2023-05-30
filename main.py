import sys
import PyQt5
from PyQt5.QtWidgets import QApplication
from analysis3dController import Contrast3dAnalysisController

if __name__ == "__main__":
  ceusApp = QApplication(sys.argv)
  ceusUI = Contrast3dAnalysisController()
  ceusUI.show()
  sys.exit(ceusApp.exec_()) 
