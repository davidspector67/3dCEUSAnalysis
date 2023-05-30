# 3D Contrast Enhanced Ultrasound (CEUS) Analysis GUI

## Overview

This program is a CEUS analyis tool which allows user to input a string of 3D B-Mode images, draw a volume of interest, and view the resulting time intensity curve (TIC).

Next, the user can choose the point to start analysis on the TIC (at point t0), modify the TIC to remove noise from the data, and fit a lognormal curve to the resulting TIC.

### Initial TIC Image

![Initial TIC Image](images/initTICImage.png "Initial TIC Image")

### Modifying TIC

![Editing TIC Image](images/midTICImage.png "Modifying TIC Image")

### Final TIC Image

![Final TIC Image](images/finalTICImage.png "Final TIC Image")

Finally, using this lognormal curve fitting, the program computes the normalized area under the curve (AUC), normalized peak value (PV), normalized time to peak (TP), normalized mean transit time (MTT), and normalizing value (TMPPV).

### Main GUI

![Main GUI Image](images/imageGUI.png "Main GUI Image")

## Dependencies

* [Python Version 3.9.13](https://www.python.org/downloads/release/python-3913/)
* [Git](https://git-scm.com/downloads)

## Building

### Mac/Linux

```shell
git clone https://github.com/davidspector67/3dCEUSAnalysis.git
cd 3dCEUSAnalysis
chmod +x init.sh
chmod +x run.sh
./init.sh
```

### Windows

```shell
git clone https://github.com/davidspector67/3dCEUSAnalysis.git
cd 3dCEUSAnalysis
pip install -r .\pyPackages.txt
```

## Running

### Mac/Linux

```shell
./run.sh
```

### Windows

```shell
python main.py
```
