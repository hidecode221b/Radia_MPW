# Radia_MPW
Magnetic field simulation codes in Radia on Mathematica and Python for development of in-vacuum hybrid Undulator/Multipole Wiggler to generate higher flux and harder X-ray in the synchrotron storage ring.

## Features
- The number of poles can be varied either in even or odd. 
- Magnetic field can be plotted in the vector field based on PyVISTA.
- Magnet geometry can be plotted by PyVISTA from VTK file.
- 

## Environment
Radia python package is available from the link below.
https://github.com/ochubar/Radia

Radia python package for macosx is available from link below.
https://github.com/hidecode221b/Radia

Python code is based on the Radia ID developed by the ESRF, so the necessary files are required to run the code in addition to Radia python.
https://gitlab.esrf.fr/IDM/radia

Mathematica code is based on the Radia Mathematica example #03.
https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia

OpenGL is not available in macosx, but vtk file generated can be displayed by ParaView.
https://www.paraview.org/

Magnetic geometry and vector field plot by PyVISTA.
https://docs.pyvista.org/

Radia jupyter cloud service provided by RadiaSoft.
https://www.radiasoft.net/jupyter-server/

## Usage
Run the mpw.py after the main parameters are setup in mpw.py.

## Youtube
https://www.youtube.com/channel/UCeSigU1yYzv8_lmQxHp3fsg


![VectorField](https://github.com/hidecode221b/Radia_MPW/blob/main/images/Screen%20Shot%202022-02-01%20at%200.43.43.png "PyVista 3D magnetic field")
