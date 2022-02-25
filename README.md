# Radia_MPW
Magnetic field simulation codes in Radia on Mathematica and Python for development of in-vacuum hybrid Undulator/Multipole Wiggler to generate higher flux and harder X-ray in the synchrotron storage ring.

## Features
- The number of poles can be varied either in even or odd. 
- Magnetic field can be plotted in the vector field based on PyVISTA.
- Magnet geometry can be displayed by PyVISTA from VTK file.
- ParaView overlays magnetic field and magnet geometry based on VTK files.
- Calculate the magnetic force on half pole array from the whole structure.

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

## Model of MPW in calculation
MPW simulation is based on the design of MPW at BL1 in SPL (SLRI, Thailand). However, the MPW was designed, contructed, and donated by the ASTeC, UK.
https://www.astec.stfc.ac.uk/

## Usage
### Python from RadiaID gitlab
Run the mpw.py after the main parameters are setup in mpw.py and radia_id_mpw_hyb.py 

Examples object and VTK files are available in examples directory.

### Mathematica from ESRF Radia website
Example#3mpw.nb

### Jupyter notebook with RadiaID
Wiggler_HYB.ipynb

### Jupyter notebook based on RadiaSoft Jupyter server
Radia_Ex03r.ipynb

## Issues and testing in progress
- 

## Youtube
https://www.youtube.com/channel/UCeSigU1yYzv8_lmQxHp3fsg


![VectorField](https://github.com/hidecode221b/Radia_MPW/blob/main/images/Screen%20Shot%202022-02-03%20at%2020.05.15.png "ParaView 3D magnetic field visualization with magnet geometry")
