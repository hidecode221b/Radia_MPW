# Radia_MPW
Magnetic field simulation codes in Radia on Mathematica and Python for development of in-vacuum hybrid Undulator/Multipole Wiggler to generate higher flux and harder X-ray in the synchrotron storage ring.

## Highlight features
- The number of poles can be varied either in even or odd (Mathematica). 
- Magnetic field can be plotted in the vector field based on PyVISTA (Python).
- Magnet geometry can be displayed by PyVISTA from VTK file (Python).
- ParaView overlays magnetic field and magnet geometry based on VTK files.
- Calculate the magnetic force on half array from the whole structure (Mathematica and Python).
- Kickmap based on the potential differentiation in horizontal and vertical kicks (Mathematica and Python).

## Abstract
The insertion devices generate a high energy and flux of synchrotron radiation without major upgrade of 3rd generation storage ring. The undulator utilizes the interference of radiations from the periodic magnet array to produce the X-rays with high flux rather than high energy. The wiggler creates the higher magnetic field than the bending magnet does to produce the X-rays with high energy rather than high flux. The permanent magnets are limited in its magnetization, so the superconducting magnets can be an option for the creation of highest magnetic field. However, the cost of cryogenic facility and its maintenance is typically expensive.

The soft magnetic materials with high permeability combined with permanent magnets enables us to achieve higher magnetic field than the pure permanent magnet does to concentrate and control the magnetic flux in a specific direction. This is a so-called hybrid type (HYB or Halbach II) named after the pure permanent magnet undulator (PPM or Halbach I). However, because of non-linearity of pole materials and beam dynamics, the magnetic structure has to be carefully designed to maximize the magnetic field under the operating conditions of storage ring.

For a 3-GeV ring in the SPS II project, the in-vacuum undulator is used to produce the soft and hard X-rays for beamlines. To obtain the know-how of production of in-vacuum undulator, we start to build the prototype of in-vacuum undulator or wiggler. To test the prototype, we choose the in-vacuum wiggler to replace the superconducting wavelength shifter (SWLS) based on the cryogenic facility which has a problem currently. However, the technical feasibility has to be taken into account for the SPS storage ring accommodated with SWLS. To simplify the designing process, the out-vacuum multipole wiggler in SPS BL1 is used as a basic design for the in-vacuum multipole wiggler to be developed. The MPW at BL1 was designed, constructed, and donated by the ASTeC, UK.

This repository stores the Radia codes on Mathematica and Python for the magnetic field simulation of the in-vacuum wiggler prototype.

## Environment
Radia python package is available from the link below.

https://github.com/ochubar/Radia

Radia python for macOS is available from link below.

https://github.com/hidecode221b/Radia

Python class is based on the Radia ID developed by the ESRF, so the necessary files are required to run the code in addition to Radia python.

https://gitlab.esrf.fr/IDM/radia

Mathematica code is based on the Radia Mathematica example #03.

https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia

OpenGL is not available in macosx, but vtk file generated can be displayed by ParaView.

https://www.paraview.org/

Magnetic geometry and vector field plot by PyVISTA.

https://docs.pyvista.org/

Radia jupyter cloud service provided by RadiaSoft.

https://www.radiasoft.net/jupyter-server/

https://github.com/radiasoft/Radia-Examples

## Model of MPW in calculation
MPW simulation is based on the design of MPW at BL1 in SPL (SLRI, Thailand). However, the MPW was designed, contructed, and donated by the ASTeC, UK.

https://www.astec.stfc.ac.uk/

## Usage
### Python from RadiaID gitlab
All the python files should be located under the radia-master/env/radia_python/.

Radia.so should be created for your python environment (See [link](https://github.com/hidecode221b/Radia) for macos).

Run the mpw.py after the main parameters are setup in mpw.py and radia_id_mpw_hyb.py 

Examples object and VTK files are available in examples directory.

### Mathematica from ESRF Radia website
Example#3mpw.nb

### Jupyter notebook with RadiaID
Wiggler_HYB.ipynb

radia_id.py

### Jupyter notebook based on RadiaSoft Jupyter server
Radia_Ex03r.ipynb

## Issues and testing in progress
- Python Radia is now working for the magnetic force calculation.

## Demo in Youtube
https://www.youtube.com/channel/UCeSigU1yYzv8_lmQxHp3fsg


![VectorField](https://github.com/hidecode221b/Radia_MPW/blob/main/images/Screen%20Shot%202022-02-03%20at%2020.05.15.png "ParaView 3D magnetic field visualization with magnet geometry")

## References
- The Science and Technology of Undulators and Wigglers, James A. Clarke (Oxford University Press, 2004) https://doi.org/10.1093/acprof:oso/9780198508557.001.0001
- Bahrdt J. (2016) Shaping Photon Beams with Undulators and Wigglers. In: Jaeschke E., Khan S., Schneider J., Hastings J. (eds) Synchrotron Light Sources and Free-Electron Lasers. Springer, Cham. https://doi.org/10.1007/978-3-319-14394-1_16
- Undulators, Wigglers and Their Applications, Edited By Hideo Onuki, Pascal Elleaume (CRC Press 2003) https://doi.org/10.4324/9780203218235 
