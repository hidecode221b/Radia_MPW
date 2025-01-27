In-vacuum multipole wiggler (IVW) simulation is based on the design of MPW at BL1W in SPL (SLRI, Thailand). The MPW was originally designed and constructed in ASTeC, UK. 

Optimizing parameters are required to initialize the parameters prior to rerunning the different optimization.

Environment: Wolfram Mathematica 13.0.0.0 on Windows 10, Radia version 4.31

Rev. 20220926
Eample#3mpw_noc.nb represents the MPW mode without pole chamfers, and exports data in text for field, field integral, and trajectories.

Rev. 20221012
Eample#3mpw_noc.nb exports data files to plot in the software. Igor procedure is uploaded for a 3D magnetic vector in a gizmo plot.

Rev. 20221018
Eample#3mpw_noc.nb exports VTK files on geometry and magnetization field vectors.

Rev. 20221118
Eample#3mpw_noc_pgn.nb can make the side chamfers on magnets based on the pgn instead of rtg.

Rev. 20230815
Eample#3mpw_noc_pgn_nos.nb can build the magnet arrays without side magnets, which enables us to set up the same width in pole and magnet. The uniform misalignment of magnetization can also be set up, but the symmetrical translation canceled the misalignment.

Rev. 20230824
Eample#3mpw_noc_rtg_nos.nb is based on the rtg instead of png without side magnets. The randomized misalignment of magnetization is implemented without symmetrical translation. std = 0 makes ideal magnetization, and std <> 0 results in the misorientation of magnetization in the standard deviation of sdr (%/100) and sdt (degrees) in magnetization amplitude (y direction) and angle (theta), respectively. Phi angle is uniformly randomized around the magnetization direction. 

Rev. 20230901
Eample#3mpw_noc_rtg_nos_for.nb is updated for for-loop in +/-y and +/-z axis with variables side=+/-1 and array=+/-1. The separation of each component in the y-axis can also be randomized to evaluate the magnetization errors.

Rev. 20240628
Eample#3mpw_noc_rtg_nos_for.nb is updated for kick map. Eample#3ivu_noc_rtg_nos_for.nb is also added for hybried pole optimization at the periodic length of 20 mm. Optimization module is revised to change the pole and magnetic lengths at the constant periodic length conditon.

Rev. 20240709
Eample#3mpw_noc_rtg_nos_for.nb is updated for the kick map in T2m2 for horizontal and vertical kicks separately exporting the data files with coordinates.

Rev. 20240815
Eample#3mpw_noc_rtg_nos_for.nb is updated for the gap dependence based on sp[[3]] in the module.

Rev. 20250103
Eample#3mpw_noc_rtg_nos_for.nb is updated for the magnetic force simulation in cut models.

Rev. 20250127
Eample#3mpw_noc_rtg_nos_for.nb is updated for the kickmap in a period and total length.

