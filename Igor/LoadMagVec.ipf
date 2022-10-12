#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

Menu "Macros"
	"Load Mag Vec"
End

Function LoadMagVec() : Panel
	String panelmagvec = "Load_Mag_Vec"
	DoWindow panelmagvec
	if (V_flag == 1)
		DoWindow/F panelmagvec
		Abort
	endif
	
	if (DataFolderExists("MG") ==0)
		String/G File_Name = "", Full = "", Path = ""
		DFREF saveDFR = GetDataFolderDFR()
		NewDataFolder/S/O root:MG
		String/G cpopupList = "HSV;BWR;Rainbow"
		SetDataFolder saveDFR
	else
		String/G File_Name, Full, Path
	endif
	
	NewPanel /W=(630,30,960,230) as panelmagvec
	
	SetVariable setpathp pos={10,10},size={210,30},title="Path:",value=Path
	SetVariable setnamep pos={10,32},size={210,30},title="File Name:",value=File_Name
	
	Button brwbudcp pos={230,10},size={80,20},title="Browse",proc=Browse_mv
	Button loaddcp pos={230,30},size={80,20},title="Load",proc=Load_mv
	
	Button plotgismov pos={230,60},size={80,20},title="Plot gismo",proc=Plot_gismo
	Button plotgismoo pos={230,80},size={80,20},title="Arrow gismo",proc=Arrow_gismo
	
	Slider/Z setlinesizev, pos={100,60},size={120,0}, vert=0, ticks=0, side=0, limits={0.02, 0.5, 0}, proc=linesizev, thumbColor=(0,1000,0)
	Slider/Z setversizev, pos={100,80},size={120,0}, vert=0, ticks=0, side=0, limits={0.1, 0.5, 0}, proc=versizev, thumbColor=(0,1000,0)
	Slider/Z setbasesizev, pos={100,100},size={120,0}, vert=0, ticks=0, side=0, limits={0.05, 0.2, 0}, proc=basesizev, thumbColor=(0,1000,0)
	Slider/Z setheadsizev, pos={100,120},size={120,0}, vert=0, ticks=0, side=0, limits={0.05, 0.5, 0}, proc=headsizev, thumbColor=(0,1000,0)
	Slider/Z setcylsizev, pos={100,140},size={120,0}, vert=0, ticks=0, side=0, limits={0.01, 0.1, 0}, proc=cylsizev, thumbColor=(0,1000,0)

	TitleBox tit_sl0,title="Line size",pos={10,60},size={20, 20},frame=0
	TitleBox tit_sl1,title="Vertex z",pos={10,80},size={20, 20},frame=0
	TitleBox tit_sl2,title="Arrow base",pos={10,100},size={20, 20},frame=0
	TitleBox tit_sl3,title="Arrow head",pos={10,120},size={20, 20},frame=0
	TitleBox tit_sl4,title="Cylinder",pos={10,140},size={20, 20},frame=0
	
	PopupMenu setcolors pos={10,165}, mode=1, value=cPopupWaveList(), proc=SetColors, title="Colors"
End

Function linesizev(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, 
				//   2: mouse up, 3: mouse moved
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	ModifyGizmo/N=$wn/Z modifyObject=scatter0, objectType=scatter, property={size,value}
			
	return 0	// other return values reserved
End

Function versizev(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, 
				//   2: mouse up, 3: mouse moved
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	ModifyGizmo/N=$wn/Z modifyObject=line0, objectType=line, property={vertex,0,0,0,0,0,value}
			
	return 0	// other return values reserved
End

Function basesizev(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, 
				//   2: mouse up, 3: mouse moved
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	Modifygizmo/N=$wn/Z modifyobject=line0,objectType=line,property={endArrowBase,value}
			
	return 0	// other return values reserved
End

Function headsizev(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, 
				//   2: mouse up, 3: mouse moved
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	Modifygizmo/N=$wn/Z modifyobject=line0,objectType=line,property={endArrowHeight,value}
			
	return 0	// other return values reserved
End

Function cylsizev(name, value, event)
	String name	// name of this slider control
	Variable value	// value of slider
	Variable event	// bit field: bit 0: value set; 1: mouse down, 
				//   2: mouse up, 3: mouse moved
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	Modifygizmo/N=$wn/Z modifyobject=line0,objectType=line,property={cylinderStartRadius,value}
	Modifygizmo/N=$wn/Z modifyobject=line0,objectType=line,property={cylinderEndRadius,value}
			
	return 0	// other return values reserved
End

Function Browse_mv(ctrlName) : ButtonControl
	String ctrlName
	String/G File_Name, Full, Path
	Variable Read_File

	String fileFilters = "Data Files (*.dat):.dat;"
	Open /R/F=fileFilters/Z=2 Read_File
	
	Variable err = V_flag
	String fullPath = S_fileName
	if (err == -1)
		return -1
	endif
	
	FStatus Read_File
	File_Name = S_filename
	Full = S_path + S_filename
	Path =S_path
	Close Read_File
	return 0
End

Function Load_mv(ctrlName) : ButtonControl
	String ctrlName
	String/G File_Name, Full
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	DFREF saveDFR = GetDataFolderDFR()
	NewDataFolder/O/S :$wn
	
	LoadWave /O/K=0/Q/N=$wn +"_a"/J/L={0,0,0,0,0} Full	//Skip loading 
	Duplicate/O $wn +"_a0", $wn +"_x0"	
	Duplicate/O $wn +"_a1", $wn +"_y0"		
	Duplicate/O $wn +"_a2", $wn +"_z0"	
	Duplicate/O $wn +"_a3", $wn +"_x1"	
	Duplicate/O $wn +"_a4", $wn +"_y1"		
	Duplicate/O $wn +"_a5", $wn +"_z1"
		
	Variable i=0
	do
		if ( !WaveExists($wn + "_a" + num2str(i)) )	// Delete all imported data.
			break
		endif
		KillWaves/Z $wn + "_a" + num2str(i)
		i += 1
	while(1)
	SetDataFolder saveDFR
End

Function/S Remove_fnmv(str, str2)
	String str, str2
	
	Variable pos= strsearch(str, str2, 0)
	if( pos < 0 )
		return str
	endif
	return str[0,pos-1]
End

Function Plot_gismo(ctrlName) : ButtonControl
	String ctrlName
	String/G File_Name, Full, Path
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	DFREF saveDFR = GetDataFolderDFR()
	SetDataFolder $wn
	
	Wave x0 = $wn + "_x0", y0 = $wn + "_y0", z0 = $wn + "_z0"
	Wave x1 = $wn + "_x1", y1 = $wn + "_y1", z1 = $wn + "_z1"
	
	Duplicate/O $wn +"_x0", $wn +"_ag", $wn +"_nm", $wn +"_cx", $wn +"_cy", $wn +"_cz"
	Duplicate/O $wn +"_x0", $wn +"_vx", $wn +"_vy", $wn +"_vz"
	Duplicate/O $wn +"_x0", $wn +"_dx", $wn +"_dy", $wn +"_dz"
	
	Wave cx = $wn + "_cx", cy = $wn + "_cy", cz = $wn + "_cz"
	Wave dx = $wn + "_dx", dy = $wn + "_dy", dz = $wn + "_dz"
	Wave vx = $wn + "_vx", vy = $wn + "_vy", vz = $wn + "_vz"
	Wave ag = $wn + "_ag", nm = $wn + "_nm"
	
	// strength
	nm=sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)
	// angle of rotation specified in the direction below
	ag=180*acos((z1-z0)/nm)/pi

	// direction of rotation, initial arrow in z 
	vx=((y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2))
	vy=(-(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2))
	vz=0

	// difference of 2 points
	dx=(x1-x0)/nm
	dy=(y1-y0)/nm
	dz=(z1-z0)/nm
	
	// center of 2 points
	cx=(x0+x1)/2
	cy=(y0+y1)/2
	cz=(z0+z1)/2
	
	Make/O/D/N=(numpnts($wn+"_x0"), 3) $wn + "_gs", $wn + "_sz"
	Make/O/D/N=(numpnts($wn+"_x0"), 4) $wn + "_cl", $wn + "_rt"
	Wave gs = $wn + "_gs", sz = $wn + "_sz", cl = $wn + "_cl", rt = $wn + "_rt"
	Variable wmax = WaveMax(nm), wmin = WaveMin(nm)
	
	// normalized in 0-1 range
	nm = (nm - wmin)/(wmax - wmin)
	Variable dmax = WaveMax(dz), dmin = WaveMin(dz)
	dx = 0.1*(dx - dmin)/(dmax - dmin)
	dy = 0.1*(dy - dmin)/(dmax - dmin)
	dz = 0.1*(dz - dmin)/(dmax - dmin)
	
	Variable i = 0
	do
		gs[i][0] = cx[i]	// position
		gs[i][1] = cy[i]
		gs[i][2] = cz[i]
		sz[i][0] = dx[i]	// size
		sz[i][1] = dy[i]
		sz[i][2] = dz[i]
		rt[i][0] = ag[i]	// angle
		rt[i][1] = vx[i]	// direction
		rt[i][2] = vy[i]
		rt[i][3] = vz[i]
		i += 1
	while (i < numpnts(x0))
	
	i = 0
	Variable v=1, s=1, j, f	// default color HSV2RGB, h: hue, s: saturation, v: brightness
	do	
		nm[i] *= 6
		j = floor(nm[i])
		f = nm[i] - j
		if (j == 0)
			cl[i][0] = v
			cl[i][1] = v*(1-s*(1-f))
			cl[i][2] = v*(1-s)
		elseif (j == 1)
			cl[i][0] = v*(1-s*f)
			cl[i][1] = v
			cl[i][2] = v*(1-s)
		elseif (j == 2)
			cl[i][0] = v*(1-s)
			cl[i][1] = v
			cl[i][2] = v*(1-s*(1-f))
		elseif (j == 3)
			cl[i][0] = v*(1-s)
			cl[i][1] = v*(1-s*f)
			cl[i][2] = v
		elseif (j == 4)
			cl[i][0] = v*(1-s*(1-f))
			cl[i][1] = v*(1-s)
			cl[i][2] = v
		elseif (j == 5)
			cl[i][0] = v
			cl[i][1] = v*(1-s)
			cl[i][2] = v*(1-s*f)
		endif
		cl[i][3] = nm[i]	// transparency alpha [0-1]
		i += 1
	while (i < numpnts(nm))
	
	DoWindow $wn
	if (V_flag == 1)
		DoWindow/K $wn
	endif
	
	// gizmo scatter plot with cone shape
	NewGizmo/I/N=$wn
	AppendTogizmo/N=$wn/Z defaultScatter=gs, name=scatter0
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={shape,9}
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={size,0.1}
	//ModifyGizmo modifyObject=scatter0, objectType=scatter, property={sizeType,1}
	//ModifyGizmo modifyObject=scatter0, objectType=scatter, property={sizeWave,sz}
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={rotationType,1}
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={rotationWave,rt}
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={scatterColorType,1}
	ModifyGizmo modifyObject=scatter0, objectType=scatter, property={colorWave,cl}
	
	SetDataFolder saveDFR
	DoWindow $"Load_Mag_Vec"	// not working?
	Print(V_flag)
	DoWindow/F panelmagvec
end

Function/S cPopupWaveList()
	SVAR cpopupList=root:MG:cpopupList
	String list = cpopupList
	return list
End

Function SetColors(ctrlName,popNum,popStr) : PopupMenuControl
	String ctrlName,popStr
	Variable popNum

	String/G File_Name, Full, Path
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	DFREF saveDFR = GetDataFolderDFR()
	SetDataFolder $wn
	Wave nm = $wn + "_nm", cl = $wn + "_cl"
	Variable wmax = WaveMax(nm), wmin = WaveMin(nm)
	
	nm = (nm - wmin)/(wmax - wmin)
	Variable i = 0, v=1, s=1, j, f
	
	if (popNum == 1)	// HSV default used in mathematica radia
		do
			nm[i] *= 6
			j = floor(nm[i])
			f = nm[i] - j
			if (j == 0)
				cl[i][0] = v
				cl[i][1] = v*(1-s*(1-f))
				cl[i][2] = v*(1-s)
			elseif (j == 1)
				cl[i][0] = v*(1-s*f)
				cl[i][1] = v
				cl[i][2] = v*(1-s)
			elseif (j == 2)
				cl[i][0] = v*(1-s)
				cl[i][1] = v
				cl[i][2] = v*(1-s*(1-f))
			elseif (j == 3)
				cl[i][0] = v*(1-s)
				cl[i][1] = v*(1-s*f)
				cl[i][2] = v
			elseif (j == 4)
				cl[i][0] = v*(1-s*(1-f))
				cl[i][1] = v*(1-s)
				cl[i][2] = v
			elseif (j == 5)
				cl[i][0] = v
				cl[i][1] = v*(1-s)
				cl[i][2] = v*(1-s*f)
			endif
			cl[i][3] = 0
			i += 1
		while (i < numpnts(nm))
	elseif (popNum == 2)	// BWR: blue-white-red
		do
			if (nm[i] < 0.5)
				cl[i][0] = 2*nm[i]
				cl[i][1] = 2*nm[i]
				cl[i][2] = 1
			else
				cl[i][0] = 1
				cl[i][1] = 1-2*(nm[i]-0.5)
				cl[i][2] = 1-2*(nm[i]-0.5)
			endif
			cl[i][3] = 0
			i += 1
		while (i < numpnts(nm))
	elseif (popNum == 3)	// rainbow
		do
			cl[i][0] = exp(-0.5*((nm[i]-0.8)^2/(0.3)^2))
			cl[i][1] = exp(-0.5*((nm[i]-0.35)^2/(0.35)^2))
			cl[i][2] = exp(-0.5*((nm[i]-0.1)^2/(0.2)^2))
			cl[i][3] = 0
			i += 1
		while (i < numpnts(nm))
	endif
	SetDataFolder saveDFR
	
End

Function Arrow_gismo(ctrlName) : ButtonControl
	String ctrlName
	String/G File_Name, Full, Path
	
	String wn = Remove_fnmv(File_Name, ".")
	wn = Remove_fnmv(wn, "_vec")
	
	// add line object with arrow setting
	AppendToGizmo/N=$wn/Z line={0,0,0,0,0,0.2},name=line0
	Modifygizmo modifyobject=line0,objectType=line,property={arrowMode,18}
	Modifygizmo modifyobject=line0,objectType=line,property={endArrowHeight,0.08}
	Modifygizmo modifyobject=line0,objectType=line,property={endArrowBase,0.07}
	Modifygizmo modifyobject=line0,objectType=line,property={cylinderStartRadius,0.02}
	Modifygizmo modifyobject=line0,objectType=line,property={cylinderEndRadius,0.02}
	
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,7}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ objectName,line0}
	
	DoWindow/F panelmagvec
End