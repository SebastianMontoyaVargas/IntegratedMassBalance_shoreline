############################################   General description  ########################################################################

The present folder contains the codes and data set used by S. Montoya-Vargas, A.F. Osorio & O. Duran in the paper entitled 
"INTEGRATED BALANCE EQUATION FOR SHORELINE EVOLUTION ANALYSIS COUPLING CROSS- AND LONG-SHORE SEDIMENT TRANSPORT PROCESSES" 
(2022). The data set consists on .nc files containing the bathymetric and hydrodynamic data for Duck Beach (NC) collected by 
the U.S. Army Military Field Research Facilities at Duck and downloaded from http://www.frf.usace.army.mil/frf_data.shtml.
The codes built in Matlab allow the user to read the data and carry the analysis presented by Montoya-Vargas et al. (2022) as
well as to reproduce the figures presented in the paper.

############################################       Instructions     ########################################################################

Execute the .m files in the DATA_PROCESSING folder in the following order (it might take several minutes):

1. shoreline_from_DEM.m
2. closure_from_DEM.m
3. waterlevel_processing.m
4. wave_read.m
5. waves_mean.m
6. Volume_from_DEM.m
7. sediment_transport.m
8. time_derivatives.m
9. compute_correlations.m

Then run the file "plots.m" for reproducing Figures 4, and 6 to 11 from paper.

############################################    Files description   ########################################################################

The following folders and files are contained inside the "POSTPROCESSING_CODES" folder:

- plots.m: Reproduce the figures presented in Montoya-Vargas et al. (2022)

- Read_me.txt: This file

- DATA_PROCESSING: Raw and processed data along with pre- and post-processing codes. It contains the following subfolders and files:
	
	- closure_from_DEM.m : Matlab code for computation of closure contour and elevation as described by Montoya-Vargas et al. (2022).
	- compute_correlations.m : Matlab code for computation of individual terms correlations for different time windows.
	- Correlations.mat : Matlab file containing the data calculated with "compute_correlations.m".
	- Derivatives.mat : Matlab file containing the data calculated with "time_derivatives.m". 
	- hc_DEM.mat :  Matlab file containing the data calculated with "closure_from_DEM.m".
	- hs.mat :  Matlab file containing the data calculated with "waterlevel_processing.m".
	- LST.mat :  Matlab file containing the data calculated with "sediment_transport.m".
	- mean_waves.mat :  Matlab file containing the data calculated with "waves_mean.m".
	- sediment_transport.m : Matlab code for computation of bulk long-shore sediment transport (Qy).
	- shoreline_from_DEM.m : Matlab code for computation of shoreline position as described by Montoya-Vargas et al. (2022).
	- time_derivatives.m : Matlab code for computation of time derivatives dL/dt, dDc/dt, db/dt, dhs/dt and dxs/dt.
	- vol_DEM.mat : Matlab file containing the data computed with "Volume_from_DEM.m".
	- waterlevel_processing.m : Matlab code for computation of sea level as described by Montoya-Vargas et al. (2022).
	- wave_read.m : Matlab code for reading raw wave data.
	- waves_mean.m :  Matlab code for computation of mean wave parameters.
	- WAVES_RAW.mat :  Matlab file containing the data readed on "wave_read.m".
	- Xc_DEM.mat : Matlab file containing the closure contours readed with "closure_from_DEM.m".
	- Xs_DEM.mat : Matlab file containing the shoreline positions readed with "shoreline_from_DEM.m".
	- 26m_waverider: Folder containing the .nc files with wave data including significant wave height Hs, wave peak period Tp, and peak direction Dir.
	- DUCK_DEM : Folder containing the .nc files with bathymetric data.
	- DUCK_TIDES : Folder containing the .nc files with sea level.

Further description of the analysis methods and data can be found elsewhere (Montoya-Vargas et a., 2022).


