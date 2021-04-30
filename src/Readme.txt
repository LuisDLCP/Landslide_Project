****************************
        README
****************************

Overview
---------
	These Python scripts generate SAR images, interferograms, 
	and statistical graphs, by using raw data acquired by the
	GB-SAR radar "IGP-ROJ".  


Input files 
------------
	*.hdf5 files 

Output files
-------------
	SAR images: *.png files
	Interferograms: *.png files
	Statistical graphs: *.png files 

Scripts
--------
	Main file: 
	----------
	"Interferometria_main.py"

	Dependencies:
	-------------
	"BP_real_main.py": The FDBP Algorithm is implemented here.
	"sarPrm.py": All the parameters(electrical, mechanical, 		
	           geometrical) are defined here.
	"drawFigures.py": Useful drawing function are found here.
	"SARdata.py": Theoric Raw Data code is developed here.

How it works?
--------------
	1) Create an INPUT folder to save the input files. 
	2) Open the main script "Interferometria_main.py"
	   and specify the INPUT folder path in the variable 
	   "dir_input_files" located in the header of the 
	   script. 
	3) The OUTPUT files will be saved in the folder 
	   ~/Results/RawData_n, where 'n' is the number of 
	   files. 
