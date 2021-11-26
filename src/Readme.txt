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
	"main.py"

	Dependencies:
	-------------
	"imagingAlgorithm.py": The FDBP Algorithm is implemented here.
	"sarPrm.py": All the parameters(electrical, mechanical, 		
	           geometrical) are defined here.
	"drawFigures.py": Useful drawing function are found here.
	*"sarData.py": Theoric Raw Data code is developed here. Just include it 
		      if you need to generate theoric data set. 

How it works?
--------------
~ : '/media/soporte/e2a2d167-bcfd-400a-91c8-f1236df2f7e4/soporte/Landslide_Project/
     Desarrollo_v2/Software/Procesamiento/'

	1) Create a new folder "new_folder", which contains the raw data, and save it
           in the following path:
	                           ~/DataSet/new_folder/

	2) Modify the script "main.py", by specifying the raw data folder path in the 
           variable "dir_input_files".

	3) The OUTPUT files will be saved in the folder (created automatically): 
	                            ~/Results/new_folder/


