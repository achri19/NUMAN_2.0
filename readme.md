# NUMAR 2.0

# --------------------------------------------------------------
**Packages Required to Install**

Make sure the user’s working environment has following packages 
- numpy
- scipy
- pandas
- spyder
- scikit-learn


The script was developed in spyder 5.5.1. python IDE 3.12.2 in anaconda navigator environment. Download and set-up instructions can be found in the respective website (https://www.anaconda.com/products/individual). It is not recommended to change the location or name of the input CSV file as the directory is set up in a certain way to access the data input. Any change to the input should be made at the user’s own discretion. 

# --------------------------------------------------------------
General Guide
Here are some general guidelines to run the model.  
* For deterministic analysis, the folder “unit_numan” should be used. Inside this folder, there are subfolders to nest the input and output files. Inside the subfolder site_parameters.csv is located, where any sites can be modified or added according to requirements. After running “numan_unit.py in the same “unit_numan” folder, the output csv gets stored in the “output_files” subfolder. User can inspect individual sites at subfolder “individual_sites” and all sites’ summary at “summary” subfolder. 

* For probabilistic analysis, access the “uncertainty_numan” folder. The folder has the same organization of input and output folder nested inside. To run the model, please open “random_numan.py” and run. If the user wants to induce uncertainties with more input variables, then an additional column in “site_parameters_random.csv” needs to be added with proper notation relating to the uncertainty. Then additional lines to be added in between line #102 ~138, which involves the importation of uncertainty parameters, modifying the type of distribution. Associated changes might need to be adjusted according to the change the user aims to make with the output. Outputs will be stored in “output_files”. To inspect any site-specific csv, go to the “individual_sites” subfolder, inside which, the site-specific analysis file is stored under their name specific nested folder. All simulation summery can be found inside “summery” subfolder. 

* Input and output variables with their unit has been discussed in the later section.  
# --------------------------------------------------------------
