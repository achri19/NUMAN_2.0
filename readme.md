# NUMAR 2.0
## Eco-geosoil Model 
The Nutrient Mangrove (NUMAN) model, initially developed by Chen & Twilley (1999) was designed to simulate soil dynamics in mangrove ecosystems. This pioneering work laid the groundwork for further research in the field. We have reformulated the NUMAN model, testing it across various sites with different nutrient gradients. Key updates include revising the equations for a more accurate mass balance and integrating recent field data, resulting in a more robust and reliable model. Now, NUMAN 2.0 can confidently track necromass from various sources, carbon density in the soil, carbon sequestration rates, and the relative volume contributions from different sources or tributaries, along with a wide range of other output variables related to soil properties. To ensure its generality, the model was tested in three distinct mangrove settings in South Florida, each representing different nutrient gradients. Additionally, this version incorporates both probabilistic and deterministic analyses to enhance its robustness.

Original paper for previous NUMAN version: Chen, R., & Twilley, R. R. (1999). A simulation model of organic matter and nutrient accumulation in mangrove wetland soils. Biogeochemistry, 44(1), 93–118. https://doi.org/10.1007/BF00993000

Related paper/manuscript for NUMAN 2.0: IN REVIEW *(The link and citation will be shared once the preprint is available to share)*
*Read the manuscript and supplementary material for more detailed understanding of the model*

## Contributors to the model
- Pradipta Biswas<sup><1>
- Robert R. Twilley
- André S. Rovai1,*,
- Alexandra Christensen2,
- Zoe I. Shribman1,**,
- Sabarethinam Kameshwar3

1 Department of Oceanography and Coastal Sciences, College of Coast and Environment, Louisiana State University, Baton Rouge, Louisiana, USA

2 Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, 91109, USA

3 Department of Civil & Environmental Engineering, Louisiana State University, Baton Rouge, Louisiana, USA
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
