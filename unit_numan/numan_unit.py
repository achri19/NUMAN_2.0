# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 18:43:04 2024

@author: pbiswa5
"""
#------------------------------------------------------------------------------
# NUMAN 2.0 Model  
# This code is for "DETERMINISTIC" analysis for NUMAN 2.0
# Original Model: Chen, R., & Twilley, R. R. (1999). A simulation model of organic matter and nutrient accumulation in mangrove wetland soils. Biogeochemistry, 44, 93-118.
# Author: Biswas, P., Twilley, R. R., Rovai, A.S., Christensen, A., Kameshwar, S.
# Python Script: Biswas P. 
# Script: December 12, 2022
# Revision 1: May 30, 2023
# Revision 2: July 21, 2023
# Revision 3: September 10, 2023
# Revision 4: December 12, 2023
# Revision 5: February 07, 2024
# Revision 6: March 21, 2024
# Revision 7: June 06, 2024
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Importing necessary libraries 
# Make sure to have the following libraries/packages installed before running the code
import numpy as np
import pandas as pd
import scipy.integrate
import os
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Define the input Excel file containing parameters for all sites
# Constructing directory 
# Step 1: find the directory of .py first 
# Step 2: make this directory the main folder path 
# Step 3: construct the path to the csv
script_dirc = os.path.dirname(os.path.abspath(__file__))
main_folder = os.path.join(script_dirc)
input_file_path = os.path.join(main_folder, "input_parameter", "site_parameters.csv")
print(f"Constructed input file path: {input_file_path}")
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# reading the excel files and store as data frame 
parameter_df = pd.read_csv(input_file_path)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# The calculation refers to from year 0 (initial cohort deposition) to the 100 years' old cohort. 
# Discrete Counting i.e. 0,1,2,3,........99 (total 100 years)
# list1 stores different years data for each site in excel files for calculation
# list4 stores end product of each sites calculation like total accretion, ultimater soil properties up to the total soil accretion calculated or predicted 
list1 = [[] for _ in range(0, 100)]
list4=[[]]
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Iterating over each site and perform calculations according to the serial provided in the excel files in the site column
for index, site_params in parameter_df.iterrows():
    try:
        # parameter extraction 
        site_name = site_params['Site'] # extracts the site name as a categorical variable 
        si = site_params['si'] # inorganic matter loading rate (g/cm2/yr)
        b0 = site_params['b0'] # self packing density of pure organic matter (g/cm3)
        bi = site_params['bi'] # self packing density of pure inorganic matter (g/cm3)
        c0 = site_params['c0'] # lignin content in the surface deposit (Dry matter) in g/g
        c1 = site_params['c1'] # ash content in the root (in dry mass) in g/g
        c2 = site_params['c2'] # cellulose content in the surface deposit (dry matter) in g/g
        c3 = site_params['c3'] # cellulose fraction in the leaf litter (g/g)
        c4 = site_params['c4'] # cellulose content in the roots (dry mass) in g/g
        e = site_params['e'] # root attenuation, or root distribution rate (/cm)
        fc1 = site_params['fc1'] # lignin content in the fine root (dry mass) in g/g
        fc2 = site_params['fc2'] # lignin content in the coarse root (dry mass) in g/g
        fc3 = site_params['fc3'] # lignin content in the large root (dry mass) in g/g
        ka = site_params['ka'] # aerobic decomposition rate of labile organic matter or decomposition rate at surface cohort (/yr)
        kb = site_params['kb'] # anaerobic decomposition rate of labile organic matter (/yr)
        kc = site_params['kc'] # cellulose decomposition rate (/yr)
        cr = site_params['cr'] # Coarse root turnover rate (/yr)
        ke = site_params['ke'] # export rate of litters 
        kl = site_params['kl'] # lignin decomposition rate (/yr)
        km = site_params['km'] # large root turnover rate (/yr)
        kr = site_params['kr'] # fine root turnover rate (/yr)
        r0 = site_params['R0'] # root biomass at the surface (g/cm2)
        LP = site_params['LP'] # leaf litter deposition rate (g/cm2/yr)
        WP = site_params['WP'] # woody debris deposition rate (g/cm2/yr)
        
        
        # Calculation Part  
        
        # Net litter production NLP, g/cm2/yr (Leaf, reproductive organs and miscellaneous)
        # Our focus is the litter which retains in site after tidal export. 
        NLP= (LP*(1.0-ke)) 
        # assumed to be constant each year
        
        
        # Net wood accumulation, NWP g/cm2/yr 
        # our concern is how much is retained after tidal export
        NWP=(WP*(1-ke))
        
        
        # Inorganic matter on the surface cohort wi_0, g/cm2/yr assuming that the deficiency for above ground production to be taken up by the fibrous root in the surface cohort to keep constant production     
        wi_0=((si) - ((c1)*(ke)*(LP+WP)))
        
        #----------------------------------------------------------------------
        def upp_depth_approx(dt, db, r, e, p, wi, bi, c): # defining function 
            # r = root biomass at the surface, g/cm2
            # e = root attenuation, or distribution rate, /cm,
            # p = OM's bulk or self-packing density (g/cm3)
            # bi = self packing density of inorganic matter (g/cm3)
            # c = Ash content (fraction) in the DM of the cohort (g/g)
            udepth = db 
            count = 0
            diff = 1
            while diff > 0.000001:
                # r(t) = integration of (r0*exp(-eD)) from dt to db 
                # D is soil depth in cm 
                # Similar iteration approach to Morris & Bowden (1986), Chen & Twilley (1999) is applied 
                # Please refer to Tom Kaiser's python version of original NUMAN model 
                # diff threshold should be set up by the user depending on how precise estimates the user want to make
                integr = scipy.integrate.quad(lambda d: r * np.exp(-e * d), dt, udepth)[0]  # total root biomass in a cohort
                ldepth = udepth
                # previous db will be dt for next cohort. this function iterates to find the next db 
                # We separated the root mass organic section from the inorganic to make sure that any addition of organic mistakenly from inorganic will give rise to the volume cohort
                udepth = db + (integr * (1 - c) / p) + ((wi + integr * c) / bi)  # p=b0
                diff = np.abs(udepth - ldepth)  # this iteration will stop if abs(i(n+1)-i(n))<0.0000001
                count += 1  # this iteration cycle will continue until the difference is 0.000001 or less between the function evaluated and the iterated value of db itself
            return udepth  # Once it fulfills the criteria imposed, the loop will be terminated automatically and the data will be stored

        # this loop for finding the depth was adapted from Tom Kaiser's python script for NUMAN, which was modified later by Alexandra Christensen. 
        #----------------------------------------------------------------------
        
        # initializing calculation setting initial values to zero
        LOM = ROM = ROM_cell = ROM_lig = LOM_a = ROM_a = ROM_a_cell = ROM_a_lig = LOM_Nr = ROM_Nr = dt = db = frp  = TOM = V = db_corr = BD = OM_Per = Age = S_m_OM = S_m_OM_C = S_m_BD = S_m_C_Seq = S_m_OC_gPcc = S_m_IM_P = S_m_IM_Con = S_necromass = S_root_Nr = 0.0
        
        # initiate loop to start calculation 
        for i in range(0, 100):
            # the range corresponds to the tenure up to which the calculation has to be carried out
            if i==0:
                # Corresponds initial cohort <-- it actually initiates the first cohort like a building block to propagate through timeline
                
                # Labile organic matter (LOM) to the cohort initially, g/cm2/yr, corresponds to LOM(0)
                LOM= NLP*(1-(c0+c3+c1))+ NWP*(1-(c2+c1+c4))
                
                
                # Refractory organic matter --> ROM (0) is comprised of cellulose (ROM_cell(0)) & Lignin (ROM_Lig (0))
                
                # ROM_cell(0)
                ROM_cell = c3*NLP+ c4*NWP
                # ROM_lig(0)
                ROM_lig = c0*NLP +c2*NWP
                
                #Refractory organic matter to cohot initially corresponds to ROM(0) is sum of cellulose and lignin pools
                ROM=ROM_cell+ROM_lig
                
                
                # Above ground contribution to a specific cohort 
                # Above ground contribution to labile pool 
                LOM_a = LOM
                # Above ground contribution to refractory cellulose pool
                ROM_a_cell = ROM_cell
                # Above ground contribution to refractory lignin pool
                ROM_a_lig = ROM_lig
                # Above ground contribution to refractory pool is the sum of cellulose and lignin pool 
                ROM_a =ROM_a_cell+ROM_a_lig
                
                
                # Rest of the thing should be from the Necro mass of root
                # Necro root contribution to labile pool 
                LOM_Nr = LOM -LOM_a
                # Necro root contribution to refractory pool 
                ROM_Nr = ROM-ROM_a
               
                
            elif i==1:  
                # Next year <-- whatever mass was in the cohort goes through aerobic decomposition 
                # New necro root adds mass to labile and refractory pools, so the cohort evolves and requires updates according to the depth profile and time 
                
                # k=ka, for surface cohort aerobic decomposition
                
                # Labile fraction of root is (1- ash fraction - cellulose fraction -lignin fraction)
                
                LOM= list1[i-1][4] - ka*list1[i-1][4]+((1-(fc1+c4+c1)))*list1[i-1][11]+((1-(fc2+c4+c1)))*list1[i-1][12]+((1-(fc3+c4+c1)))*list1[i-1][13]#+list1[i-1][6]*kc.value*(1-f2.value)
                
                # Refractiory cellulose pool 
                ROM_cell = list1[i-1][5]-kc*list1[i-1][5]+(list1[i-1][11]*(c4)+list1[i-1][12]*(c4)+list1[i-1][13]*(c4))
                
                # Refractory lignin pool 
                ROM_lig= list1[i-1][6]-kl*list1[i-1][6]+(list1[i-1][11]*fc1+list1[i-1][12]*fc2+list1[i-1][13]*fc3) #ka.value*list1[i-1][5]*f3.value*(1-f2.value)
                
                # Total Refractory pool for the specific cohort 
                ROM=ROM_cell+ROM_lig
                
                #Above ground contribution can be calculated by just tracking the retaining  labile and refractoy pool as it was initially using the decay constants each steps
                
                # Remaining labile from above ground litter 
                LOM_a = list1[i-1][26] - ka*list1[i-1][26]
                
                # Remaining cellulose from the aboveground litter 
                ROM_a_cell= list1[i-1][27]-kc*list1[i-1][27] 
                
                # Remaining lignin from the aboveground litter 
                ROM_a_lig= list1[i-1][28]-kl*list1[i-1][28]
                
                # Total remaining refractory from aboveground litter (above ground litter residue)
                ROM_a =ROM_a_cell+ROM_a_lig
                
                # Necro mass of root is the the remaining other than the aboveground remaining contribution after decay
                
                # Labile Necro root mass 
                LOM_Nr = LOM -LOM_a
                
                # Refractory necro root mass 
                ROM_Nr = ROM-ROM_a
                
                
                
            else:
                
                # k=kb, for surface cohort anaerobic decomposition
                
                # Labile fraction of root is (1- ash fraction - cellulose fraction -lignin fraction)
                
                LOM= list1[i-1][4] - kb*list1[i-1][4]+((1-(fc1+c4+c1)))*list1[i-1][11]+((1-(fc2+c4+c1)))*list1[i-1][12]+((1-(fc3+c4+c1)))*list1[i-1][13]#+list1[i-1][6]*kc.value*(1-f2.value)
                
                # Refractiory cellulose pool 
                ROM_cell = list1[i-1][5]-kc*list1[i-1][5]+(list1[i-1][11]*(c4)+list1[i-1][12]*(c4)+list1[i-1][13]*(c4))
                
                # Refractory lignin pool 
                ROM_lig= list1[i-1][6]-kl*list1[i-1][6]+(list1[i-1][11]*fc1+list1[i-1][12]*fc2+list1[i-1][13]*fc3) #ka.value*list1[i-1][5]*f3.value*(1-f2.value)
                
                # Total Refractory pool for the specific cohort 
                ROM=ROM_cell+ROM_lig
                
                #Above ground contribution can be calculated by just tracking the retaining  labile and refractoy pool as it was initially using the decay constants each steps
                
                # Remaining labile from above ground litter 
                LOM_a = list1[i-1][26] - kb*list1[i-1][26]
                
                # Remaining cellulose from the aboveground litter 
                ROM_a_cell= list1[i-1][27]-kc*list1[i-1][27] 
                
                # Remaining lignin from the aboveground litter 
                ROM_a_lig= list1[i-1][28]-kl*list1[i-1][28]
                
                # Total remaining refractory from aboveground thing
                ROM_a =ROM_a_cell+ROM_a_lig
                
                #Necro mass of root is the the remaining other than the aboveground remaining contribution after decay
                
                # Labile Necro root mass 
                LOM_Nr = LOM -LOM_a
                
                # Refractory necro root mass 
                ROM_Nr = ROM-ROM_a
            
            
            # Cohort's lower depth calculation 
            db=upp_depth_approx(dt,dt + ((LOM + ROM) /b0),r0,e,b0,wi_0,bi,c1)
            #db calls the derived function to perform the iteration on the basis of the data available 
            
            #Root biomass in a cohort at any specific year 
            rt=(scipy.integrate.quad(lambda d: r0 * np.exp(-e * d),dt,db)[0])
            
            # Root turnover is already set in such a way so that it captures relative contribution of different class size in the biomass stock
            # fine root production = Total root x fine root turnover 
            frp=rt*kr
            
            # Coarse root production = Total root x Course root turnover 
            crp = rt*cr
            
            # large root production = Total root x large root turnover 
            lrp=rt*km
            
            #Total Organic Matter
            TOM=(LOM+ROM+rt*(1.0-c1))
            
            # inorganic matter in a specific cohort specific time 
            wi_t= wi_0+rt*c1
            
            # volume of the cohort per unit surface area
            V=((TOM)/b0)+(wi_t/bi) 
            
            # corrected depth 
            db_corr=dt+V 
            
            # bulk density 
            BD= (TOM+wi_t)/(db_corr-dt) 
            m_BD= BD *(db_corr-dt) # bulk density x cohort volume or height 
            S_m_BD = S_m_BD+m_BD # Arithmatic sum of (Bulk Density x cohort volume or height)
            
            #ASH FREE DRY WEIGHT (AFDW) & Organic Matter % (OM%) are same 
            OM_Per= (TOM/(TOM+wi_t)) * 100.00
            m_OM = OM_Per * (db_corr-dt) # OM% x cohort volume or height 
            S_m_OM = S_m_OM+m_OM # sum of (OM%x cohort volume or height) 
            
            #OC is the fraction organic matter 
            # Holmquist et al. (2017)
            OC = 0.074*((OM_Per/100)**2) +0.421*(OM_Per/100)-0.0080 
            
            
            #Organic matter content in g/cc
            OC_gPcc = BD*OC 
            m_OC_gPcc = OC_gPcc*(db_corr-dt)*1000 #1000 is to convert it to mg # OC x cohort volume or height 
            S_m_OC_gPcc= S_m_OC_gPcc + m_OC_gPcc # Sun of (OC X cohort volume or height) 
            
            
            #carbon sequestration for each cohort 
            C_Seq=(db_corr-dt)*OC_gPcc*(10**4) # g/m2 for every year's cohort 
            m_C_Seq= C_Seq*(db_corr-dt) # C_seq x cohort volume or height 
            S_m_C_Seq=S_m_C_Seq+m_C_Seq #sum of (C_seq x cohort volume or height)
            
            #Organic matter concentration
            OM_C= BD*OM_Per*1000/100# organic matter content, mg/cm3
            m_OM_C=OM_C*(db_corr-dt) # OM concentration x cohort volume or height 
            S_m_OM_C=S_m_OM_C+m_OM_C # Sum of (OM concentration x cohort volume or height)
            
            #inorganic content%
            IM_P =(100-OM_Per) #%
            m_IM_P=IM_P*(db_corr-dt) # Inorganic matter% x cohort volume or height 
            S_m_IM_P = S_m_IM_P + m_IM_P # sum of (Inorganic matter% x cohort volume or height)
            
            #Inorganic content mg/m3
            IM_Con = (IM_P*BD)*1000/100 # In organic matter concentration (converted to mg/cm3)
            m_IM_Con=IM_Con*(db_corr-dt) # Inorganic matter concentration x cohort volume or height 
            S_m_IM_Con = S_m_IM_Con + m_IM_Con # sum of (Inorganic matter concentration x cohort volume or height )
            
            # Necro mass:
            # Total necromass
            necro_mass = LOM+ROM # necromass of a specific cohort
            S_necromass = S_necromass+necro_mass # sum of necromass up to certain depth after certain year simulation 
            
            
            # Root necromass
            Root_Nr = LOM_Nr + ROM_Nr # root necromass of a specific cohort
            S_root_Nr = S_root_Nr +Root_Nr # sum of root necromass up to certain depth after certain year simulation 
            
            #volume contribution for 1 cm2 area
            Va=(LOM_a+ROM_a)/b0 # aboveground litter contribution to cohort
            Vr=(rt*(1-c1))/b0 # root contribution organic portion 
            VNr = (LOM_Nr+ROM_Nr)/b0 # Necromass contribution 
            Vi=wi_t/bi # inorganic contribution
            
            #Volume contribution percentage
            # Volume contribution % = (specific element's contribution in a cohort / total volume of the cohort for all elements) x 100%
            Va_P = Va*100/V
            Vr_P = Vr*100/V
            VNr_P = VNr*100/V
            Vi_P = Vi*100/V
            
            
            Age= i
        
            # concept of array or matrices has been used for calculation.
            # i-1 corresponds to the immediate previous year's value
            # 4,5,6....... etc. are corresponding to the specific property as I made entry.
            # the index starts from 0, not from 1.
            # e.g.: Like LOM is [1], ROM is [2]
        
            list2 = []
        
            list2 += [Age,LP,NLP,NWP,LOM,ROM_cell,ROM_lig,ROM,dt, db,rt,frp,crp,lrp,TOM,wi_t,V,db_corr,BD,OM_Per,OC_gPcc*1000,C_Seq,OM_C,IM_P,IM_Con,necro_mass*100, LOM_a, ROM_a_cell, ROM_a_lig, ROM_a, LOM_Nr,ROM_Nr,S_root_Nr*100,Va, Vr,VNr,Vi,Va_P, Vr_P,VNr_P,Vi_P]

            dt = db
            
            # previous year's Db set as new year's dt
            
            list1[i] = list2  # stores data
        
        # stores individual files in the output folder
        output_folder = os.path.join(main_folder, "output_files", "individual_site")
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Define the output file path & naming convention
        output_file_path = os.path.join(output_folder, f"site {site_name} output.csv")
        
        # Convert list1 to a pandas DataFrame
        df = pd.DataFrame(list1, columns=[
            'Age(yr)','LP','NLP','NWP','LOM (g/cm2)','ROM_cellulose (g/cm2)','ROM_lignin(g/cm2)','ROM (g/cm2)','dt (cm)','db (cm)','rt (g/cm2)','frp (g/cm2)','crp(g/cm2)','lrp(g/cm2)','TOM (g/cm2)','wi_t (g/cm2)','V(cm3)','db_corr (cm)','BD (g/cm3)','AFDW or OM %','OC Density (mg/cm3)','C_Seq (g/m2/yr)','OM Density (mg/cm3)','Inorganic Matter(%)', 'Inorganic Matter Density (mg/cm3)','Necro_Mass (Mg/ha)','LOM_A (g/cm2)','ROM_A_cell (g/cm2)','ROM_A_lig(g/cm2)','ROM_A (g/cm2)','LOM_Nr (g/cm2)','ROM_Nr (g/cm2)','Cumulative Nr (Mg/ha)','Vol cont to cohort from litter','Vol cont to cohort from live root','Vol cont to cohort from necro root','Vol cont to cohort from inorganic','Vol cont % to cohort from litter','Vol cont % to cohort from live root','Vol cont %to cohort from necro root','Vol cont % to cohort from inorganic'
        ])
        df.drop(['LP', 'NLP', 'NWP'], axis=1, inplace=True)
        # Save the DataFrame to a CSV file
        df.to_csv(output_file_path, index=False)
        
        print(f"site {site_name} completed")
        
        list4.append([site_name,list1[-1][17],list1[-1][17]/100,S_m_OM/list1[-1][17],S_m_OM_C/list1[-1][17],S_m_BD/list1[-1][17],S_m_C_Seq/list1[-1][17],S_m_OC_gPcc/list1[-1][17],S_m_IM_P/list1[-1][17],S_m_IM_Con/list1[-1][17],S_necromass*100,S_root_Nr*100,(list1[-1][32]+((list1[-1][32]-list1[-2][32])*(40-list1[-1][17])/(list1[-1][17]-list1[-2][17])))])
        
    except Exception as e:
        print(f"An error occurred during calculation for {site_name}: {e}")
list4 = list4[1:]
# Convert list4 to a pandas DataFrame
summary_df = pd.DataFrame(list4, columns=['site', 'Total Accretion (cm)', 'Mean Accretion (cm/yr)', 'Mean OM%', 'Mean OM Density (mg/cm3)', 'Mean BD (g/cm3)', 'C_Seq (g/m2/yr)', 'OC Density (mg/cm3)', 'IM(%)', 'IM Density (mg/cm3)', 'Necromass (Mg/ha)', 'Root Necromass (Mg/ha)', 'Root Necromass Projected Up to 40 cm (Mg/ha)'])

# Define the output folder
output_folder = os.path.join(main_folder, "output_files", "summary")
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Save the DataFrame to a CSV file
output_file_path = os.path.join(output_folder, "summary_output.csv")
summary_df.to_csv(output_file_path, index=False)

print(parameter_df.columns)
