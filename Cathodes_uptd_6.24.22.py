# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 13:51:29 2022

@author: donmc
"""
from simmate.shortcuts import setup
from simmate.database.third_parties import MatProjStructure, JarvisStructure, CodStructure, OqmdStructure
import pandas as pd

#%% ONLY DO ONCE! LOAD IN ALL MATPROJ STRUCTURES & CALC ENERGY COLUMNS
#MatProjStructure.load_remote_archive()
#MatProjStructure.update_all_stabilities()

#%% MATPROJ FILTER TOOLKIT

# FILTER MATPROJ DATABASE TO INCLUDE F STRUCTURES ON HULL
data_mp_list = MatProjStructure.objects.filter(
    elements__icontains='"F"', 
    #elements__icontains='"F"',
    is_stable=True,
    nelements__gt="1",
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_toolkit()


#%% MATPROJ FILTER DATAFRAME
data_mp_table = MatProjStructure.objects.filter(
    elements__icontains='"F"', 
    #elements__icontains='"F"',
    is_stable=True,
    nelements__gt="1",
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_dataframe()

#%% FILTER JARVIS AND COD

# FILTER JARVIS DATABASE TO INCLUDE F STRUCTURES ON HULL
data_jar = JarvisStructure.objects.filter(
    elements__icontains='"F"', 
    #elements__icontains='"F"',
    #is_stable=True,
    energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_dataframe()

# COD HAS NO FINAL ENERGY OR ENERGY ABOVE HULL
data_cod = CodStructure.objects.filter(
    elements__icontains='"F"', 
    #elements__icontains='"F"',
    #is_stable=True,
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_dataframe()
#%% DEF DATAFRAME

#Examples.objects.filter(...).filter(...)
#to_toolkit() --> gives list of structure objects


#for structure in structures:
 #   structure.remove_species("F") 
  #  MatProjStructure.objects.filter(chemical_system=structure.composition.chemical_system).to_dataframe()


# MAKE DATAFRAME WITH DEFLUORINATED STRUCTURES: DATA MUST BE IN TOOKLIT FORM TO USE THIS
print("Len of data_mp_list:", len(data_mp_list))
deF_table2 = pd.DataFrame()
mp_index = 0
for struc in data_mp_list:
    struc.remove_species("F") 
    s = MatProjStructure.objects.filter(chemical_system=struc.composition.chemical_system, 
    is_stable=True).to_dataframe()
    #print("#1: data_mp index:", mp_index, ": # of de-F struc", len(s))
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            #print(len(s))
            new_row = pd.Series([s.structure_string[i], s.nsites[i],s.nelements[i],s.elements[i], s.chemical_system[i], s.density[i], s.density_atomic[i], s.volume[i], s.volume_molar[i], s.formula_full[i], s.formula_reduced[i], s.formula_anonymous[i], s.spacegroup[i], s.energy[i], s.energy_per_atom[i],s.energy_above_hull[i],s.is_stable[i],s.decomposes_to[i],s.formation_energy[i],s.formation_energy_per_atom[i],s.id[i]], 
                                index=[s.columns[0], s.columns[1], s.columns[2], s.columns[3], s.columns[4], s.columns[5], s.columns[6], s.columns[7], s.columns[8], s.columns[9], s.columns[10], s.columns[11],s.columns[12], s.columns[13], s.columns[14], s.columns[15],s.columns[16], s.columns[17], s.columns[18], s.columns[19],s.columns[20]])
            deF_table2 = pd.concat([deF_table2, new_row.to_frame().T], ignore_index=True)
            i = i+1
            #print("len of dataframe:",len(deflo.index), "i value:", i)
            if i == len(s):
                break
        #print("done1")
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1
    
#%% DEF LIST

# MAKE LIST/TOOLKIT WITH DEFLUORINATED STRUCTURES: DATA MUST BE IN TOOKLIT FORM TO USE THIS
print("Len of data_mp_list:", len(data_mp_list))
mp_index = 0
deF_list2 = []
for struc in data_mp_list:
    struc.remove_species("F") 
    s = MatProjStructure.objects.filter(chemical_system=struc.composition.chemical_system, 
    is_stable=True).to_toolkit()
    #print("#1: data_mp index:", mp_index, ": # of de-F struc", len(s))
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            #print(len(s))
            deF_list2.append(x)
            #i = i+1
            #print("len of dataframe:",len(deflo.index), "i value:", i)
            #if i == len(s):
             #   break
        #print("done1")
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1


#%% DEF DATAFRAME W ATOM #'s 

# MAKE DATAFRAME WITH DEFLUORINATED STRUCTURES: DATA MUST BE IN TOOKLIT FORM TO USE THIS
print("Len of data_mp_list:", len(data_mp_list))
deF_table = pd.DataFrame()
mp_index = 0
for struc in data_mp_list:
    struc.remove_species("F") 
    s = MatProjStructure.objects.filter(formula_full=struc.composition.formula, 
    is_stable=True).to_dataframe()
    #print("#1: data_mp index:", mp_index, ": # of de-F struc", len(s))
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            #print(len(s))
            new_row = pd.Series([s.structure_string[i], s.nsites[i],s.nelements[i],s.elements[i], s.chemical_system[i], s.density[i], s.density_atomic[i], s.volume[i], s.volume_molar[i], s.formula_full[i], s.formula_reduced[i], s.formula_anonymous[i], s.spacegroup[i], s.energy[i], s.energy_per_atom[i],s.energy_above_hull[i],s.is_stable[i],s.decomposes_to[i],s.formation_energy[i],s.formation_energy_per_atom[i],s.id[i]], 
                                index=[s.columns[0], s.columns[1], s.columns[2], s.columns[3], s.columns[4], s.columns[5], s.columns[6], s.columns[7], s.columns[8], s.columns[9], s.columns[10], s.columns[11],s.columns[12], s.columns[13], s.columns[14], s.columns[15],s.columns[16], s.columns[17], s.columns[18], s.columns[19],s.columns[20]])
            deF_table = pd.concat([deF_table, new_row.to_frame().T], ignore_index=True)
            i = i+1
            #print("len of dataframe:",len(deflo.index), "i value:", i)
            if i == len(s):
                break
        #print("done1")
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1

#%% DEF LIST W ATOM #'s

# MAKE LIST/TOOLKIT WITH DEFLUORINATED STRUCTURES: DATA MUST BE IN TOOKLIT FORM TO USE THIS
print("Len of data_mp_list:", len(data_mp_list))
mp_index = 0
deF_list = []
for struc in data_mp_list:
    struc.remove_species("F") 
    s = MatProjStructure.objects.filter(formula_full=struc.composition.formula, 
    is_stable=True).to_toolkit()
    #print("#1: data_mp index:", mp_index, ": # of de-F struc", len(s))
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            #print(len(s))
            deF_list.append(x)
            #i = i+1
            #print("len of dataframe:",len(deflo.index), "i value:", i)
            #if i == len(s):
             #   break
        #print("done1")
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index + 1



#%% DEF DATAFRAME W ATOM #'S, GET DELTA VOLUME TOO

# 5/4/22 MAKE DATAFRAME WITH DEFLUORINATED STRUCTURES: DATA MUST BE IN TOOKLIT FORM TO USE THIS
print("Len of data_mp_list:", len(data_mp_list))
deF_table2 = pd.DataFrame()
deltav_table2 = pd.DataFrame()
mp_index = 0
for struc in data_mp_list:
    mv = data_mp_table.volume_molar[mp_index]
    ff = data_mp_table.formula_full[mp_index]
    f_id = data_mp_table.id[mp_index]
    f_energy = data_mp_table.energy[mp_index]
    f_energyperatom = data_mp_table.energy_per_atom[mp_index]
    f_nsites = data_mp_table.nsites[mp_index]
    struc.remove_species("F")
    s = MatProjStructure.objects.filter(formula_full=struc.composition.formula, 
    is_stable=True).to_dataframe()
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            deltav = round(mv-s.volume_molar[i],2)
            print(ff, "volume diff:", deltav)
            new_row = pd.Series([s.structure_string[i], s.nsites[i],s.nelements[i],s.elements[i], s.chemical_system[i], s.density[i], s.density_atomic[i], s.volume[i], s.volume_molar[i], s.formula_full[i], s.formula_reduced[i], s.formula_anonymous[i], s.spacegroup[i], s.energy[i], s.energy_per_atom[i],s.energy_above_hull[i],s.is_stable[i],s.decomposes_to[i],s.formation_energy[i],s.formation_energy_per_atom[i],s.id[i]], 
                                index=[s.columns[0], s.columns[1], s.columns[2], s.columns[3], s.columns[4], s.columns[5], s.columns[6], s.columns[7], s.columns[8], s.columns[9], s.columns[10], s.columns[11],s.columns[12], s.columns[13], s.columns[14], s.columns[15],s.columns[16], s.columns[17], s.columns[18], s.columns[19],s.columns[20]])
            deF_table2 = pd.concat([deF_table2, new_row.to_frame().T], ignore_index=True)
            
            # Add row to dataframe with fluorinated struc and difference in volume 
            new_deltav_row = pd.Series([ff, deltav, f_id, s.id[i], f_energy, s.energy[i], f_energyperatom, s.energy_per_atom[i], f_nsites, s.nsites[i]], index=['Formula', 'Delta V', 'F id', 'De-F id', 'F energy','De-F energy','F energy per atom','De-F energy per atom','F nsites','De-F nsites'])
            deltav_table2 = pd.concat([deltav_table2, new_deltav_row.to_frame().T], ignore_index=True)
            i = i+1
            if i == len(s):
                break
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1


#%% MAKE COPY OF DEF LIST/TABLE
deF_list_cp = deF_list
deF_table_cp = deF_table
deltav_table_cp = deltav_table
#deF_list2_cp = deF_list2
#deF_table2_cp = deF_table2
    
#%% CELL VOLTAGE VS LI/LI+  
li_energy = -1.909
liF_energy = -9.6903
e = 0
voltage_table = pd.DataFrame()
#print(deltav_table.["F nsites"].values[1])
print(len(deltav_table2))
while e<len(deltav_table2):
    formula = deltav_table2.loc[e][0]
    deF_en = deltav_table2.loc[e][5]
    F_en = deltav_table2.loc[e][4]
    F_transfer = deltav_table2.loc[e][8] - deltav_table2.loc[e][9]
    deltaG = ((deF_en+(F_transfer*liF_energy)) - (F_en+(F_transfer*li_energy)))
    voltage = deltaG/(-F_transfer)
    
    new_row = pd.Series([formula, voltage], index=['Formula', 'Voltage'])
    voltage_table = pd.concat([voltage_table, new_row.to_frame().T], ignore_index=True)
    e = e+1

print(type(deF_en))
print(type(F_transfer))
print(type(liF_energy))
print("F_en",type(F_en))
print(type(li_energy))
    
#%%
for i in deF_table_cp.formula_full:
    #print(i)
    #if data_mp_table['formula_full'].str.contains(i).any():   
    filter[data_mp_table['formula_full'].str.contains(i)]    
        #print(data_mp_table.formula_full[])
    print(filter)
        


if deF_table_cp['formula_full'].str.contains('Li1 Au1').any():
    print("yes")
    
if data_mp_table['formula_full'].str.contains('Rh1').any():
    print("yes again")


#%% KBI TESTING 5/4/22

kbif_table = MatProjStructure.objects.filter(
    #elements__icontains='"Bi"', 
    #elements__icontains='"Bi","F"',
    is_stable=True,
    id="mp-557724"
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_dataframe()

kbif_list = MatProjStructure.objects.filter(
    #elements__icontains='"Bi"', 
    #elements__icontains='"Bi","F"',
    is_stable=True,
    id="mp-557724"
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_toolkit()

kbi_table = pd.DataFrame()
mp_index = 0
for struc in kbif_list:
    struc.remove_species("F")
    s = MatProjStructure.objects.filter(formula_full=struc.composition.formula, 
    is_stable=True).to_dataframe()
    print(struc.composition.formula)
    if len(s)>0:
        i = 0
        print("# of matches:",len(s))
        for x in s:
            new_row = pd.Series([s.structure_string[i], s.nsites[i],s.nelements[i],s.elements[i], s.chemical_system[i], s.density[i], s.density_atomic[i], s.volume[i], s.volume_molar[i], s.formula_full[i], s.formula_reduced[i], s.formula_anonymous[i], s.spacegroup[i], s.energy[i], s.energy_per_atom[i],s.energy_above_hull[i],s.is_stable[i],s.decomposes_to[i],s.formation_energy[i],s.formation_energy_per_atom[i],s.id[i]], 
                                index=[s.columns[0], s.columns[1], s.columns[2], s.columns[3], s.columns[4], s.columns[5], s.columns[6], s.columns[7], s.columns[8], s.columns[9], s.columns[10], s.columns[11],s.columns[12], s.columns[13], s.columns[14], s.columns[15],s.columns[16], s.columns[17], s.columns[18], s.columns[19],s.columns[20]])
            kbi_table = pd.concat([kbi_table, new_row.to_frame().T], ignore_index=True)
            
            i = i+1
            if i == len(s):
                break
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1



#%% 5/4/22 FIND KBI FROM KBIF6




kbif = MatProjStructure.objects.filter(
    #formula_full="K Bi", 
    formula_reduced="KBiF6",
    #id="mp-557724",
    is_stable=True).to_dataframe()

#kbiDEf = kbif[0].remove_species("F")

kbi = MatProjStructure.objects.filter(
    #formula_full="K Bi", 
    id="mp-31104",
    is_stable=True).to_dataframe()

print(kbif.formula_full[0])
#print(kbiDEf)
print(kbi.formula_full[0])
print(kbif.formula_reduced[0])
print(kbi.formula_reduced[0])
#print(kbi)


#%% TABLE 3 CELL 1: 6/13/22 test of struc.composition.reduced_formula instead of struc.composition.formula

print("Len of data_mp_list:", len(data_mp_list))
deF_table3 = pd.DataFrame()
deltav_table3 = pd.DataFrame()
mp_index = 0
for struc in data_mp_list:
    mv = data_mp_table.volume_molar[mp_index]
    ff = data_mp_table.formula_full[mp_index]
    rf = data_mp_table.formula_reduced[mp_index]
    f_id = data_mp_table.id[mp_index]
    f_energy = data_mp_table.energy[mp_index]
    f_energyperatom = data_mp_table.energy_per_atom[mp_index]
    f_nsites = data_mp_table.nsites[mp_index]
    struc_deF = struc
    struc_deF.remove_species("F")   # Remove F from each entry in data_mp_list
    s = MatProjStructure.objects.filter(formula_reduced=struc_deF.composition.reduced_formula, # For each entry
    is_stable=True).to_dataframe()          # check if there is another MP entry with reduced formula (formula_reduced)
                                            # identical to the reduced fromula after removing F (struc.composition.reduced_formula)
    if len(s)>0:
        i = 0
        print("chem system= ", s.formula_reduced)
        print("formula= ", struc_deF.composition.reduced_formula)
        print("# of matches:",len(s))
        for x in s:
            deltav = round(mv-s.volume_molar[i],2)
            #print(ff, "volume diff:", deltav)
            
            #Make the deF_table3 table with new rows
            new_row = pd.Series([s.structure_string[i], s.nsites[i],s.nelements[i],s.elements[i], s.chemical_system[i], s.density[i], s.density_atomic[i], s.volume[i], s.volume_molar[i], s.formula_full[i], s.formula_reduced[i], s.formula_anonymous[i], s.spacegroup[i], s.energy[i], s.energy_per_atom[i],s.energy_above_hull[i],s.is_stable[i],s.decomposes_to[i],s.formation_energy[i],s.formation_energy_per_atom[i],s.id[i]], 
                                index=[s.columns[0], s.columns[1], s.columns[2], s.columns[3], s.columns[4], s.columns[5], s.columns[6], s.columns[7], s.columns[8], s.columns[9], s.columns[10], s.columns[11],s.columns[12], s.columns[13], s.columns[14], s.columns[15],s.columns[16], s.columns[17], s.columns[18], s.columns[19],s.columns[20]])
            deF_table3 = pd.concat([deF_table3, new_row.to_frame().T], ignore_index=True)
            
            # Make the deltav_table3 table with new rows 
            new_deltav_row = pd.Series([ff, rf, deltav, f_id, s.id[i], f_energy, s.energy[i], f_energyperatom, s.energy_per_atom[i], f_nsites, s.nsites[i]], index=['Formula', 'Reduced Formula','Delta V', 'F id', 'De-F id', 'F energy','De-F energy','F energy per atom','De-F energy per atom','F nsites','De-F nsites'])
            deltav_table3 = pd.concat([deltav_table3, new_deltav_row.to_frame().T], ignore_index=True)
            i = i+1
            if i == len(s):
                break
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1

#%% TABLE 3 CELL 2: 6/13/22 CELL VOLTAGE FOR 323 STRUCTURES
li_energy = -1.909
liF_energy = -9.6903
e = 0
voltage_table3 = pd.DataFrame()
#print(deltav_table.["F nsites"].values[1])
print(len(deltav_table3))
while e<len(deltav_table3):
    formula = deltav_table3.loc[e][0]
    deF_energy = deltav_table3.loc[e][6]
    F_energy = deltav_table3.loc[e][5]
    deF_energy_per_atom = deltav_table3.loc[e][8]
    F_energy_per_atom = deltav_table3.loc[e][7]
    # START HERE FOR NEW F_TRANSFER CODE
    F_reduced_form = deltav_table3.loc[e][1]       #reduced formula fluorinated
    deF_reduced_form = deF_table3.loc[e][10]       #reduced formual defluorinated
    
    F_capital = sum(1 for c in F_reduced_form if c.isupper())

    for x in F_reduced_form:     # take reduced formula and replace all letters with spaces so that                   
        if x.isdigit()==False:      # only numbers are left
            F_reduced_form = F_reduced_form.replace(x, " ")
        

    F_reduced_form_split = F_reduced_form.split() # make workable string list of #'s from previous "for" loop
    
    F_sites_list = []
    total=0
    for c in F_reduced_form_split:   # convert string list to integer list
        integer = int(c)
        F_sites_list.append(integer) 
        #total = total + integer

    if F_capital!=len(F_sites_list):
        F_total_sites = F_capital-len(F_sites_list)+sum(F_sites_list)
    else:
        F_total_sites = sum(F_sites_list) 
    
    #print(formula,";", ":total", F_total_sites)
    
    #START DEF 
    deF_capital = sum(1 for c in deF_reduced_form if c.isupper())

    for x in deF_reduced_form:
        if x.isdigit()==False:
            deF_reduced_form = deF_reduced_form.replace(x, " ")

    deF_reduced_form_split = deF_reduced_form.split() #make list of #'s from formula string
     
    deF_sites_list = []
    total=0
    for c in deF_reduced_form_split:               # convert string list to integer list
        integer = int(c)
        deF_sites_list.append(integer) 
        #total = total + integer

    if deF_capital!=len(deF_sites_list):
        deF_total_sites = deF_capital-len(deF_sites_list)+sum(deF_sites_list)
    else:
        deF_total_sites = sum(deF_sites_list)
            
    #print("DEF",deF_table3.loc[e][10],":", deF_total)
    
    F_transfer = F_total_sites - deF_total_sites
    
    print(deltav_table3.loc[e][1], ";",F_transfer)
    
    # END NEW F_TRANSFER CODE
    
    
    #deltaG = ((deF_energy+(F_transfer*liF_energy)) - (F_energy+(F_transfer*li_energy)))
    
    deltaG = ((deF_total_sites*deF_energy_per_atom+(F_transfer*liF_energy)) - 
              (F_total_sites*F_energy_per_atom+(F_transfer*li_energy)))
    voltage = deltaG/(-F_transfer)
    
    new_row = pd.Series([formula, voltage], index=['Formula', 'Voltage'])
    voltage_table3 = pd.concat([voltage_table3, new_row.to_frame().T], ignore_index=True)
    e = e+1



#%% PARTIAL 1 TABLE 4: 6/13/22 --> MAKE DELTAV_TABLE4 WHICH HAS ALL POTENTIAL PARTIAL DEF'S AND THE ORIGINAL ENTRY

data_mp_list = MatProjStructure.objects.filter(
    elements__icontains='"F"', 
    #elements__icontains='"F"',
    is_stable=True,
    nelements__gt="1",
    #energy_above_hull__lte=0.01,
    #formula_reduced="LuSF"
    ).to_toolkit()

print("Len of data_mp_list:", len(data_mp_list))
deF_table4 = pd.DataFrame()
deltav_table4 = pd.DataFrame()
mp_index = 0
for struc in data_mp_list:
    mv = data_mp_table.volume_molar[mp_index]
    ff = data_mp_table.formula_full[mp_index]
    rf = data_mp_table.formula_reduced[mp_index]
    f_id = data_mp_table.id[mp_index]
    f_energy = data_mp_table.energy[mp_index]
    f_energyperatom = data_mp_table.energy_per_atom[mp_index]
    f_nsites = data_mp_table.nsites[mp_index]
    
    # How to deal with partial defluorination? Likely problems:
        # 1. Identify entry that contains F but less than a reference entry: count # of total atoms
            # struc.composition.num_atoms gives number of atoms in entry
        # 2. If entry has less atoms than reference entry, make sure the "missing" atoms are F and not other atoms: count # of elements
            # len(struc.composition.reduced_composition) gives # of different elements in entry
        # 3. Get entries with same elemental composition: struc.composition.elements 
        # 4. Criteria to start logic pathway of full or partial defluorination? If ___?: 
    
    #print("formula: ", rf, "# elements", len(struc.composition.reduced_composition))
    
    p = MatProjStructure.objects.filter(    # NOTE: THIS EXCLUDES ALL FULLY DEFLUORINATED STRUCTURES
        chemical_system=struc.composition.chemical_system, #Check that partial deF entries have same chemical system as original
        is_stable=True, 
        nsites__lt=struc.composition.num_atoms, #Check that # atoms of partial deF is less than original 
        ).to_dataframe()                        # some matches could slip through here if the true partial deF structure
                                                # had more atoms in the non-reduced form
                                                # Ex. Rb2CuF6 (nsites=9) as ref would miss Rb4Cu2F8 (nsites=14)
                 
                                 
    if len(p)>0:
        i = 0
        
        for x in p:
            
            print("formula:", p.formula_full[i], "and:", p.columns[0])
            
            deltav = round(mv-p.volume_molar[i],2)  #get difference in molar volume
            #print(ff, "volume diff:", deltav)
            
            #Make the deF_table4 table with new rows
            new_row = pd.Series([p.nsites[i],
                                 p.nelements[i],p.elements[i], 
                                 p.chemical_system[i], p.density[i], 
                                 p.density_atomic[i], p.volume[i], 
                                 p.volume_molar[i], p.formula_full[i], 
                                 p.formula_reduced[i], p.formula_anonymous[i], 
                                 p.spacegroup[i], p.energy[i], 
                                 p.energy_per_atom[i],p.energy_above_hull[i],
                                 p.is_stable[i],p.decomposes_to[i],
                                 p.formation_energy[i], p.formation_energy_per_atom[i],
                                 p.id[i]], 
                                index=[p.columns[1], 
                                       p.columns[2], p.columns[3], 
                                       p.columns[4], p.columns[5], 
                                       p.columns[6], p.columns[7], 
                                       p.columns[8], p.columns[9], 
                                       p.columns[10], p.columns[11],
                                       p.columns[12], p.columns[13],
                                       p.columns[14], p.columns[15],
                                       p.columns[16], p.columns[17], 
                                       p.columns[18], p.columns[19],
                                       p.columns[20]])
            deF_table4 = pd.concat([deF_table4, new_row.to_frame().T], ignore_index=True)
            
            # Make the deltav_table3 table with new rows 
            new_deltav_row = pd.Series([ff, p.formula_full[i], rf, p.formula_reduced[i], deltav, f_id, p.id[i], f_energy, p.energy[i], 
                                        f_energyperatom, p.energy_per_atom[i],
                                        f_nsites, p.nsites[i]], 
                                       index=['Original','Formula', 'Original Reduced','Reduced Formula','Delta V', 'F id', 'De-F id', 'F energy','De-F energy','F energy per atom','De-F energy per atom','F nsites','De-F nsites'])
            deltav_table4 = pd.concat([deltav_table4, new_deltav_row.to_frame().T], ignore_index=True)
            i = i+1
            if i == len(p):
                break
    else:
        print("Zero matching for index:", mp_index)
    mp_index = mp_index +1


#%% PARTIAL 2 TABLE 4: 6/14 READ FROM DELTAV_TABLE4 AND FILTER DOWN POSSIBLE PARTIALS TO THE ONES THAT ARE CORRECT

partial_deF_table1 = pd.DataFrame()

i = 0
k=0
while i<len(deltav_table4):
    ff_orig = deltav_table4.loc[i][0]   #full fromula of original entry
    ff_new = deltav_table4.loc[i][1]    # full formula of partial deF entry
    rf_orig = deltav_table4.loc[i][2]
    rf_new = deltav_table4.loc[i][3]
    deltav = deltav_table4.loc[i][4]
    id_orig = deltav_table4.loc[i][5] 
    id_new = deltav_table4.loc[i][6] 
    energy_orig = deltav_table4.loc[i][7] 
    energy_new = deltav_table4.loc[i][8] 
    energyperatom_orig = deltav_table4.loc[i][9] 
    energyperatom_new = deltav_table4.loc[i][10] 
    nsites_orig = deltav_table4.loc[i][11] 
    nsites_new = deltav_table4.loc[i][12] 
    
    #print(deltav_table4.loc[i][3])
    
    red_form_orig = deltav_table4.loc[i][2] # will become reduced form with numbers only
    red_form_new = deltav_table4.loc[i][3]  # will become reduced form with numbers only
    
    #for x in red_form_orig:     # take reduced formula and replace all letters with spaces so that                   
     #   if x.isdigit()==False and x!="F":      # only numbers are left
      #      red_form_orig = red_form_orig.replace(x, " ")
    
    partition0_orig = red_form_orig.partition("F")[0]  #partition takes a string and splits into three parts:
    partition1_orig = red_form_orig.partition("F")[1]  # the part before quotes, part in quotes, and part after quotes
    partition2_orig = red_form_orig.partition("F")[2]  # access 0, 1, or 2 index with [] after quotes 
    #print("orig",deltav_table4.loc[i][2], ":",red_form_orig, "partition:", partition2_orig)
    
  
    #for x in red_form_new:     # take reduced formula and replace all letters with spaces so that                   
     #   if x.isdigit()==False and x!="F":      # only numbers are left
      #      red_form_new = red_form_new.replace(x, " ")
            
    partition0_new = red_form_new.partition("F")[0]
    partition1_new = red_form_new.partition("F")[1]
    partition2_new = red_form_new.partition("F")[2]
    
    #print(i,"new",deltav_table4.loc[i][3], ":",red_form_new, "partition:", partition2_new) 
    
    
    # IF STATEMENT BELOW: Compare if partition index 0 for orig and new structure are equal. Index 0 is 
    # elements before F. If they are equal, then non-F elements did not change and it will be added to list
    if partition0_orig==partition0_new and partition2_orig!=partition2_new:
        print(i,"orig",deltav_table4.loc[i][2], ":",red_form_orig, "partition:", partition2_orig)
        print(i,"new",deltav_table4.loc[i][3], ":",red_form_new, "partition:", partition2_new) 
        print(k, "Big fat win")
        print("")
        
        new_deltav_row = pd.Series([ff_orig, ff_new, rf_orig, rf_new, deltav, id_orig, id_new, energy_orig, 
                                    energy_new, energyperatom_orig, energyperatom_new, nsites_orig, nsites_new], 
                                   index=['Original','De-F Formula', 'Original Reduced','De_F Reduced','Delta V', 'F id', 'De-F id', 'F energy','De-F energy','F energy/atom','De-F energy/atom','F nsites','De-F nsites'])
        partial_deF_table1 = pd.concat([partial_deF_table1, new_deltav_row.to_frame().T], ignore_index=True)
        k=k+1    
        
    #print("")        
    
    i = i+1

#%% PARTIAL 3 TABLE 4: 6/15 CELL VOLTAGES FOR PARTIAL DEF ENTRIES

li_energy = -1.909
liF_energy = -9.6903
e = 0
partial_deF_voltage_table1 = pd.DataFrame()
#print(deltav_table.["F nsites"].values[1])
print(len(partial_deF_table1))
while e<len(partial_deF_table1):
    
    rf_orig = partial_deF_table1.loc[e][2]
    rf_new = partial_deF_table1.loc[e][3]
    
    deF_energy = partial_deF_table1.loc[e][8]
    F_energy = partial_deF_table1.loc[e][7]
    deF_energy_per_atom = partial_deF_table1.loc[e][10]
    F_energy_per_atom = partial_deF_table1.loc[e][9]
    # START HERE FOR NEW F_TRANSFER CODE
    F_reduced_form = partial_deF_table1.loc[e][2]       #reduced formula fluorinated
    deF_reduced_form = partial_deF_table1.loc[e][3]       #reduced formual defluorinated
    
    F_capital = sum(1 for c in F_reduced_form if c.isupper())

    for x in F_reduced_form:     # take reduced formula and replace all letters with spaces so that                   
        if x.isdigit()==False:      # only numbers are left
            F_reduced_form = F_reduced_form.replace(x, " ")
        

    F_reduced_form_split = F_reduced_form.split() # make workable string list of #'s from previous "for" loop
    
    F_sites_list = [] # empty list that will have integer values of subscripts
    total=0
    for c in F_reduced_form_split:   # convert string list ['2', '8'] to integer list [2, 8]
        integer = int(c)
        F_sites_list.append(integer) 
        #total = total + integer

    print(F_reduced_form_split, ":",F_sites_list)

    if F_capital!=len(F_sites_list):
        F_total_sites = F_capital-len(F_sites_list)+sum(F_sites_list)
    else:
        F_total_sites = sum(F_sites_list) 
    
    #print(formula,";", ":total", F_total_sites)
    
    #START DEF 
    deF_capital = sum(1 for c in deF_reduced_form if c.isupper())

    for x in deF_reduced_form:
        if x.isdigit()==False:
            deF_reduced_form = deF_reduced_form.replace(x, " ")

    deF_reduced_form_split = deF_reduced_form.split() #make list of #'s from formula string
     
    deF_sites_list = []
    total=0
    for c in deF_reduced_form_split:               # convert string list to integer list
        integer = int(c)
        deF_sites_list.append(integer) 
        #total = total + integer

    if deF_capital!=len(deF_sites_list):
        deF_total_sites = deF_capital-len(deF_sites_list)+sum(deF_sites_list)
    else:
        deF_total_sites = sum(deF_sites_list)
            
    #print("DEF",deF_table3.loc[e][10],":", deF_total)
    
    F_transfer = F_total_sites - deF_total_sites
    
    print(partial_deF_table1.loc[e][2],F_total_sites,"-->",partial_deF_table1.loc[e][3] , deF_total_sites,";",F_transfer)
    
    # END NEW F_TRANSFER CODE
    
    
    #deltaG = ((deF_energy+(F_transfer*liF_energy)) - (F_energy+(F_transfer*li_energy)))
    
    if F_transfer>0:    # accounts for structure pairs where initial has less F than final
                        # ( CaPdF4 --> CaPdF6 ) so F_transfer would be nefative and mess up voltage calculation
    
        deltaG = ((deF_total_sites*deF_energy_per_atom+(F_transfer*liF_energy)) - 
                  (F_total_sites*F_energy_per_atom+(F_transfer*li_energy)))
        voltage = deltaG/(-F_transfer)
    
        new_row = pd.Series([rf_orig, rf_new, voltage], index=['Original','Partial De-F', 'Voltage'])
        partial_deF_voltage_table1 = pd.concat([partial_deF_voltage_table1, new_row.to_frame().T], ignore_index=True)
    
    print(e)
    e = e+1
    
    
    # DO SAMPLE CALCULATIONS TO MAKE SURE VOLTAGES ARE CORRECT


#%% PARTIAL 4 TABLE 4: 6/15 HOW MANY POTENTIAL STRUCS DID "nsites < struc.composition.num_atoms" EXCLUDE?
    # THIS WOULD HAPPEN IF PARTIAL DEF STRUC HAS MORE NSITES THAN ORIGINAL EVEN IF REDUCED FORMULA IS CORRECT
    # Ex. Rb2CuF6 (nsites=9) as ref would miss Rb4Cu2F8 (nsites=14)


i = 0
n=0
while i<len(deF_table4):

    ff_to_rf = deF_table4.loc[i][8]
    rf = deF_table4.loc[i][9]
    #print(ff_to_rf)

    for x in ff_to_rf:
        if x==" " or x=="1":
            ff_to_rf = ff_to_rf.replace(x,"")
    
    #print(ff_to_rf, rf)
    
    if ff_to_rf==rf:
        n = n+1        
            
    i = i+1    
    #print(ff_to_rf)    

print(n)    # only 223 of 711 structures have the same full formula as reduced formula
            # this means that 488 might or might not correctly work for the nsites<struc.comp.num_atoms criteria





#%% Stuff that didn't work
# STUFF THAT DIDNT WORK
#data = MatProjStructure.objects.to_dataframe()[:15]

#data = data.reset_index()  # make sure indexes pair with number of rows
elements = []
elements2 = []
elements3 = []

# GET ELEMENTS LISTED IN EACH STRUCTURE AND PUT IN elements LIST
for index, row in data.iterrows():
    #print(row['elements'], row['chemical_system'])
    elements.append(row['chemical_system'])
    
# REMOVE "F-" FROM EACH LIST ENTRY AND PUT IN NEW elements2 LIST
# TAKES CARE OF F's AT BEGINNING AND MIDDLE OF ENTRY
for i in elements:
#if any("F-" in s for s in elements):
        #print(elements[i])
        replace1 = i.replace('F-','')
        #replace2 = i.replace('-F','')
        #print(elements[i])
        elements2.append(replace1)
        #elements2.append(replace2)
        
# REMOVE "-F" FROM elements2 LIST AND PUT IN NEW elements3 LIST
# TAKES CARE OF F's AT END OF ENTRY BUT ALSO TURNS -Fe INTO e
for i in elements2:
    replace3 = i.replace('-F','') #this will also replace -Fe with e
    elements3.append(replace3)
      

#if "F-" in elements:
 #   print(elements)
#new = elements[1].replace('F-','')


# GET VOLUME/ATOM FOR DEFLUORINATED STRUCTURES
# USE MOLAR VOLUME INSTEAD
deF_vpera = []
i = 0
for struc in deF_table[:40]:
    v = round(deF_table.volume[i]/deF_table.nsites[i], 2)
    deF_vpera.append(v)
    i = i + 1
    print("i =",i)
    if i == len(deF_table[:40]):
        break
    


#%% TEST CODE


#df = pd.DataFrame({'vals': [1, 2, 3, 4], 'ids': [u'aball', u'bball', u'cnut', u'fball']})
#df = df[df['ids'].str.contains("ball")]

#print[df]

data = {'month': ['January','February','March','April','May','June','July','August','September','October','November','December'],
        'days_in_month': [31,28,31,30,31,30,31,31,30,31,30,31]
        }

df = pd.DataFrame(data, columns = ['month', 'days_in_month'])

contain_values = df[df['month'].str.contains('Ju')]
print (contain_values)


#%% 5/6/22 FULL FORMULA F_TRANSFER NUMBER  

message = "HeLLo"

print("Capital Letters: ", sum(1 for c in message if c.isupper()))

deF_sites = deF_table3.loc[0][9]
F_sites = deltav_table3.loc[20][0]

print(F_sites)

print("Capital Letters: ", sum(1 for c in text if c.isupper()))

txt = "h3110 23 cat 444.4 rabbit 11 2 dog"
print([int(s) for s in txt.split() if s.isdigit()])

for element in deF_sites:
    if element.isdigit()==True:
        print(element)

for element in F_sites:
    if element.isdigit()==True:
        print(element)

# THIS IS WHAT WAS USED
#Take full formula and delete elements so that only spaces and numbers are left
for x in F_sites:
    if x==" ":
        g=0
    elif x.isdigit()==False:
        print("not number")
        F_sites = F_sites.replace(x, "")

F_sites = F_sites.split() #make list of #'s from formula string

print(F_sites)        
F_sites_list = []
total=0
for c in F_sites:               # convert string list to integer list
    integer = int(c)
    F_sites_list.append(integer) 
    #total = total + integer

print("line 580",F_sites_list)

f = []
for i in range(len(F_sites_list)):      #check if values in F_sites_list are even, if yes add entry to "f"
    if F_sites_list[i] % 2 == 0:
       f.append(1)
       
    
F_sites_list2 = []    
if len(f)==len(F_sites_list):   #if all vlaues = even, divide each by lowest value in list
    low = min(F_sites_list)
    high = max(F_sites_list)
    if high/low % 2 !=0:
        for i in F_sites_list:
            F_sites_list2.append(i/2)
    else:
        for i in F_sites_list:
            F_sites_list2.append(i/low)
    
print(F_sites_list2)
        
# END OF WHAT WAS USED

        
        
#%% 5/6/22 REDUCED FORMULA F_TRANSFER NUMBER  

message = "HeLLo"

print("Capital Letters: ", sum(1 for c in message if c.isupper()))

print("Capital Letters: ", sum(1 for c in text if c.isupper()))

txt = "h3110 23 cat 444.4 rabbit 11 2 dog"
print([int(s) for s in txt.split() if s.isdigit()])


# THIS IS WHAT WAS USED
#Take full formula and delete elements so that only spaces and numbers are left
deF_sites = deF_table3.loc[0][9]
F_sites = deltav_table3.loc[66][1]

F_capital = sum(1 for c in F_sites if c.isupper())
print("line 641", F_capital)

for x in F_sites:
    if x.isdigit()==False:
        F_sites = F_sites.replace(x, " ")

print("line 647",F_sites)

F_sites = F_sites.split() #make list of #'s from formula string

print("line 651",F_sites)        
F_sites_list = []
total=0
for c in F_sites:               # convert string list to integer list
    integer = int(c)
    F_sites_list.append(integer) 
    #total = total + integer

print("line 659",F_sites_list)

if F_capital!=len(F_sites_list):
    F_total = F_capital-len(F_sites_list)+sum(F_sites_list)
else:
    F_total = sum(F_sites_list)
        
print("line 666", F_total)
        
        
        
        
        
        
        
