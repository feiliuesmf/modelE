#!/bin/bash
# Regridding IPCC emissions files to the C90 cubed sphere grid
#cp /discover/nobackup/projects/giss/prod_input_files/IPCC_emissions*_1850_0.5x0.5_v1_20_04_2009.nc .

#- cleanup
rm -f paraffin_anthropogenic.nc propane+pentanes_anthropogenic.nc propane+pentanes_Biomass.nc propane+pentanes+butanes+hexanes_anthropogenic.nc
rm -f propane+pentanes+butanes+hexanes_ships.nc propane+pentanes_ships.nc propene+oalkenes_anthropogenic.nc propene+oalkenes_Biomass.nc
rm -f propene+oalkenes_ships.nc alkenes_anthropogenic.nc alkenes_Biomass.nc alkenes_ships.nc 
rm -f butanes+hexanes_* ethane+ketones_*nc propane+pentanes+butanes+hexanes_Biomass.nc paraffin_Biomass.nc

#-- CH4 Natural sources - have been adjusted and rescaled during pre-processing step in /discover/nobackup/gfaluveg/PRE/make_1x1_nat_sources/CH4
./remap.pl -par regridCH4natural.par -in CH4WETL+TUNDRA_1X1_temp -out CH4WETL+TUNDRA_C90_Dec_2009
./remap.pl -par regridCH4natural.par -in CH4TRMITE_1X1_temp -out CH4TRMITE_C90_Dec_2009
./remap.pl -par regridCH4natural.par -in CH4SOILABS_1X1_temp -out CH4SOILABS_C90_Dec_2009

#-- NOx natural sources - have been adjusted and rescaled during pre-processing step in /discover/nobackup/gfaluveg/PRE/make_1x1_nat_sources/NOx
./remap.pl -par regridNOxnatural.par -in NOx_Soil_GEIA_1x1_half -out NOx_Soil_GEIA_C90_Dec_2009

#-- CO
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CO_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CO_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CO_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CO_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_CO_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_CO_decadalmonthlymean1850_C90_Dec_2009.nc

#-- CH4
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CH4_anthropogenic_1850_0.5x0.5_v0_14_10_2009.nc -out IPCC_emissions_CH4_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_CH4_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_CH4_ships_1850_C90_Dec_2009.nc  
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_CH4_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_CH4_decadalmonthlymean1850_C90_Dec_2009.nc 

#-- NOx
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_NO_anthropogenic_1850_0.5x0.5_v1_07_05_2009.nc -out IPCC_emissions_NO_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_NO_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_NO_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_NOx_decadalmonthlymean1850_C90_Dec_2009.nc

#-- Alkenes
# first perform Alkenes = propene/42.0 + other_alkenes_and_alkynes/67.0 + ethene/28.0
ncflint -c -v emiss_shp -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_ships.nc 
ncflint -c -v emiss_shp -w 1.0,0.03571428571428571428 propene+oalkenes_ships.nc IPCC_emissions_ethene_ships_1850_0.5x0.5_v1_20_04_2009.nc -o alkenes_ships.nc 

ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02380952380952380952,0.01492537313432835820 IPCC_emissions_propene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o propene+oalkenes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,0.03571428571428571428 propene+oalkenes_anthropogenic.nc IPCC_emissions_ethene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o alkenes_anthropogenic.nc

ncflint -c -v grassfire,forestfire -w 0.02380952380952380952,0.01492537313432835820 IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean1850_v1.nc -o propene+oalkenes_Biomass.nc 
ncflint -c -v grassfire,forestfire -w 1.0,0.03571428571428571428 propene+oalkenes_Biomass.nc IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean1850_v1.nc -o alkenes_Biomass.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_anthropogenic.nc -out IPCC_emissions_alkenes_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_ships.nc -out IPCC_emissions_alkenes_ships.nc

./remap.pl -par ncregrid-ijl.par -in alkenes_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_alkenes_decadalmonthlymean1850_C90_Dec_2009.nc

#-- Paraffin
# first perform Paraffin = propane/44.0 + pentanes/72.0 + butanes/57.8 + hexanes_and_higher_alkanes/106.8 + ethane/30.0 + ketones/75.3
ncflint -c -v emiss_shp -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_ships.nc
ncflint -c -v emiss_shp -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_ships.nc
ncflint -c -v emiss_shp -w 0.03333333333333333333,0.0 IPCC_emissions_ethane_ships_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ethane_ships_1850_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_ships.nc  #- note the weight == 0.0 there is no ketones_ships file
ncflint -c -v emiss_shp -w 1.0,1.0 propane+pentanes_ships.nc butanes+hexanes_ships.nc -o propane+pentanes+butanes+hexanes_ships.nc
ncflint -c -v emiss_shp -w 1.0,1.0 propane+pentanes+butanes+hexanes_ships.nc ethane+ketones_ships.nc -o paraffin_ships.nc

ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.02272727272727272727,0.01388888888888888888 IPCC_emissions_propane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_pentanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o propane+pentanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.01730103806228373702,0.00936329588014981273 IPCC_emissions_butanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o butanes+hexanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 0.03333333333333333333,0.01328021248339973439 IPCC_emissions_ethane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc IPCC_emissions_ketones_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -o ethane+ketones_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes_anthropogenic.nc butanes+hexanes_anthropogenic.nc -o propane+pentanes+butanes+hexanes_anthropogenic.nc
ncflint -c -v emiss_dom,emiss_ind,emiss_wst -w 1.0,1.0 propane+pentanes+butanes+hexanes_anthropogenic.nc ethane+ketones_anthropogenic.nc -o paraffin_anthropogenic.nc

ncflint -c -v grassfire,forestfire -w 0.02272727272727272727,0.01388888888888888888 IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean1850_v1.nc -o propane+pentanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 0.01730103806228373702,0.00936329588014981273 IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean1850_v1.nc -o butanes+hexanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 0.03333333333333333333,0.01328021248339973439 IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean1850_v1.nc IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean1850_v1.nc -o ethane+ketones_Biomass.nc
ncflint -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes_Biomass.nc butanes+hexanes_Biomass.nc -o propane+pentanes+butanes+hexanes_Biomass.nc
ncflint -c -v grassfire,forestfire -w 1.0,1.0 propane+pentanes+butanes+hexanes_Biomass.nc ethane+ketones_Biomass.nc -o paraffin_Biomass.nc

# then do regridding
./remap.pl -par ncregrid-ijl.par -in paraffin_ships.nc -out IPCC_emissions_paraffin_shups_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in paraffin_anthropogenic.nc -out IPCC_emissions_paraffin_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in paraffin_Biomass.nc -out IPCC_GriddedBiomassBurningEmissions_paraffin_decadalmonthlymean1850_C90_Dec_2009.nc
