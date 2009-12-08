#!/bin/bash
# Regridding IPCC emissions files to the C90 cubed sphere grid

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
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_propene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_propene_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_other_alkenes_and_alkynes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_propene_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_ethene_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_ethene_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_propene_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_propene_ships_1850_C90_Dec_2009.nc 
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_other_alkenes_and_alkynes_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_other_alkenes_and_alkynes_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_ethene_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_ethene_ships_1850_C90_Dec_2009.nc 

./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_propene_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_other_alkenes_and_alkynes_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean1850_v1.nc -out  IPCC_GriddedBiomassBurningEmissions_ethene_decadalmonthlymean1850_C90_Dec_2009.nc

#-- Paraffin
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_propane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_propane_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_pentanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_pentanes_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_butanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_butanes_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_hexanes_and_higher_alkanes_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_ethane_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_ethane_anthropogenic_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_ketones_anthropogenic_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_ketones_anthropogenic_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_propane_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_propane_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_pentanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_pentanes_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_butanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_butanes_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_hexanes_and_higher_alkanes_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_hexanes_and_higher_alkanes_ships_1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_emissions_ethane_ships_1850_0.5x0.5_v1_20_04_2009.nc -out IPCC_emissions_ethane_ships_1850_C90_Dec_2009.nc

./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_propane_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_pentanes_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_butanes_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_hexanes_and_higher_alkanes_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_ethane_decadalmonthlymean1850_C90_Dec_2009.nc
./remap.pl -par ncregrid-ijl.par -in IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean1850_v1.nc -out IPCC_GriddedBiomassBurningEmissions_ketones_decadalmonthlymean1850_C90_Dec_2009.nc
