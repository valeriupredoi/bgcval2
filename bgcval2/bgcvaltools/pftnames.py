#
# Copyright 2014, Plymouth Marine Laboratory
#
# This file is part of the bgc-val library.
#
# bgc-val is free software: you can redistribute it and/or modify it
# under the terms of the Revised Berkeley Software Distribution (BSD) 3-clause license.

# bgc-val is distributed in the hope that it will be useful, but
# without any warranty; without even the implied warranty of merchantability
# or fitness for a particular purpose. See the revised BSD license for more details.
# You should have received a copy of the revised BSD license along with bgc-val.
# If not, see <http://opensource.org/licenses/BSD-3-Clause>.
#
# Address:
# Plymouth Marine Laboratory
# Prospect Place, The Hoe
# Plymouth, PL1 3DH, UK
#
# Email:
# ledm@pml.ac.uk
#
"""
.. module:: pftnames
   :platform: Unix
   :synopsis: A list of names used for makeing text on plots pretty.
.. moduleauthor:: Lee de Mora <ledm@pml.ac.uk>

"""

from calendar import month_name
from ..UKESMpython import AutoVivification, AutoVivToYaml, folder, YamlToDict
#from itertools import product
from os.path import exists
import numpy as np

#####
#

regions = [
    'Surface',
    '200m',
    '100m',
    '500m',
    '1000m',
    'Transect',
    'All',
    '',
]

MaredatTypes = ['chl', 'diatoms', 'bac', 'mesozoo', 'picophyto', 'microzoo']

Ocean_names = [
    'SouthPacificOcean',
    'ArcticOcean',
    'AntarcticOcean',
    'NorthAtlanticOcean',
    'SouthAtlanticOcean',
    'NorthPacificOcean',
    'IndianOcean',
    'EquatorialPacificOcean',
    'EquatorialAtlanticOcean',
]

IFREMERTypes = [
    'mld',
    'mld_DT02',
    'mld_DR003',
    'mld_DReqDTm02',
]

WOATypes = [
    'silicate', 'nitrate', 'phosphate', 'salinity', 'temperature', 'oxygen'
]

CMIP5models = [
    'MEDUSA',
    'ERSEM',
    'BNU-ESM',
    'IPSL-CM5A-LR',
    'CESM1-BGC',
    'IPSL-CM5A-MR',
    'CMCC-CESM',
    'IPSL-CM5B-LR',
    'CNRM-CM5',
    'MPI-ESM-LR',
    'GFDL-ESM2G',
    'MPI-ESM-MR',
    'GFDL-ESM2M',
    'MRI-ESM1',
    'HadGEM2-CC',
    'NorESM1-ME',
    'HadGEM2-ES',
]

TAKAHASHITypes = [
    'pCO2',
]

GEOTRACESTypes = [
    'iron',
]

BGCmodels = [
    'Diat-HadOCC',
    'ERSEM',
    'HadOCC',
    'MEDUSA',
    'PlankTOM6',
    'PlankTOM10',
]

Seasons = ['JFM', 'AMJ', 'JAS', 'OND']

Hemispheres = [
    'NorthHemisphere',
    'SouthHemisphere',
]

months = [m for m in month_name
          if m]  # Because months starts at 1, and 0 is empty.
OceanMonth_names = [o + m for o in Ocean_names for m in months]
OceanSeason_names = [o + s for o in Ocean_names for s in Seasons]
HemispheresMonths = [h + m for h in Hemispheres for m in months]
SouthHemispheresMonths = [
    h + m for h in [
        'SouthHemisphere',
    ] for m in months
]
NorthHemispheresMonths = [
    h + m for h in [
        'NorthHemisphere',
    ] for m in months
]


def makeLongNameDict():
    lnd = {}
    SameSame = [
        'temperature',
        "salinity",
        "nitrate",
        "phosphate",
        'silicate',
        'oxygen',
        'iron',
        'chlorophyll',
        'chl',
        'alkalinity',
        'surface',
        'oceans',
        'ocean',
        'months',
        'month',
        'depth',
        'depths',
        'biomass',
        'Tropics',
        'Temperate',
        'Arctic',
        'Depth',
        'Overestimate',
        'Underestimate',
        'Matched',
        'Seasons',
    ]

    for txt in SameSame:
        title = txt.title()
        lnd[txt] = title

    #####
    #Specific

    #####
    # Data names:
    lnd['intpp'] = 'Integrated Primary Production'
    lnd['ppint'] = 'Integrated Primary Production'
    lnd['PP'] = 'Primary Production'
    lnd['netPP'] = 'Net Primary Production'
    lnd['IntegratedPrimaryProduction'] = "Integrated Primary Production"
    lnd['IntegratedPrimaryProduction_OSU'] = "Integrated Primary Production (OSU)"
    lnd['TotalIntegratedPrimaryProduction'] = "Total Integrated Primary Production"

    lnd['ZMI'] = 'Microzooplankton'
    lnd['ZME'] = 'Mesozooplankton'
    lnd['PHD'] = 'Diatoms'
    lnd['bac'] = 'Bacteria'

    lnd['CHL'] = 'Chlorophyll'
    lnd['Chlorophyll'] = "Chlorophyll"
    lnd['Chlorophylla'] = 'Chlorophyll'
    lnd['Chlorophyll_cci'] = "Chlorophyll (CCI)"
    lnd['Chlorophyll_pig'] = "Chlorophyll (Pigments)"
    lnd['TotalChlorophyll'] = 'Total Chlorophyll' 
    lnd['CHD'] = "Diatom Chlorophyll"
    lnd['CHN'] = "Non-Diatom Chlorophyll"
    lnd['DiatomChlorophyll'] = "Diatom Chlorophyll"
    lnd['NonDiatomChlorophyll'] = "Non-Diatom Chlorophyll"
    lnd['DiaFrac'] = "Diatom Fraction"
    lnd['Dust'] = "Dust"

    lnd['picophyto'] = 'Picophytoplankton'
    lnd['microzoo'] = 'Microzooplankton'
    lnd['mesozoo'] = 'Mesozooplankton'
    lnd['diatoms'] = 'Diatoms'

    lnd['Seawifs-micro'] = 'Seawifs Microphyto. chl.'
    lnd['Seawifs-nano'] = 'Seawifs Nanophyto. chl.'
    lnd['Seawifs-pico'] = 'Seawifs Picophyto. chl.'
    lnd['SeawifsBM-micro'] = 'Seawifs Microphyto. Biomass'
    lnd['SeawifsBM-nano'] = 'Seawifs Nanophyto. Biomass'
    lnd['SeawifsBM-pico'] = 'Seawifs Picophyto. Biomass'
    lnd['Seawifs-biomass'] = 'Phytoplankton Biomass'

    lnd['AMOC'] = "AMOC"
    lnd['AMOC_26N'] = "AMOC 26N"
    lnd['AMOC_26N_nomexico'] = "AMOC 26N (excluding Gulf of Mexico)"
    lnd['AMOC_32S'] = "AMOC 32S"
    lnd['ADRC_26N'] = "Atlantic Deep Return Current at 26N"

    lnd['ZonalCurrent'] = "Zonal Current"
    lnd['MeridionalCurrent'] = "Meridional Current"
    lnd['VerticalCurrent'] = "Vertical Current"
    lnd['WindStress'] = "Wind Stress"

    lnd['GlobalMeanTemperature'] = "Global Volume-weighted Mean Temperature"
    lnd['VolumeMeanTemperature'] = "Volume-weighted Mean Temperature"
    lnd['GlobalMeanSalinity'] = "Global Volume-weighted Mean Salinity"

    lnd['VolumeMeanTemperature'] = "Volume-weighted Mean Temperature"
    lnd['GlobalMeanSalinity'] = "Global Volume-weighted Mean Salinity"
   

    lnd['TotalHeatFlux'] = "Global Total Heat Flux"
    lnd['HeatFlux'] = "Heat Flux"

    lnd['scvoltot'] = "Sea Water Volume"
    lnd['soga'] = "Global Average Sea Water Salinity"
    lnd['thetaoga'] = "Global Average Sea Water Potential Temperature"
    lnd['scalarHeatContent'] = "Global heat content"

    lnd['sowaflup'] = "Net Upward Water Flux"
    lnd['sohefldo'] = "Net Downward Heat Flux"
    lnd['sofmflup'] = "Water flux due to freezing/melting"
    lnd['sosfldow'] = "Downward salt flux"
    lnd['soicecov'] = "Ice Fraction"
    lnd['sossheig'] = "Sea Surface Height"
    lnd['FreshwaterFlux'] = "Freshwater Flux"

    lnd['exportRatio'] = "Export Ratio"
    lnd['LocalExportRatio'] = "Export Ratio"

    lnd['PCO2_SW'] = 'pCO2'
    lnd['pCO2'] = 'pCO2'

    lnd['iron'] = "Iron"
    lnd['Fe_D_CONC_BOTTLE'] = "Iron (Dissolved)"

    lnd['AirSeaFluxCO2'] = "Air Sea CO2 Flux"
    lnd['TotalAirSeaFluxCO2'] = "Total Air Sea CO2 Flux"
    lnd['NoCaspianAirSeaFluxCO2'] = "Air Sea CO2 Flux"

    lnd['DIC'] = "DIC"
    lnd['pH'] = "pH"

    lnd['TotalDust'] = 'Total Dust'

    lnd['OMZ'] = 'Oxygen minimum zone'
    lnd['OMZExtent'] = 'Oxygen minimum zone Extent'
    lnd['ExtentMaps'] = 'Extent Maps'
    lnd['TotalOMZVolume'] = 'Total Oxygen minimum zone volume (<20 mmol O2/m^3)'
    lnd['TotalOMZVolume50'] = 'Total Oxygen minimum zone volume (<50 mmol O2/m^3)'
    lnd['OMZThickness'] = 'Oxygen minimum zone thickness (<20 mmol O2/m^3)'
    lnd['OMZMeanDepth'] = 'Oxygen minimum zone mean depth (<20 mmol O2/m^3)'
    lnd['AOU'] = 'Apparent Oxygen Usage'
    lnd['VolumeMeanOxygen'] = "Volume-weighted Mean Oxygen"

    lnd['DrakePassageTransport'] = 'Drake Passage Transport'
    lnd['DPT'] = 'Drake Passage Transport'

    lnd['t_mn'] = 'Mean Temperature'
    lnd['t_an'] = 'Temperature'
    lnd['s_mn'] = 'Mean Salinity'
    lnd['s_an'] = 'Salinity'
    lnd['n_mn'] = 'Mean Nitrate'
    lnd['n_an'] = 'Nitrate'
    lnd['p_mn'] = 'Mean Phosphate'
    lnd['p_an'] = 'Phosphate'
    lnd['p_an'] = 'Oxygen'
    lnd['i_mn'] = 'Mean Silicate'
    lnd['i_an'] = 'Silicate'

    lnd['NorthernTotalIceExtent'] = 'Northern Hemisphere Ice Extent'
    lnd['SouthernTotalIceExtent'] = 'Southern Hemisphere Ice Extent'
    lnd['WeddelIceExent'] = 'Weddel Sea Ice Extent'
    lnd['TotalIceExtent'] = 'Total Ice Extent'

    lnd['AMM'] = 'Atlantic Margins'
    lnd['AMM_Shelf'] = 'Atlantic Margins Shelf'
    lnd['AMM_OffShelf'] = 'Atlantic Margins Off Shelf'
    lnd['AMM7_O3c'] = 'DIC'
    lnd['AMM7_talk'] = 'TA'
    lnd['AMM7_SST'] = 'Temperature'
    lnd['AMM7_SSS'] = 'Salinity'

    lnd['MA_SST'] = 'Sea Surface Temperature'
    lnd['MA_SSS'] = 'Sea Surface Salinty'
    lnd['SSS'] = 'Sea Surface Salinty'

    lnd['MA_DraKE'] = 'Drake Passage Current'
    lnd['MA_AMOC_26N'] = 'AMOC at 26.5N' 
    lnd['MA_ZOS'] = 'Sea Surface Height'
    lnd['MA_MLD'] = 'Mixed Layer Depth'
    lnd['MA_MLD_Sigma'] = 'Mixed Layer Depth (Sigma)'
    lnd['MA_Nitrate'] = 'Nitrate'
    lnd['MA_Phosphate'] = 'Phosphate'
    lnd['MA_Iron'] = 'Iron'
    lnd['MA_Silicate'] = 'Silicate'
    lnd['MA_TotalChlorophyll'] = 'Total Chlorophyll'
    lnd['MA_TotalPhytoC'] = 'Total Phytoplankton Carbon'
    lnd['MA_TotalZooC'] = 'Total Zooplankton Carbon'
    lnd['MA_NPP'] = 'Net Primary Production'
    lnd['MA_GPP'] = 'Gross Primary Production'
    lnd['MA_pH'] = 'pH'
    lnd['MA_O2'] = 'Oxygen'
    lnd['MA_Oxygen'] = 'Oxygen'
    lnd['MA_DIC'] = 'Dissolved Inorganic Carbon'
    
    lnd['MA_SouthernTotalIceExtent'] = 'Southern Total Ice Extent'
    lnd['MA_NorthernTotalIceExtent'] = 'Northern Total Ice Extent'
    lnd['MA_TotalIceExtent'] = 'Total Ice Extent'

    lnd['SouthernTotalIceArea'] = 'Southern Hemisphere Ice Area'
    lnd['TotalIceArea'] = 'Total Ice Area'

    lnd['MIZ'] = 'Marginal Ice Zone'
    lnd['NorthernMIZArea'] = 'Northern Hemisphere Marginal Ice Zone Area'
    lnd['SouthernMIZArea'] = 'Southern Hemisphere Marginal Ice Zone Area'
    lnd['TotalMIZArea'] = 'Total Marginal Ice Zone Area'
    lnd['NorthernMIZfraction'] = 'Northern Hemisphere Marginal Ice Zone fraction'
    lnd['SouthernMIZfraction'] = 'Southern Hemisphere Marginal Ice Zone fraction'
    lnd['TotalMIZfraction'] = 'Total Marginal Ice Zone fraction'

    lnd['JFM'] = 'JFM'
    lnd['AMJ'] = 'AMJ'
    lnd['JAS'] = 'JAS'
    lnd['OND'] = 'OND'

    lnd['mld'] = 'Mixed Layer Depth'
    lnd['MaxMonthlyMLD'] = 'Annual Maximum Mixed Layer Depth'
    lnd['MinMonthlyMLD'] = 'Annual Minimum Mixed Layer Depth'

    lnd['mld_DT02'] = 'MLD:Fixed Threshold Temperature '
    lnd['mld_DR003'] = 'MLD:Fixed Threshold Density'
    lnd['mld_DReqDTm02'] = 'MLD:Variable Threshold Density'

    #####
    # Depth layers/transects.
    lnd['AtlanticTransect'] = "Atlantic Transect"
    lnd['Atlantic28W'] = "Atlantic (28 W)"
    lnd['ArcTransect'] = "Arctic Transect"
    lnd['AntTransect'] = "Antarctic Transect"
    lnd['CanRusTransect'] = "Canada-Siberia Arctic Transect"
    lnd['PacificTransect'] = "Pacific Transect"
    lnd['SouthernTransect'] = "Southern Transect"
    lnd['SOTransect'] = "Southern Ocean Transect"
    lnd['Transect'] = "Atlantic Transect"
    lnd['PTransect'] = "Pacific Transect"
    lnd['Pacific135W'] = "Pacific (135 W)"

    lnd['NordicSea'] = "Nordic Sea"
    lnd['LabradorSea'] = "Labrador Sea"
    lnd['NorwegianSea'] = "Norwegian Sea"
    lnd['Cornwall'] = "Cornwall"

    lnd['100m'] = "100m deep"
    lnd['200m'] = "200m deep"
    lnd['500m'] = "500m deep"
    lnd['1000m'] = "1000m deep"
    lnd['10N'] = "10 degree North Transect"
    lnd['10S'] = "10 degree South Transect"

    #####
    # Names from plots.
    lnd['hist'] = "Histogram"
    lnd['scatter'] = "scatter diagram"
    lnd['percentiles'] = "Time series"
    lnd['mean'] = "mean"
    lnd['median'] = "median"
    lnd['robinquad'] = "Maps"
    lnd['robinquad-cartopy'] = "Interpolated Map"
    lnd['hov'] = "Hovmoeller"
    lnd['10-90pc'] = ""
    lnd['Target'] = "Target Diagram"
    lnd['Taylor'] = "Taylor Diagram"
    lnd['SummaryTargets'] = "Summary Diagrams"
    lnd['RobustTarget'] = "Robust Statisitcs Target Diagram"
    lnd['regionless'] = ''
    lnd['layerless'] = ''
    lnd['metricless'] = ''
    lnd['together'] = ''
    lnd['dataonly'] = ''
    lnd['wcvweighted'] = ''

    lnd['so'] = 'NEMO Diagnostics'

    lnd['RegionLegend'] = "Region Legend"
    lnd['TransectsLegend'] = "Transects Legend"
    lnd['TransectsLegendBoth'] = "Transects Legend"

    #####
    # DMS fields:
    lnd['anderson'] = 'Anderson et al.'
    lnd['dms_and'] = 'Anderson et al.'
    lnd['dms_andSurface'] = 'Anderson et al.'
    lnd['dms_p_and'] = 'Anderson et al.'
    lnd['dms_p_andSurface'] = 'Anderson et al.'
    lnd['dms_p_and1'] = 'DMS (Anderson - all CHL)'
    lnd['dms_p_and1Surface'] = 'DMS (Anderson - all CHL)'
    lnd['dms_p_and2'] = 'DMS (Anderson - CHN only)'
    lnd['dms_p_and2Surface'] = 'DMS (Anderson - CHN only)'

    lnd['aranamit'] = 'Aranami et al.'
    lnd['dms_ara'] = 'Aranami et al.'
    lnd['dms_araSurface'] = 'Aranami et al.'
    lnd['dms_p_ara'] = 'Aranami et al.'
    lnd['dms_p_araSurface'] = 'Aranami et al.'
    lnd['dms_p_ara1'] = 'DMS (Aranamit - all CHL)'
    lnd['dms_p_ara1Surface'] = 'DMS (Aranamit - all CHL)'
    lnd['dms_p_ara2'] = 'DMS (Aranamit - CHN only)'
    lnd['dms_p_ara2Surface'] = 'DMS (Aranamit - CHN only)'

    lnd['halloran'] = 'Halloran et al.'
    lnd['dms_hal'] = 'Halloran et al.'
    lnd['dms_halSurface'] = 'Halloran et al.'
    lnd['dms_p_hal'] = 'Halloran et al.'
    lnd['dms_p_halSurface'] = 'Halloran et al.'
    lnd['dms_p_hal1'] = 'DMS (Halloran - all CHL)'
    lnd['dms_p_hal1Surface'] = 'DMS (Halloran - all CHL)'
    lnd['dms_p_hal2'] = 'DMS (Halloran - CHN only)'
    lnd['dms_p_hal2Surface'] = 'DMS (Halloran - CHN only)'

    lnd['simodach'] = 'Simo & Dach'
    lnd['dms_sim'] = 'Simo & Dach'
    lnd['dms_simSurface'] = 'Simo & Dach'
    lnd['dms_p_sim'] = 'Simo & Dach'
    lnd['dms_p_simSurface'] = 'Simo & Dach'
    lnd['dms_p_sim1'] = 'DMS (Simodach - all CHL)'
    lnd['dms_p_sim1Surface'] = 'DMS (Simodach - all CHL)'
    lnd['dms_p_sim2'] = 'DMS (Simodach - CHN only)'
    lnd['dms_p_sim2Surface'] = 'DMS (Simodach - CHN only)'

    lnd['LANA'] = 'Lana et al. (extrapolated)'
    lnd['LANA_p'] = 'Lana et al. (pixels)'
    lnd['lanaetal'] = 'DMS extrapolated (Lana et al. 2011)'
    lnd['DMS'] = 'DMS'
    lnd['DMS_p'] = 'DMS (pixels)'
    lnd['DMS_e'] = 'DMS (extrapolated)'

    #####
    # Specific regions and slices
    lnd['Top40m'] = "Top 40m"
    lnd['Top200m'] = "Top 200m"
    lnd['Top40mNoArtics'] = "Top 40m (No Arctics)"
    lnd['Top200mNoArtics'] = "Top 200m (No Arctics)"

    lnd['NoShelf'] = "No Shelf"
    lnd['NoShelfTop40'] = "No Shelf (Top 40m)"
    lnd['NoShelfSurface'] = "No Shelf (Surface)"

    lnd['ArcticOcean'] = "Arctic Ocean"
    lnd['AntarcticOcean'] = "Antarctic Ocean"
    lnd['NorthAtlanticOcean'] = "North Atlantic Ocean"
    lnd['SouthAtlanticOcean'] = "South Atlantic Ocean"
    lnd['NorthPacificOcean'] = "North Pacific Ocean"
    lnd['SouthPacificOcean'] = "South Pacific Ocean"
    lnd['EquatorialAtlanticOcean'] = "Equatorial Atlantic Ocean"
    lnd['EquatorialPacificOcean'] = "Equatorial Pacific Ocean"
    lnd['IndianOcean'] = "Indian Ocean"
    lnd['NorthHemisphere'] = "North Hemisphere"
    lnd['SouthHemisphere'] = "South Hemisphere"
    lnd['26N'] = "26N"
    lnd['32S'] = "32S"

    lnd['WeddelSea'] = "Weddel Sea"
    lnd['Enderby'] = "Enderby Region"  # Regions from Pierce 1995 - https://doi.org/10.1175/1520-0485(1995)025<2046:CROHAF>2.0.CO;2
    lnd['Wilkes'] = "Wilkes Region"
    lnd['Ross'] = "Ross Region"
    lnd['Amundsen'] = "Amundsen Region"
    lnd['Weddel'] = "Weddel Region"

    lnd['Global'] = "Global"
    lnd['Equator10'] = "Equator (+/-10)"
    lnd['Remainder'] = "Oligotrophic Gyres"
    lnd['ArcticOcean'] = "Arctic Ocean"
    lnd['NorthernSubpolarAtlantic'] = "Northern Subpolar Atlantic"
    lnd['NorthernSubpolarPacific'] = "Northern Subpolar Pacific"

    lnd['SouthernOcean'] = "Southern Ocean"
    lnd['AtlanticSOcean'] = 'Atlantic S. Ocean'
    lnd['NorthTemperate'] = "North Temperate"
    lnd['SouthTemperate'] = "South Temperate"
    lnd['NorthTropics'] = "North Tropics"
    lnd['SouthTropics'] = "South Tropics"
    lnd['NorthArctic'] = "North Arctic"
    lnd['Antarctic'] = "Antarctic"
    lnd['Equatorial'] = "Equatorial"
    lnd['AMT'] = "AMT"
    lnd['AMTTop40m'] = "AMT (Top 40m)"
    lnd['AMTTop200m'] = "AMT (Top 200m)"

    lnd['BlackSea'] = "Black Sea"
    lnd['RedSea'] = "Red Sea"
    lnd['BalticSea'] = "Baltic Sea"
    lnd['PersianGulf'] = "Persian Gulf"

    lnd['ignoreInlandSeas'] = "No Inland Seas"
    lnd['ignoreBlackSea'] = "No Black Sea"
    lnd['ignoreRedSea'] = "No Red Sea"
    lnd['ignoreBalticSea'] = "No Baltic Sea"
    lnd['ignorePersianGulf'] = "No Persian Gulf"
    lnd['ignoreInlandSeas'] = "No Inland Seas"
    lnd['ignoreMediteranean'] = "No Mediteranean"
    lnd['ignoreExtraArtics'] = "No Arctic Oceans (50 degrees)"
    lnd['ignoreMoreArtics'] = "No Arctic Oceans (60 degrees)"
    lnd['ignoreMidArtics'] = "No Arctic Oceans (65 degrees)"
    lnd['ignoreArtics'] = "No Arctic Oceans (70 degrees)"
    lnd['ignoreCaspian'] = "Ignoring Caspian Sea"

    lnd['CHL_JJA'] = '(with CCI JJA mask)'
    lnd['CHL_DJF'] = '(with CCI DJF mask)'

    lnd['Depth_0-10m'] = 'Depth <10m'
    lnd['Depth_10-20m'] = '10m <= Depth < 20m'
    lnd['Depth_20-50m'] = '20m <= Depth < 50m'
    lnd['Depth_50-100m'] = '50m <= Depth < 100m'
    lnd['Depth_100-500m'] = '100m <= Depth < 500m'
    lnd['Depth_500m'] = 'Depth > 500m'
    lnd['Depth_700m'] = 'Depth < 700m'
    lnd['Depth_1000m'] = 'Depth > 1000m'

    lnd['Depth_0-50m'] = 'Depth <50m'
    lnd['Depth_50-100m'] = '50m <= Depth < 100m'
    lnd['Depth_100-200m'] = '100m <= Depth < 200m'
    lnd['Depth_200-500m'] = '200m <= Depth < 500m'
    lnd['Depth_500-1000m'] = '500m <= Depth < 1000m'
    lnd['Depth_700-2000m'] = '700m <= Depth < 2000m'
    lnd['Depth_1000-2000m'] = '1000m <= Depth < 2000m'
    lnd['Depth_2000m'] = 'Depth < 2000m'

    lnd['0-1pc'] = '0-1 percentiles'
    lnd['1-5pc'] = '1-5 percentiles'
    lnd['5-25pc'] = '5-25 percentiles'
    lnd['25-40pc'] = '25-40 percentiles'
    lnd['40-60pc'] = '40-60 percentiles'
    lnd['60-75pc'] = '60-75 percentiles'
    lnd['75-95pc'] = '75-95 percentiles'
    lnd['95-99pc'] = '95-99 percentiles'
    lnd['99-100pc'] = '99-100 percentiles'
    lnd['nonZero'] = 'Non zero'
    lnd['aboveZero'] = ''
    lnd['1-99pc'] = '1-99 percentiles'
    lnd['5-95pc'] = '5-95 percentiles'
    lnd['0-99pc'] = 'up to 99th percentile'

    lnd['All'] = 'Global'
    lnd['Best'] = 'Best'
    lnd['Standard'] = 'Standard'
    lnd['SalArtifact'] = 'Salinity Artifact (<15psu)'
    lnd['NitArtifact'] = 'Nitrogen Artifact'
    lnd['TypicalIron'] = 'Iron < 4 umol m^-3'
    lnd['HighLatWinter'] = 'High Latitude Winter'
    lnd['OffAxis'] = 'Off Axis'
    lnd['Depth'] = 'Depth >200m'
    lnd['Shallow'] = 'Depth <200m'

    lnd['Sum'] = 'Total'

    lnd['Overestimate_2sig'] = "Overestimate 2 sigma"
    lnd['Overestimate_3sig'] = "Overestimate 3 sigma"
    lnd['Underestimate_2sig'] = "Overestimate 2 sigma"
    lnd['Underestimate_3sig'] = "Overestimate 3 sigma"

    lnd['OceansMonths'] = 'Oceans Months'
    lnd['maskBelowBathy'] = 'Masked Below Bathymetery'
    lnd['OnShelf'] = 'On Shelf'
    lnd['OffShelf'] = 'Off Shelf'

    lnd['movingav30years'] = '30 year moving average'
    lnd['movingav100years'] = '100 year moving average'

    #####
    # Models
    lnd['NEMO'] = 'NEMO'
    lnd['CICE'] = 'CICE'
    lnd['ERSEM'] = 'ERSEM'
    lnd['MAREDAT'] = 'MAREDAT'
    lnd['MEDUSA'] = 'MEDUSA'

    #####
    # Data sets
    lnd['Takahashi'] = 'Takahashi 2009'
    lnd['Takahashi2009'] = 'Takahashi 2009'
    lnd['Seawifs'] = 'Seawifs'
    lnd['cci'] = "CCI"
    lnd['pig'] = "pigments"
    lnd['OSU'] = "(OSU)"
    lnd['1x1'] = "(iMarNet)"
    lnd['GEOTRACES'] = "GEOTRACES"
    lnd['WOA'] = 'WOA'
    lnd['IFREMER'] = "IFREMER"
    lnd['InitialConditions'] = 'Initial Conditions'
    lnd['GLODAP'] = 'GLODAP'
    lnd['GLODAPv2'] = 'GLODAP v2'

    for txt in list(lnd.keys()):
        longname = lnd[txt]
        lnd[longname] = longname
        lnd[txt.lower()] = longname
        lnd[txt.upper()] = longname
        lnd[txt.title()] = longname

    ThreeDFields = [
        'DIC',
        'temperature',
        "salinity",
        "nitrate",
        "phosphate",
        'silicate',
        'oxygen',
        'iron',
        'chl',
    ]
    depthLevels = [
        'Surface',
        '100m',
        '200m',
        '500m',
        '1000m',
        'Transect',
        'PTransect',
        'SOTransect',
        'ArcTransect',
        'AntTransect',
        'CanRusTransect',
        '10N',
        '10S',
    ]

    for f in ThreeDFields:
        for d in depthLevels:
            if d in [
                    '100m',
                    '200m',
                    '500m',
                    '1000m',
            ]:
                longname = lnd[f] + ' (' + lnd[d] + ')'
            else:
                longname = lnd[d] + ' ' + lnd[f]

            lnd[f + d] = longname
            lnd[d + f] = longname
            lnd[(d + f).upper()] = longname
            lnd[(f + d).upper()] = longname
            lnd[(d + f).lower()] = longname
            lnd[(f + d).lower()] = longname
            lnd[f.title() + d.title()] = longname
            lnd[d.title() + f.title()] = longname

    return lnd


longNameDict = makeLongNameDict()


def getLongName(text, debug=False):
    if debug: print("Getting long name:", text)


    if type(text) in [type([
            'a',
            'b',
    ]), type((
            'a',
            'b',
    ))]:
        return ' '.join([getLongName(t) for t in text])
        #out = ''
    if text in [None, '', ' ',]: return ''

    try:
        return longNameDict[text]
    except:
        if debug: print("text not in dict:", text)
    try:
        return longNameDict[text.lower()]
    except:
        if debug: print("text.lower() not in dict:", text.lower())
    print("text not in dict:", text)

    return text


def fancyUnits(units, debug=False):
    """	
	Converts ascii units string into latex style formatting.
	"""
    units = units.replace('[', '').replace(']', '')

    #if units in ['mg C/m^3','mg C/m^2',]:		return 'mg C m'+r'$^{-3}$'
    if units in [
            'umol/l, uM, mo/l, ug/l, ',
    ]:
        return 'mg m' + r'$^{-3}$'  # silly nitrates multi units
    if units in [
            'mg C/m^3',
    ]: return 'mg C m' + r'$^{-3}$'
    if units in [
            'mg Chl/m3',
            'ng/L',
            'mgCh/m3',
    ]:
        return 'mg Chl m' + r'$^{-3}$'
    if units in [
            'mg C/m^3/d',
    ]: return 'mg C m' + r'$^{-3}$/day'
    if units in [
            'mg N/m^3',
    ]: return 'mg N m' + r'$^{-3}$'
    if units in [
            'mg P/m^3',
    ]: return 'mg P m' + r'$^{-3}$'
    if units in ['mmol N/m^3', 'mmol-N/m3']: return 'mmol N m' + r'$^{-3}$'
    if units in [
            'mmol P/m^3',
    ]: return 'mmol P m' + r'$^{-3}$'
    if units in [
            'mmol C/m^3',
    ]: return 'mmol C m' + r'$^{-3}$'
    if units in [
            'umol F/m^3',
    ]: return r'$\mu$' + 'mol m' + r'$^{-3}$'
    if units in [
            'umol /m^3',
            'umol / m3',
    ]:
        return r'$\mu$' + 'mol m' + r'$^{-3}$'
    if units in [
            'mmol S/m^3',
    ]: return 'mmol S m' + r'$^{-3}$'
    if units in [
            'mmolSi/m3',
            'mmol Si/m^3',
    ]: return 'mmol Si m' + r'$^{-3}$'
    if units in [
            'mmolFe/m3',
    ]: return 'mmol Fe m' + r'$^{-3}$'

    if units in [
            'mg/m2/s',
    ]: return 'mg m' + r'$^{-2}$' + ' s' + r'$^{-1}$'
    if units in ['J', 'Joules']: return 'Joules'

    if units in [
            'ug/l',
            'mg/m^3',
            'ug/L',
    ]: return 'mg m' + r'$^{-3}$'
    if units in [
            '10^12 g Carbon year^-1',
    ]:
        return r'$10^{12}$' + ' g Carbon/year'
    if units in [
            'mol C/m^',
    ]: return 'mol C/m' + r'$^{2}$'
    if units in [
            'mmmol/m^3', 'mmol/m^3', 'umol/l', 'micromoles/l', 'mmolO2/m3'
    ]:
        return 'mmol m' + r'$^{-3}$'
    if units in ['mmol/m^2']: return 'mmol m' + r'$^{-2}$'
    #if units in ['mmol/m^3']:			return 'mmol m'+r'$^{-3}$'
    if units in [
            'degrees Celsius',
            'degreesC',
            'C',
            'degC',
            'degrees_celsius',
    ]:
        return r'$\,^{\circ}\mathrm{C}$'
    if units in [
            'psu',
            'PSU',
    ]: return 'psu'
    #if units in ['umol/l',]:			return r'$\mu$'+'mol/l'
    if units in [
            'm',
            'meters',
            'meter',
    ]: return 'm'
    if units in [
            '1/m',
    ]: return r'$\mathrm{m}^{-1}$'
    if units in [
            'm/s',
    ]: return r'$\mathrm{ms}^{-1}$'
    #if units in ['ug/l']:			#	return 'mg m'+r'$^{-3}$'
    if units in ['W/m^2']: return 'W m' + r'$^{-2}$'
    if units in [
            'umol/kg',
    ]: return r'$\mu$' + 'mol kg' + r'$^{-1}$'
    if units in [
            'nmol/kg',
    ]: return 'nmol kg' + r'$^{-1}$'
    if units in [
            'tons C/d',
    ]: return 'tons C/day'
    if units in ['ug/L/d', 'ug                  ']:
        return 'mg m' + r'$^{-3}$' + '/day'  #yes, there are lots of spaces
    if units.replace(' ', '') in [
            'ug',
    ]: return r'$\mu$' + 'g'  #r'$\mu$'+
    if units in [
            '1',
    ]:
        print('fancyUnits:\tWarning:\tStrange units:', units)
        return ''
    if units in [
            'uatm',
    ]: return r'$\mu$' + 'atm'
    if units in [
            'ppmv',
    ]: return 'ppm'
    if units in [
            'milliliters_per_liter',
    ]: return 'ml/l'
    print('fancyUnits:\tERROR:\t', units, ' not found in fancyUnits.')
    if debug:
        assert False
    return units
