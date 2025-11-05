# RiOMar

A collection of code used to process and analyse data relating to the [RiOMar project](https://riomar.lsce.ipsl.fr/).

Note that the majority of the code in the __/func__ folder was originally written by [Louis Terrats](https://github.com/louis-terrats) under an MIT license and may be found on [GitHub](https://github.com/louis-terrats/myRIOMAR_dev) where it may be installed and used as a standalone module.

The code is being repurposed in this repository as a continuation of the work previously performed, and with an eye towards the next two years of development. Terminating with the completion of the RiOMar project itself.

The primary objective of this repository is to provide a clear layout of the workflows used to verify the integrity of a number of data products, as well as any analyses used for publications.

The workflow may be found in the __/code__ folder, with the order of operatioons given in the names of the scripts.

## Data sources

### In situ data

#### SOMLIT

#### REPHY

### Wind

The dedicated CMEMS wind reanalysis product was used in this study: https://data.marine.copernicus.eu/product/WIND_GLO_PHY_L4_NRT_012_004/description

### Tides

Tide gauge data were downloaded via web form at: https://data.shom.fr/donnees/refmar/download

The codes used were:
- LE_HAVRE (Seine)
- SAINT-NAZAIRE (Loire)
- PORT-BLOC (Gironde)
- FOS-SUR-MER + PORT_DE_BOUC (Rhône)
- MARSEILLE (Rhoône)

Note that the tide gauge choice for the Rhône river is not as clear as it is for the other rivers. This is because the nearest tide gauge to the mouth _FOS-SUR-MER_ has a limited time series, which is extended to the present by the next nearest tide gauge at _PORT_DE_BOUC_. However, if we re willing to search a bit further down the coast we find the tide gauge for _MARSEILLE_, which has a much longer time series.

### Oceanographic variables

The other variables used to investigate the plumes were taken from the GLORYS12V1 reanalysis product: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

