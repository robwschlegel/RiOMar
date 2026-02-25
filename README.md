# RiOMar

A collection of code used to process and analyse data relating to the [RiOMar project](https://riomar.lsce.ipsl.fr/).

Note that the majority of the code in the __/func__ folder was originally written by [Louis Terrats](https://github.com/louis-terrats) under an MIT license and may be found on [GitHub](https://github.com/louis-terrats/myRIOMAR_dev) where it may be installed and used as a standalone module.

The code is being repurposed in this repository as a continuation of the work performed by Louis, and with an eye towards the next two years of development. Terminating with the completion of the RiOMar project itself.

The primary objective of this repository is to provide a clear layout of the workflows used to verify the integrity of a number of data products, as well as any analyses used for publications.

The workflow may be found in the __/code__ folder, with the order of operations given in the names of the scripts.

# Output

The output of the code in this repository is primarily a series of figures and data files used in publications.

## Data sources

### _In situ_ networks

Multiple relevant variables are accessible via the following _in situ_ networks.

#### SOMLIT

There is a public request for access form available directly on the [SOMLIT](https://www.somlit.fr/demande-de-donnees/) website.
The following data were requested for all available SOMLIT sites:
- Type de série: SÉRIE HYDRO
- Paramètres: T (temperature), S (salinity), COP (POC; particulate organic matter), MES (SPM; suspended particulate matter), CHLA (chlorophyll a)
- Choix de la marée: PLEINE MER, BASSE MER, PAS DE MARÉE
- Profondeur: SURFACE, NIVEAU INTERMÉDIARE, FOND
- Dates: 1998-01-01 to 2025-12-31

#### REPHY

The full REPHY dataset from 1987 - 2022 is openly available on [SEANOE](https://www.seanoe.org/data/00361/47248/).
These data were downloaded in their entirety then processed for further use in this project.

#### OSR

The [Observatoire des Sediments du Rhône](https://bdoh.inrae.fr/OBSERVATOIRE-DES-SEDIMENTS-DU-RHONE/) 
makes available a long list of variables for water measurements made at stations along the Rhône river and its tributaries.

One must request access to these data via their online form.
Once granted, the data used in this project were:
- ARLES/CMES-2
- ARLES/FMES
- BEAUCAIRE/FMES
- MESURHO/CMES
- MESURHO/CMES-2
- MESURHO/TURB
- MESURHO/TURB-2
- MESURHO/TURB-3

### River discharge

#### HydroPortail

The primary public source of river discharge data is [HydroPortail](https://hydro.eaufrance.fr/).
On this site one can create lists of river flow gauges of interest, query the sites, and download data.
It is also possible to access this database via an API at [hub'eau](https://hubeau.eaufrance.fr/page/api-ecoulement#/ecoulement/getResultats).

Note that while many of the river gauges have data going back decades, 
we only use data from 1998-01-01 onward as this is the beginning of the SEXTANT data product used for the plume analyses.

NB: In order to get an idea of the range of dates available for a station, click the 'calendrier' tab..

#### La Seine

- La Seine à Vernon - H320 0001 01: 1998 - 2006
  - This appears to be one of the closest flow gauges to the mouth of the Seine.
- La Seine à Vernon - H320 0001 04: 2006 - 2025
  - Active flow gauge to present date.
- La Seine à Poses - H322 0110 03: 1998-2006
  - Closer than Vernon, but only data to 2006.
- L'Orne à May-sur-Orne - I362 1010 01: 1998-2025
  - The Orne river is a notable contributor to the Bay of Seine.
  - There are three flow gauges on this river, this is the closest to the mouth.

#### La Loire

- La Loire à Saumur - L800 0010 20: 1998 - 2025
  - Far from the mouth, but the only available flow gauge on the river.
- Le Lay à la Bretonnière-la-Claye - N351 1610: 2003 - 2025
  - A non-negligble input into the waters of Southern Brittany.
  - This area frequently connects to the Loire plume.
- Le Lay à Mareuil-sur-Lay-Dissais - N330 1610 10: 1998 - 2025
  - Further from the mouth of the river, but with a complete time series.
- La Vilaine à Rieux - J930 0611 01: 2002 - 2025
  - Another important river contributing to the southern coast of Brittany.
  - This is the closest flow gauge to the mouth.
- La Vilaine à Guipry - J770 0610 02: 1998 - 2025
  - Further from the mouth, but complete time series.

#### Gironde

- La Garonne à Marmande - O909 0010 01: 1998  - 2025
  - One of the two major contributing rivers to the Gironde.
  - This gauge is far from the mouth, but appears to be the closest.
- La Garonne à Tonneins - O900 0010 02: 1998 - 2025
  - Another complete time series for comparison.
- La Dordogne à Pessac-sur-Dordogne - P555 0010 01: 1998 - 2025
  - The second of the two major contributors to the Gironde.

##### Rhône

Multiple stations were accessed to summarise the river flow as close to the mouth of the Rhône as possible.

- Le Rhône à Arles - V730 0003 02: 1998 - 2023
  - This was used as the closest approximation of the river flow for the Grand Rhône.
- Le Rhône à Fourques - V730 0002 02: 1998 - 2023
  - This was used as the closest approximation of the river flow for the Petit Rhône.
- Le Rhône à Tarascon - V720 0010 02: 2024-2025
  - Used to update the time series from 2024 to present.
- Le Rhône à Beaucaire - V720 0005 01:
  - Not used as it does not appear to have any data.

Note that if using the station for the Rhône at Beaucaire or Tarascon, in order to calculate the discharge at the mouths of the Grand Rhône and 
Petit Rhône one should apportion 11% of the river flow to the Petit Rhône up until 2012-05-28, after which this should be reduced to 10%
of the river flow up until the present date.

#### Ifremer

Another source of river discharge data is made available by [Ifremer](https://co-en.ifremer.fr/eulerianNetwork?lang=en&contextId=386).

Unfortunately this site is largely unresponsive and I've never managed to download data from it,
but if one were to attempt to do so, these are the product IDs to use.

- Rhône : IF000132 & EXSC0042
- Seine : IF000134 & IF000509 & EXSC0048
- Loire :  IF000125 & EXSC0002
- Gironde = Dordogne + Garonne
  - Garonne :  IF000234 & EXSC0003
  - Dordogne :  IF000507 & EXSC0018

Note that the IF* products are historic data, and EXS* are the current/ongoing data streams.

### Waves

_In situ_ measurements for waves may be found at [Candhis](https://candhis.cerema.fr/index.php).

### Wind

The dedicated CMEMS wind reanalysis product was used in this study: https://data.marine.copernicus.eu/product/WIND_GLO_PHY_L4_NRT_012_004/description

### Tides

Tide gauge data were downloaded via web form at: https://data.shom.fr/donnees/refmar/download

The id : name used were:
- 4 : LE_HAVRE (Seine)
- 37 : SAINT-NAZAIRE (Loire)
- 15 : PORT-BLOC (Gironde)
- 719 : FOS-SUR-MER + 720 : PORT_DE_BOUC (Rhône)
- 524 : MARSEILLE (Rhône)

Note that the tide gauge choice for the Rhône river is not as clear as it is for the other rivers. 
This is because the nearest tide gauge to the mouth _FOS-SUR-MER_ has a limited time series, which is extended to the present by the next nearest tide gauge at _PORT_DE_BOUC_. 
However, if we re willing to search a bit further down the coast we find the tide gauge for _MARSEILLE_, which has a much longer time series. 
An explanation for how to access these data is available [here](https://refmar.shom.fr/donnees-refmar-sur-data.shom.fr/telechargement-des-donnees)

### Oceanographic variables

The other variables used to investigate the plumes were taken from the GLORYS12V1 reanalysis product: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

