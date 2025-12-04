# Code and data for "Movement and space use patterns reveal greater vulnerability to environmental changes and poaching in the Beale’s eyed turtle (*Sacalia bealei*)"

This repository contains the data and code used in the analysis of the spatial patterns, habitat selection and microhabitat uses of Beale's eyed turtles (*Sacalia bealei*) in Hong Kong Special Administrative Region:

> Chan, W. S., Fong J. J., Bonebrake, T. C., & Sung, Y. H.  (2024) Movement and space use patterns reveal greater vulnerability to environmental changes and poaching in the Beale’s eyed turtle (*Sacalia bealei*). *In submission*.

## Data

All data are contained in the directory `data/`. The original data containing the spatial coordinates of the animals are deposited on Movebank (ID: 5056848767). Due to conservation concerns, accessing to the original data on Movebank is restricted to the project team only. All data in the repository were pre-processed to spatial estimates used to fit subsequent models.

-   `SABE_biometry.csv` contains biometric measurements of 25 monitored individuals

-   `SABE_movement_mcp.csv` contains the calculated home range sizes of each individual across the entire study period

-   `SABE_movement_stdist.csv` contains the 100% stream distances travelled by each individual across the entire study period

-   `SABE_movement_mcp_season.csv` contains the calculated home range sizes of each individual for each respective season

-   `SABE_movement_dist.csv` contains movement trajectories of each individual *S. bealei*

-   `SABE_habitatProp_studyarea.csv` summarises composition of aquatic habitat types found within the study area over the entire study period

-   `SABE_habitatProp_studyarea_season.csv` summarises composition of aquatic habitat types found within the study area during each respective season

-   `SABE_habitatProp_mcp.csv` summarises composition of aquatic habitat types found within area encompassing the home range of each individual turtle over the entire study period

-   `SABE_habitatProp_mcp_season.csv` summarises composition of aquatic habitat types found within area encompassing the home range of each individual turtle during each respective season

-   `SABE_habitatProp_loc.csv` summarises composition of aquatic habitat types among all used locations of each individual turtle over the entire study period

-   `SABE_habitatProp_loc_season.csv` summarises composition of aquatic habitat types among all used locations of each individual turtle during each respective season

-   `SABE_habitatselection.csv` contains paired samples of the microhabitat features used by and available to each individual *S. bealei*

Description of the variables are listed in the table below:

| variable    | description                                                                                                                                                                                                            |
|---------------------|---------------------------------------------------|
| event.id    | unique id assigned for each relocation event during the radiotelemetry process                                                                                                                                         |
| turtle.id   | unique id assigned for each *Sacalia bealei*                                                                                                                                                                           |
| frequency   | radio-transmitting frequency of the respective radio transmitter deployed                                                                                                                                              |
| radio.id    | unique id denoted by `turtle.id`\_`frequency`                                                                                                                                                                          |
| sex         | reproductive sex of the respective individual                                                                                                                                                                          |
| cl          | straight-line carapace length (mm) of the respective individual                                                                                                                                                        |
| pl          | straight-line plastron length (mm) of the respective individual                                                                                                                                                        |
| wt          | body mass (g) of the respective individual                                                                                                                                                                             |
| date        | date when the relocation event was recorded                                                                                                                                                                            |
| time        | datetime (UTC + 08:00) when the relocation event was recorded                                                                                                                                                          |
| session     | time of day (`Day` or `Night`) when the relocation event was recorded                                                                                                                                                  |
| month       | month for which the last point of the trajectory path was sampled from / the relocation event was recorded at                                                                                                          |
| year        | year for which the last point of the trajectory path was sampled from / the relocation event was recorded at                                                                                                           |
| season      | season for which the home range size was characterized for / the last point of the trajectory path was sampled from / the composition of the aquatic habitat was summarised for / the relocation event was recorded at |
| HR_Type     | method used for the estimation of home range type                                                                                                                                                                      |
| MCP_Lv      | percentage of the points used for the drawing of the MCP                                                                                                                                                               |
| Area_m2     | home range size (m<sup>2</sup>) for the corresponding individuals in the respective season                                                                                                                                    |
| Pts         | number of relocation points recorded / used for the estimation of the home range size                                                                                                                                  |
| dx          | change of position in the x direction between two successive relocations                                                                                                                                               |
| dy          | change of position in the y direction between two successive relocations                                                                                                                                               |
| dist        | displacement distance between two successive relocations                                                                                                                                                               |
| dt          | time interval (s) between successive relocations                                                                                                                                                                       |
| dailyDist   | daily displacement distance                                                                                                                                                                                            |
| Pool        | proportion of pool among other aquatic habitats sampled within the selected areas in the predefined period                                                                                                             |
| Riffle      | proportion of riffle among other aquatic habitats sampled within the selected areas in the predefined period                                                                                                           |
| Run         | proportion of run among other aquatic habitats sampled within the selected areas in the predefined period                                                                                                              |
| visible     | boolean value (`TRUE` or `FALSE`) identifying if the tracking individual was visually identified during the relocation event                                                                                           |
| behaviour   | activity identified in the tracking animal during the relocation event                                                                                                                                                 |
| elev        | elevation of the sampled location                                                                                                                                                                                      |
| habitat     | habitat (`pool`, `run`, `riffle`, `forested upland`, `open upland`) of the sampled location                                                                                                                            |
| waterDepth  | distance (cm) between water surface and bottom of the stream in the sampled location                                                                                                                                   |
| streamWidth | width (m) of the flowing channel in the sampled location                                                                                                                                                               |
| gravel      | percentage covers of gravel in the sampled location                                                                                                                                                                    |
| pebble      | percentage covers of pebble in the sampled location                                                                                                                                                                    |
| cobble      | percentage covers of cobble in the sampled location                                                                                                                                                                    |
| boulder     | percentage covers of boulder in the sampled location                                                                                                                                                                   |
| litter      | percentage of area covered by leaf litter in the sampled location                                                                                                                                                      |
| canopy      | percentage of forest canopy in the sampled location                                                                                                                                                                    |
| status      | binary expression of whether the listed microhabitat feature was selected by the respective individual at the corresponding time period                                                                                |

## Scripts

All R codes used to analyse the data are available in the directory `script/`.

-   `SABE_SpatialAnalysis_01_Preprocessing.R` contains codes to clean the original data and summarise the raw data containing the spatial coordinates into used data formats for the subsequent analyses

-   `SABE_SpatialAnalysis_02_Movement.R` contains code to fit several Generalised Linear Mixed Models (GLMMs) regarding the summarised movement characterisitcs

-   `SABE_SpatialAnalysis_03_CompositionalAnalysis.R` contains code to fit compositional analysis to investigate habitat selection pattern of *S. bealei*

-   `SABE_SpatialAnalysis_04_RSF.R` contains code to fit Resource Selection Functions (RSFs) pertaining the microhabitat selection of *S. bealei*


## Contact Information

Questions regarding the code or the data should be directed to the project author Wing Sing Chan at [wschan1021\@gmail.com](mailto:wschan1021@gmail.com).

## License

This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution - Non-Commercial 4.0 International License</a>.<br />
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commos License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a>
