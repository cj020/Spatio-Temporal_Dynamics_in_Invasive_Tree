## Data Description

Tree cover estimates at a 30 m $\times$ 30 m spatial resolution grid from 1985 to 2020 are obtained from the Rangeland Analysis Platform.

### Spatial Coverage

-   Longitude range: -100.6998 to -100.1980\
-   Latitude range: 42.80019 to 42.99907

### Data Structure

Each CSV file corresponds to one year and contains observations on a fixed spatial grid with **1,371,134 observations per file**.

-   From **1985 to 2018**, each CSV file contains **20 columns**.\
-   From **2019 to 2020**, each CSV file contains **24 columns**.

### Variables

The dataset includes the following variables:

-   **Longitude, Latitude**: Spatial coordinates of each grid cell\
-   **trees**: Percentage of tree cover within each grid cell\
-   **DEM_30m**: Elevation (Digital Elevation Model at 30 m resolution)\
-   **Stream**: Indicator of whether a stream is present in the grid cell\
-   **Climatic variables**: Additional environmental covariates (varying slightly across years, which accounts for the increase in the number of columns after 2019)
