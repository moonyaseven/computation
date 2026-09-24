Ecosystem Service Computation Codes

This repository contains the Google Earth Engine (GEE) codes used to generate annual ecosystem service datasets for China from 2005 to 2024 at a target spatial resolution of 30 m.

The repository currently includes codes for:

Net Primary Productivity (NPP)
Net Ecosystem Productivity (NEP)
Carbon Sequestration (CS)
Soil Retention (SR)
Water Yield (WY)
Code environment

The calculations were implemented in the Google Earth Engine JavaScript API.

The scripts can be copied directly into the Google Earth Engine Code Editor:

https://code.earthengine.google.com/

Input data

The calculations use remote-sensing, meteorological, land-cover, soil, and other ancillary datasets described in the associated publication.

Some input datasets were uploaded to user-specific Google Earth Engine assets. Therefore, asset paths such as:

ee.Image('users/username/...')

or

ee.Image('projects/project_name/assets/...')

may need to be replaced with the corresponding datasets in the user's own GEE account.

The solar radiation dataset used in the calculations is subject to data-use restrictions and is therefore not redistributed through this repository. Users should obtain authorized access to the dataset before reproducing calculations that require solar radiation input.

Spatial and temporal coverage
Study area: China
Period: 2005–2024
Temporal resolution: Annual
Target spatial resolution: 30 m
Output coordinate reference system: EPSG:4326
Output format: GeoTIFF

Continuous variables were resampled using bilinear interpolation, while categorical variables were resampled using nearest-neighbor interpolation where applicable.

Running the codes
Open the Google Earth Engine Code Editor.
Copy the required .js script into the Code Editor.
Replace user-specific GEE asset paths where necessary.
Set the required year range and study region.
Run the script.
Start the generated export tasks in the Tasks panel.

The scripts generally process the data by administrative region and year to reduce the computational load in Google Earth Engine.

Citation

If you use these codes, please cite the associated dataset publication.

Citation information will be added after publication.

Data availability

This repository contains the computational codes only. The generated ecosystem service datasets are archived separately. Please refer to the associated publication for detailed information on data access and source datasets.

Contact

For questions regarding the codes or reproduction of the dataset, please contact the authors of the associated publication.
