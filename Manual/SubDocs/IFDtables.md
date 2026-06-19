# Creating the IFD tables

The work below is stored in the ```BureauIFD``` folder where Sunwater scripts are on the server. Repositories needed include: json, geopandas, rasterio, xarray, rioxarray, , and 

Grids of the IFD depths have been downloaded (on the 26th June 2024) from the Bureau's IFD portal across all of Sunwater's dam cacthments - except for the big ones: Burdekin, Paradise, Fairbairn, and Beardmore. The ASCII files were merged into a single IFD database (NetCDF file) called, ```IFDdatabase/BureauIFD.nc```, using the ```CreateNC.py``` script. The NetCDF file containes three-dimentional dataframes -- ```x``` (degrees), ```y``` (degrees), and ```duration``` (minutes) -- for each AEP (1 in X AEP). The projection is: WGS 84 (EPSG: 4326). The ASCII files have not been stored on the network due to the large number of files. Nevertheless, the NetCDF file can be loaded into QGIS. Storm durations span 30 minutes to 5 days and storm magnitudes span 1 in 2 to 1 in 2000 AEP. 

The ```CollateIFDtable.py``` script creates the IFD tables for reading into Bryan. The IFD database is used in this script to extract the average rainfall depth across each subcatchment in the model. The path to the shapefile of subcatchments must be specified in the ```paths.json``` file; an example is shown below:
```json
{
  "BureauIFDdata": "IFDdatabase\\BureauIFD.nc",
  "callide": "Dams\\Callide\\callide-SubCatchments.shp"
}
```
Any number of dam paths can be appended to the dictionary. Then the path is applied in the ```CollateIFDtable.py``` scipt as follws, where the ```subcatchmentcolumn``` specifies the header name in the shapefile with the subcatchment names (values must be integers like used in URBS):
```python
 # Set the dam to assess
dam = 'callide'
subcatchmentcolumn = 'ID_Number'
```
Then running the script will create the tables in the same folder as the shapefile. As always, do some spot checks on the results to make sure the averaged depths are as expected.