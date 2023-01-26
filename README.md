# FluxCT

Welcome to the batch code for the FluxCT webtool. The following downloadable script allows the user to download the script for to batch run Kepler ID data. 

Inputs: 
--> A .csv file of KIC IDs, using a comma as the separator, with the title 'id' in the column.
--> A directory path for where the script and id file will be kept. 
--> A directory path for the fits files, plots, and final data file to be sent. 

Outputs: 
--> A .fits file for each of the stars downloaded from lightkurve. 
--> A .png file for the plot produced by FluxCT. 
--> A .csv file containing the output data from FluxCT. 
