# FluxCT V2.0

Welcome to the batch code for the FluxCT webtool. The following downloadable script allows the user to download the script for to batch run Kepler ID data. 

Inputs: 
<ul>
<li> A .csv file of KIC IDs, using a comma as the separator, with the title 'id' in the column. </li>
<li> A directory path for where the script and id file will be kept. </li>
<li> A directory path for the fits files, plots, and final data file to be sent. </li>
</ul>
  
Outputs: 
<ul>
<li> A .fits file for each of the stars downloaded from lightkurve. </li>
<li> A .png file for the plot produced by FluxCT. </li>
<li> A .csv file containing the output data from FluxCT. </li> 
</ul>

Any problems, concerns or suggestions should be sent to: Zilin.Dong@vanderbilt.edu and jessica.s.stasik@vanderbilt.edu for the most timely review and reply. 

Update: 
>1. TESS datasets are enormously big, around 60MB (Way larger than 1~2 Mb, the size of Kepler dataset)
>2. Does MastDownload dir data need to be clean?
>3. kic_fits.fits dataset is removed after plot and r1~r7 is calculated/generated.
>4. plot2.jpg is removed after reading it from memory
