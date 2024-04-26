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
  We upgraded brand new version FluxCT V2.0, a web tool for identifying contaminating flux in
Kepler and TESS target pixel files due to secondary visual sources. FluxCT V2.0 focuses on enhancing
functionality, user experience, and data processing capabilities. We resolved existing issues in using
target pixel files from TESS, allowing the user to search any point object observed by TESS and acquire
parameters for the star and any possible contaminating visual pairs. Issues relating to the running8
of the web tool have been resolved, and a batch code for TESS is now available on the companion9
GitHUb.   
