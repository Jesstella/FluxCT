'''
This script allows the user to create the data products and plots from FluxCT for multiple targets without using the webtool. 

INPUTS: 
 --> The user will need to specify the directory path for code and plots. 
 --> An input table of either KIC IDs with the format 'id' as the column name and comma separators. 

OUTPUTS: 
 --> Downloaded fits file for the star from lightkurve to the plots directory. 
 --> Downloaded plot created by FluxCT to the plots directory. 
 --> Downloaded FluxCT data products to a .csv file.
'''

# Import modules
import time 
import pandas as pd  
from lightkurve import search_targetpixelfile
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.io import fits 
import numpy as np
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import math
from ast import literal_eval
import ssl

try:
     _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default,
    pass
else:
    # Handle target environment that doesn't support HTTPS verification,
    ssl._create_default_https_context = _create_unverified_https_context

# Paths and files 
code_file_path = '/users/jess/FluxCT/' # USER INPUT - Code directory 
plot_path = '/users/jess/FluxCT/plots/' # USER INPUT - Plot directory 
identifiers = pd.read_csv(code_file_path + 'test_kic_list.csv') # USER INPUT - KIC ID list  
identifiers = list(identifiers['id'])
print('There are ' + str(len(identifiers)) + ' stars in the sample.') 

# Creating open lists for collecting data
kic_list = []
gaia_source_list = []
ra_list = []
dec_list = []
g_mag_list = []
ruwe_list = []
flux_list = []
not_found = []

# Starting time counter 
t0 = time.time()
a = 0

# Beginning process for i in the list 
for i in identifiers:
   
    # Printing explanatory data  
    kic = i[4:]
    print('\n********** KIC ' + str(kic) + ' – Star Number ' + str(a) + ' **********') 
   
    # Searching for tpf with lightkurve
    tpf = search_targetpixelfile('KIC ' + str(kic), author='Kepler', cadence='short').download()
    tpf_one = tpf[0]
    tpf_one.to_fits(plot_path + 'kic_' + str(kic) + '_fits.fits')

    # Finding tpa
    array = tpf_one.pipeline_mask 

    # Finding the count of each pixel in the rows and columns of the mask. 
    row_count = np.count_nonzero(array, axis=1)
    col_count = np.count_nonzero(array, axis=0) 

    # Columns = col are counting the x coordinate from left to right. 
    # This is the lower limit in pixel number for the left hand side of the array
    # By adding 1 to the result I am removing the index to 0.
    count_x = 0
    if col_count[0] != 0:
        count_x = 0
    else:
        for i in range(len(col_count)):
            if col_count[i] == 0:
                count_x = i
            else:
                count_x = count_x + 1
                break
     
    # Rows are counting the y coordinate from bottom to top.
    # This is the lower limit in pixel number for the bottom of the array
    # By adding 1 to the result I am removing the index to 0.
    count_y = 0
    if row_count[0] != 0:
        count_y = 0
    else:
        for j in range(len(row_count)):
            if row_count[j] == 0:
                count_y = j 
            else:
                count_y = count_y + 1
                break
    
    # Finding the top right pixel
    array_2  = array[np.ix_(~np.all(array == False, axis=1), ~np.all(array == False, axis=0))]
    tr_y = array_2.shape[0]  
    tr_x = array_2.shape[1]  

    # Top left pixel position 
    tl_x = count_x - 0.5
    tl_y = count_y + tr_y - 1 + 0.5

    # Bottom right pixel position 
    br_x = count_x + tr_x - 1 + 0.5
    br_y = count_y - 0.5

    # Bottom left pixel position
    # Starts from the first pixel with a non-zero value so only include count_x and count_y 
    bl_x = count_x - 0.5
    bl_y = count_y - 0.5

    # Top right pixel position
    # Already calaculated these so just add the 0.5
    tr_x = count_x + tr_x - 1 + 0.5
    tr_y = count_y + tr_y - 1 + 0.5 

    # Pulling image to plot
    tpf_data = fits.open(plot_path + 'kic_' + str(kic) + '_fits.fits')
    image = tpf_data[1].data
    image = image['FLUX'][0]
    wcs = WCS(tpf_data[2].header)

    # Finding the corners of the aperture mask
    tl = tpf_one.wcs.pixel_to_world(tl_x, tl_y)
    tr = tpf_one.wcs.pixel_to_world(tr_x, tr_y)
    bl = tpf_one.wcs.pixel_to_world(bl_x, bl_y)
    br = tpf_one.wcs.pixel_to_world(br_x, br_y)

    # Converting the corners of the aperture mask to coordinates for Gaia
    top_left = wcs.world_to_pixel(tl)
    top_right = wcs.world_to_pixel(tr)
    bottom_left = wcs.world_to_pixel(bl)
    bottom_right = wcs.world_to_pixel(br)

    # Coordinates to search in Gaia 
    tr_ra = tr.ra.deg
    tr_dec = tr.dec.deg
    tl_ra = tl.ra.deg
    tl_dec = tl.dec.deg
    br_ra = br.ra.deg
    br_dec = br.dec.deg
    bl_ra = bl.ra.deg
    bl_dec = bl.dec.deg

    # Creating Gaia call
    polygon = str(br_ra) + ', ' + str(br_dec) + ', ' + str(bl_ra) + ', ' + str(bl_dec) + ', ' + str(tl_ra) + ', ' + str(tl_dec) + ', ' + str(tr_ra) + ', ' + str(tr_dec)
    columns = 'source_id, ra, dec, phot_g_mean_mag, ruwe, phot_g_mean_flux'
    polygon_top10query_base = """SELECT
    {columns}
    FROM gaiaedr3.gaia_source
    WHERE 1=CONTAINS(
            POINT(ra, dec), 
            POLYGON({polygon}))
    """

    # Querying Gaia
    polygon_top10query = polygon_top10query_base.format(columns=columns, 
                          polygon=polygon)

    polygon_top10query_job = Gaia.launch_job_async(polygon_top10query)

    polygon_top10query_results = polygon_top10query_job.get_results()

    # Pulling data from Gaia
    plot_source = list(polygon_top10query_results['source_id'])
    plot_ra = list(polygon_top10query_results['ra'])
    plot_dec = list(polygon_top10query_results['dec'])
    phot = list(polygon_top10query_results['phot_g_mean_mag'])
    ruwe = list(polygon_top10query_results['ruwe']) 
    flux = list(polygon_top10query_results['phot_g_mean_flux']) 

    # Initializing lists
    plot_source_order = []
    plot_ra_order = []
    plot_dec_order = []
    plot_ruwe_order = []
    plot_phot_order = np.sort(phot)
    plot_flux_order = []

    # Putting the data in the correct order for plotting
    for i in range(len(phot)):
        if plot_phot_order[i] in phot:
            index = phot.index(plot_phot_order[i])
            plot_source_order.append(plot_source[index]) 
            plot_ra_order.append(plot_ra[index])
            plot_dec_order.append(plot_dec[index])
            plot_ruwe_order.append(ruwe[index]) 
            plot_flux_order.append(flux[index]) 

    # Saving final data lists for the output file
    kic_list.append(kic)
    ra_list.append(plot_ra_order) 
    dec_list.append(plot_dec_order) 
    gaia_source_list.append(plot_source_order) 
    ruwe_list.append(plot_ruwe_order) 
    flux_list.append(plot_flux_order) 

    # Removing nans from the photometry
    plot_phot_order = [-999 if math.isnan(x) else x for x in plot_phot_order]
    g_mag_list.append(plot_phot_order) 

    # Determining the coordinates of the companions for plotting
    companions = SkyCoord(plot_ra_order, plot_dec_order, unit='deg')
    companions_to_plot = wcs.world_to_pixel(companions)

    # Making beautiful plot! 
    if len(companions_to_plot[0]) == 0: 
        not_found.append(kic) 
    else: 

        # Setting figure
        fig = plt.figure(figsize=(18, 15))
        fig.add_subplot(111, projection = wcs) 

        # Plotting corners and box around TPA 
        top_line_x = [tl_x, tr_x]
        top_line_y = [tl_y, tr_y]
        plt.plot(top_line_x, top_line_y, linewidth=3, color='white')
        bottom_line_x = [bl_x, br_x]
        bottom_line_y = [bl_y, br_y]
        plt.plot(bottom_line_x, bottom_line_y, linewidth=3, color='white')
        left_line_x = [bl_x, tl_x]
        left_line_y = [bl_y, tl_y]
        plt.plot(left_line_x, left_line_y, linewidth=3, color='white')
        right_line_x = [br_x, tr_x]
        right_line_y = [br_y, tr_y]
        plt.plot(right_line_x, right_line_y, linewidth=3, color='white')

        # Plotting magnitudes 
        for i in range(len(phot)):
            try:
                plt.text(companions_to_plot[0][i], companions_to_plot[1][i], str(round(plot_phot_order[i], 3)), color='#dd1c77', fontsize=25)
            except:
                continue

        # Plotting target star and companions
        plt.scatter(companions_to_plot[0][1:], companions_to_plot[1][1:], marker='*', s=2000, color='white', edgecolor='black')
        plt.scatter(companions_to_plot[0][0], companions_to_plot[1][0], marker='*', s=2000, color='pink', edgecolor='black')
        
        # Plotting corner text 
        plt.text(tr_x+0.2, tr_y+0.2, 'RA = ' + str(round(tr.ra.deg, 4)) + '\nDec = ' + str(round(tr.dec.deg, 4)), fontsize=15, backgroundcolor='white') 
        plt.text(tl_x-1.8, tl_y+0.2, 'RA = ' + str(round(tl.ra.deg, 4)) + '\nDec = ' + str(round(tl.dec.deg, 4)), fontsize=15, backgroundcolor='white') 
        plt.text(br_x+0.2, br_y+0.2, 'RA = ' + str(round(br.ra.deg, 4)) + '\nDec = ' + str(round(br.dec.deg, 4)), fontsize=15, backgroundcolor='white') 
        plt.text(bl_x-1.8, bl_y+0.2, 'RA = ' + str(round(bl.ra.deg, 4)) + '\nDec = ' + str(round(bl.dec.deg, 4)), fontsize=15, backgroundcolor='white') 
       
        # Plotting corners
        plt.scatter(tr_x, tr_y, s=200, marker='X', color='white') # Top right
        plt.scatter(tl_x, tl_y, s=200, marker='X', color='white') # Top left
        plt.scatter(br_x, br_y, s=200, marker='X', color='white') # Bottom right
        plt.scatter(bl_x, bl_y, s=200, marker='X', color='white') # Bottom left

        # Setting axes for ticks
        ax = plt.gca()

        # Plotting axes labels, titles, and images
        plt.ylabel('DEC [degrees]', fontsize=20)
        plt.xlabel('RA [hourangle]', fontsize=20)
        plt.imshow(image, origin='lower', cmap='RdPu_r', alpha=1)
        plt.imshow(array, origin='lower',  cmap='binary_r', alpha=0.2)
        plt.title('KIC ' + str(kic), fontsize=20)

        # Setting tick parameters 
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)

        # Plotting grid and saving the file
        plt.grid(axis = 'both', color='grey', ls = ':', linewidth=6)
        plt.savefig(plot_path + str(kic) + '.png')
        a = a + 1
        plt.close()

# Timing code 
t1 = time.time()
total = t1 - t0
print('The total time to create these plots is ' + str(total/60) + ' minutes.') 

# Saving data file
data = zip(kic_list, ra_list, dec_list, gaia_source_list, g_mag_list, ruwe_list, flux_list) 
header = ['kic', 'ra', 'dec', 'gaia_source_list', 'g_mag_list', 'ruwe', 'flux']
df = pd.DataFrame(data=data, columns=header) 
df.to_csv(code_file_path + 'fluxct_data.csv', index=False) 
