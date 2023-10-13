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
from io import StringIO
from io import BytesIO
import sys
import io
import lightkurve
#kic = 'KIC 3430893'

def do_calculation(kic):
    # Set up lists for parameters
    kic_list = []
    gaia_source_list = []
    ra_list = []
    dec_list = []
    not_found = []
    companion_mag_list = []
    flux_list = []
    companion_flux_list = []

    if kic[:3] == 'KIC':
        tpf = search_targetpixelfile(kic, author='Kepler').download()
    elif kic[:3] == 'TIC':
        tpf = search_targetpixelfile(kic, author='SPOC').download()
    else:
        print('nope')

    m1 = tpf[0]
    m2 = m1.to_fits('kic_fits.fits', overwrite=True)

    array = m1.pipeline_mask

    row_count = np.count_nonzero(array, axis=1)
    col_count = np.count_nonzero(array, axis=0)

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

    array_2  = array[np.ix_(~np.all(array == False, axis=1), ~np.all(array == False, axis=0))]
    tr_y = array_2.shape[0]
    tr_x = array_2.shape[1]

    tl_x = count_x - 0.5
    tl_y = count_y + tr_y - 1 + 0.5

    br_x = count_x + tr_x - 1 + 0.5
    br_y = count_y - 0.5

    bl_x = count_x - 0.5
    bl_y = count_y - 0.5

    tr_x = count_x + tr_x - 1 + 0.5
    tr_y = count_y + tr_y - 1 + 0.5

    m2 = fits.open('kic_fits.fits')
    image = m2[1].data
    image = image['FLUX'][0]
    wcs = WCS(m2[2].header)

    tl = m1.wcs.pixel_to_world(tl_x, tl_y)
    tr = m1.wcs.pixel_to_world(tr_x, tr_y)
    bl = m1.wcs.pixel_to_world(bl_x, bl_y)
    br = m1.wcs.pixel_to_world(br_x, br_y)

    top_left = wcs.world_to_pixel(tl)
    top_right = wcs.world_to_pixel(tr)
    bottom_left = wcs.world_to_pixel(bl)
    bottom_right = wcs.world_to_pixel(br)

    tr_ra = tr.ra.deg
    tr_dec = tr.dec.deg
    tl_ra = tl.ra.deg
    tl_dec = tl.dec.deg
    br_ra = br.ra.deg
    br_dec = br.dec.deg
    bl_ra = bl.ra.deg
    bl_dec = bl.dec.deg

    polygon = str(br_ra) + ', ' + str(br_dec) + ', ' + str(bl_ra) + ', ' + str(bl_dec) + ', ' + str(tl_ra) + ', ' + str(tl_dec) + ', ' + str(tr_ra) + ', ' + str(tr_dec)
    columns = 'source_id, ra, dec, phot_g_mean_mag, ruwe, phot_g_mean_flux'

    polygon_top10query_base = """SELECT
    {columns}
    FROM gaiaedr3.gaia_source
    WHERE 1=CONTAINS(
            POINT(ra, dec),
            POLYGON({polygon}))
    """
    polygon_top10query = polygon_top10query_base.format(columns=columns,
                          polygon=polygon)

    polygon_top10query_job = Gaia.launch_job_async(polygon_top10query)

    polygon_top10query_results = polygon_top10query_job.get_results()
    polygon_top10query_results

    plot_source = list(polygon_top10query_results['source_id'])
    plot_ra = list(polygon_top10query_results['ra'])
    plot_dec = list(polygon_top10query_results['dec'])
    phot = list(polygon_top10query_results['phot_g_mean_mag'])
    ruwe = list(polygon_top10query_results['ruwe'])
    flux = list(polygon_top10query_results['phot_g_mean_flux'])

    plot_source_order = []
    plot_ra_order = []
    plot_dec_order = []
    plot_ruwe_order = []
    plot_flux_order = []
    plot_phot_order = np.sort(phot)

    for i in range(len(phot)):
        if plot_phot_order[i] in phot:
            index = phot.index(plot_phot_order[i])
            plot_source_order.append(plot_source[index])
            plot_ra_order.append(plot_ra[index])
            plot_dec_order.append(plot_dec[index])
            plot_ruwe_order.append(ruwe[index])
            plot_flux_order.append(flux[index])

    kic_list.append(kic)
    ra_list.append(plot_ra_order)
    dec_list.append(plot_dec_order)
    gaia_source_list.append(plot_source_order)
    ruwe_list = plot_ruwe_order
    flux_list = plot_flux_order

    plot_phot_order = [-999 if math.isnan(x) else x for x in plot_phot_order]
    g_mag_list = plot_phot_order
    primary_g_mag = g_mag_list[0]
    companion_mag_list = list(g_mag_list[1:])
    try:
        mag_diff_list = [abs(primary_g_mag - k) for k in companion_mag_list]
    except:
        mag_diff_list = []

    flux_list = [-999 if math.isnan(x) else x for x in flux_list]
    primary_g_flux = flux_list[0]
    companion_flux_list = list(flux_list[1:])
    try:
        flux_ratio_list = [(k / (primary_g_flux + k)) for k in companion_flux_list]
        percentage_flux = [h*100 for h in flux_ratio_list]
    except:
        flux_ratio_list = []
    flux_contamination_total = (sum(companion_flux_list) / (primary_g_flux + sum(companion_flux_list))) * 100

    companions = SkyCoord(plot_ra_order, plot_dec_order, unit='deg')
    companions_to_plot = wcs.world_to_pixel(companions)

    if len(companions_to_plot[0]) == 0:
        not_found.append(kic)
    else:

        fig = plt.figure(figsize=(18, 15))
        fig.add_subplot(111, projection = wcs)

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

        for i in range(len(phot)):
            try:
                plt.text(companions_to_plot[0][i], companions_to_plot[1][i], str(round(plot_phot_order[i], 3)), color='#dd1c77', fontsize=25)
            except:
                continue
        plt.scatter(companions_to_plot[0][1:], companions_to_plot[1][1:], marker='*', s=2000, color='white', edgecolor='black')
        plt.scatter(companions_to_plot[0][0], companions_to_plot[1][0], marker='*', s=2000, color='pink', edgecolor='black')
        plt.text(tr_x+0.2, tr_y+0.2, 'RA = ' + str(round(tr.ra.deg, 4)) + '\nDec = ' + str(round(tr.dec.deg, 4)), fontsize=15, backgroundcolor='white')
        plt.text(tl_x-1.8, tl_y+0.2, 'RA = ' + str(round(tl.ra.deg, 4)) + '\nDec = ' + str(round(tl.dec.deg, 4)), fontsize=15, backgroundcolor='white')
        plt.text(br_x+0.2, br_y+0.2, 'RA = ' + str(round(br.ra.deg, 4)) + '\nDec = ' + str(round(br.dec.deg, 4)), fontsize=15, backgroundcolor='white')
        plt.text(bl_x-1.8, bl_y+0.2, 'RA = ' + str(round(bl.ra.deg, 4)) + '\nDec = ' + str(round(bl.dec.deg, 4)), fontsize=15, backgroundcolor='white')
        plt.scatter(tr_x, tr_y, s=200, marker='X', color='white') # Top right
        plt.scatter(tl_x, tl_y, s=200, marker='X', color='white') # Top left
        plt.scatter(br_x, br_y, s=200, marker='X', color='white') # Bottom right
        plt.scatter(bl_x, bl_y, s=200, marker='X', color='white') # Bottom left

        ax = plt.gca()

        plt.ylabel('DEC [degrees]', fontsize=20)
        plt.xlabel('RA [hourangle]', fontsize=20)
        plt.imshow(image, origin='lower', cmap='RdPu_r', alpha=1)
        plt.imshow(array, origin='lower',  cmap='binary_r', alpha=0.2)
        plt.title(str(kic), fontsize=20)

        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)

        plot = plt.grid(axis = 'both', color='grey', ls = ':', linewidth=6)

        plt.savefig('plot2.jpg')
        # Remove datasets
        os.remove("kic_fits.fits")

    return ruwe_list, primary_g_mag, g_mag_list, mag_diff_list, flux_ratio_list, percentage_flux, flux_contamination_total
