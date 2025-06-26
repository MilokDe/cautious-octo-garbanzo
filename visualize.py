#made my Delvin Demke for the VOSS2025 

import numpy as np
from reproject import reproject_interp
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
import numpy as np
import pandas as pd
import csv



from astropy.wcs.utils import pixel_to_skycoord
import astropy.units as u
#wcs
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales

# Use the Qt5Agg backend for interactive window (outside Jupyter)
matplotlib.use('Qt5Agg')

#include here your images
#in our case, we use 480m as the reference image
image360m = fits.open(get_pkg_data_filename('Data/jw03222-o001_t002_nircam_clear-f360m_i2d.fits'))[1]
image480m = fits.open(get_pkg_data_filename('Data/jw03222-o001_t002_nircam_clear-f480m_i2d.fits'))[1] #project to this
image200w = fits.open(get_pkg_data_filename('Data/jw03222-o001_t002_nircam_clear-f200w_i2d.fits'))[1] 
image212n = fits.open(get_pkg_data_filename('Data/jw03222-o001_t002_nircam_clear-f212n_i2d.fits'))[1]

#Define file where to save the clicked points
csv_file = 'targets_list.csv'
#change orientation of the images
wcs_orig = WCS(image480m.header)
# Determine pixel scale (arcsec/pixel)
pixel_scales = proj_plane_pixel_scales(wcs_orig)  # in degrees/pixel
# Center of original image in sky coords
ny, nx = image480m.data.shape
center = wcs_orig.pixel_to_world(nx / 2, ny / 2)
# Define new WCS with North up, RA on x-axis
new_wcs_480m = WCS(naxis=2)
new_wcs_480m.wcs.crpix = [nx / 2, ny / 2]  # center in pixel coords
new_wcs_480m.wcs.cdelt = [-pixel_scales[0], pixel_scales[1]]  # RA left to right (negative), Dec bottom to top
new_wcs_480m.wcs.crval = [center.ra.deg, center.dec.deg]  # center in sky coords
new_wcs_480m.wcs.ctype = ["RA---TAN", "DEC--TAN"]
rotated_data480m, footprint = reproject_interp((image480m.data, wcs_orig), new_wcs_480m, shape_out=image480m.data.shape)


#reproject other images to the new WCS
reprojected_200w, footprint = reproject_interp(image200w, new_wcs_480m, shape_out=image480m.data.shape)
reprojected_212n, footprint = reproject_interp(image212n, new_wcs_480m, shape_out=image480m.data.shape)
reprojected_360m, footprint = reproject_interp(image360m, new_wcs_480m, shape_out=image480m.data.shape)


catalog = pd.read_csv('3222_source-catalog_VOSS2025.csv')
primary_catalog = catalog[catalog['Area'] == 1]
filler_catalog = catalog[catalog['Area'] == -1]
Ta_candidate_catalog = catalog[catalog['Reference'] == True]

# convert to WCS coordinates
ra_prim = catalog['RA'].values
dec_prim = catalog['DEC'].values
x_prim, y_prim = new_wcs_480m.world_to_pixel_values(ra_prim, dec_prim)
# convert to WCS coordinates for filler catalog
ra_fill = filler_catalog['RA'].values
dec_fill = filler_catalog['DEC'].values
x_fill, y_fill = new_wcs_480m.world_to_pixel_values(ra_fill, dec_fill)
# convert to WCS coordinates for Ta candidate catalog
ra_Ta = Ta_candidate_catalog['RA'].values
dec_Ta = Ta_candidate_catalog['DEC'].values
x_Ta, y_Ta = new_wcs_480m.world_to_pixel_values(ra_Ta, dec_Ta)

#include the core coordinates as a refrence point
core_central_coords = SkyCoord('16h32m30.51s', '−24d28m53.7s', frame='icrs')
x_core, y_core = new_wcs_480m.world_to_pixel_values(core_central_coords.ra.deg, core_central_coords.dec.deg)



# Select filter
filters = [reprojected_200w, reprojected_212n, reprojected_360m, rotated_data480m]
filter_name = ['F200W', 'F212N', 'F360M', 'F480M']
current_filter = 3 # Default to F480M because it is the reference image(and looks the best)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=new_wcs_480m)

im = ax.imshow(filters[current_filter], cmap='viridis', vmin=0, vmax=2, origin='lower')
ax.scatter(x_Ta, y_Ta, s=20, facecolors='none', edgecolors='black', label='Catalog Sources', alpha=1, marker='o')
ax.scatter(x_prim, y_prim, s=10, facecolors='none', edgecolors='blue', label='Primary Catalog', alpha=1, marker='o')
ax.scatter(x_fill, y_fill, s=10, facecolors='none', edgecolors='red', label='Filler Catalog', alpha=1, marker='o')
ax.scatter(x_core, y_core, s=50, c='r', label='Core', alpha=1, marker='x')

plt.colorbar(im, ax=ax)
ax.legend(loc='upper right')
ax.set_title(f'Catalog Sources on {filter_name[current_filter]} Image')
ax.set_xlabel('RA (pixels)')
ax.set_ylabel('Dec (pixels)')
ax.grid(True)
plt.tight_layout()

# Hover annotation setup (could be useful, is disabled right now)
annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)


def update_annot(event):
    if event.inaxes == ax and event.xdata is not None and event.ydata is not None:
        x, y = event.xdata, event.ydata
        skycoord = pixel_to_skycoord(x, y, new_wcs_480m)
        ra_str = skycoord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
        dec_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4, alwayssign=True, pad=True)
        annot.xy = (x, y)
        annot.set_text(f"RA = {ra_str}\nDec = {dec_str}")
        annot.set_visible(True)
        fig.canvas.draw_idle()

#if wished, add here more buttons 
def key_press_event(event):
    if event.key == 't':
        if event.inaxes == ax and event.xdata is not None and event.ydata is not None:
            #mark the clicked point
            x, y = event.xdata, event.ydata
            ax.scatter(x, y, s=50, c='yellow', label='Clicked Point', alpha=1, marker='*')
            skycoord = pixel_to_skycoord(x, y, new_wcs_480m)
            ra_str = skycoord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
            dec_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4, alwayssign=True, pad=True)
    if event.key == 'o':
        if event.inaxes == ax and event.xdata is not None and event.ydata is not None:
            #mark the clicked point
            x, y = event.xdata, event.ydata
            ax.scatter(x, y, s=50, c='yellow', label='Clicked Point', alpha=1, marker='*')
            skycoord = pixel_to_skycoord(x, y, new_wcs_480m)
            ra_str = skycoord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
            dec_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4, alwayssign=True, pad=True)

#save clicked point to a file. This function can be improved to import more information
def key_release_event(event):
    if event.key == 't':
        x, y = event.xdata, event.ydata
        ax.scatter(x, y, s=50, c='yellow', label='Clicked Point', alpha=1, marker='*')
        skycoord = pixel_to_skycoord(x, y, new_wcs_480m)
        ra_str = skycoord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
        dec_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4, alwayssign=True, pad=True)
        print(f'Got coordinates')
        with open(csv_file, 'r', newline='') as f:
            reader = csv.reader(f)
            header = next(reader)
        new_row = [np.nan] * len(header) #Empty row for new point
        #include things
        new_row[header.index('RA')] = ra_str
        new_row[header.index('DEC')] = dec_str
        #standard values
        new_row[header.index('Redshift')] = 0.0
        new_row[header.index('NRS_F140X')] = -1
        new_row[header.index('Reference')] = False
        new_row[header.index('FWHM')] = 0.0
        new_row[header.index('R50')] = 0.0
        new_row[header.index('Area')] = 1
        new_row[header.index('Ellipticity')] = 0.0
        new_row[header.index('Theta')] = 0.0

        new_row[header.index('Stellarity')] = 1 #star
        #Write the new row
        with open(csv_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(new_row)

    if event.key == 'o':
        x, y = event.xdata, event.ydata
        ax.scatter(x, y, s=50, c='yellow', label='Clicked Point', alpha=1, marker='*')
        skycoord = pixel_to_skycoord(x, y, new_wcs_480m)
        ra_str = skycoord.ra.to_string(unit=u.hourangle, sep=':', precision=4, pad=True)
        dec_str = skycoord.dec.to_string(unit=u.deg, sep=':', precision=4, alwayssign=True, pad=True)

        with open(csv_file, 'r', newline='') as f:
            reader = csv.reader(f)
            header = next(reader)
        new_row = [np.nan] * len(header) #Empty row for new point
        #include things
        new_row[header.index('RA')] = ra_str
        new_row[header.index('DEC')] = dec_str
        #standard values
        new_row[header.index('Redshift')] = 0.0
        new_row[header.index('NRS_F140X')] = -1
        new_row[header.index('Reference')] = False
        new_row[header.index('FWHM')] = 0.0
        new_row[header.index('R50')] = 0.0
        new_row[header.index('Area')] = 1
        new_row[header.index('Ellipticity')] = 0.0
        new_row[header.index('Theta')] = 0.0
        new_row[header.index('Stellarity')] = 0 #not a star (check it)
        print('Clicked point for not a star at RA = {ra_str}, Dec = {dec_str}')
        #Write the new row
        with open(csv_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(new_row)    

# Radio buttons to select filter
rax = plt.axes([0.01, 0.5, 0.1, 0.2])  # [left, bottom, width, height]
radio = RadioButtons(rax, filter_name, active=current_filter)

def update_filter(label):
    global current_filter, im
    idx = filter_name.index(label)
    current_filter = idx
    im.set_data(filters[current_filter])
    ax.set_title(f'Catalog Sources on {filter_name[current_filter]} Image')
    fig.canvas.draw_idle()

radio.on_clicked(update_filter)

#fig.canvas.mpl_connect("motion_notify_event", update_annot) #not activated
fig.canvas.mpl_connect("key_press_event", key_press_event)
fig.canvas.mpl_connect("key_release_event", key_release_event)


# Add info text to the plot
ax.text(0.01, 0.01, 'Press "t" to mark stars\nPress "o" for objects',
        transform=ax.transAxes, fontsize=10, color='white',
        bbox=dict(facecolor='black', alpha=0.5))
plt.show()


