import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.stats import mad_std
from matplotlib.colors import LogNorm
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture, CircularAnnulus

RR = '43'   # the number of the RR Lyrae in question in my list (01-94)
RA = 149.08947254547863   # RA and dec from Gaia ID
dec = 33.84335914950056
epoch = str(input('Enter the epoch of observation as either "xc##" (include zeros) or a two digit number, not including zeros: '))
if "xc" in epoch:
    filename = '/u/m/amw58/dos/MPHYS_PROJECT/spitzer_images/'+RR+'/'+RR+'__e1_'+epoch+'/'+RR+'__e1_'+epoch+'_3p6um.fits'
else:
    filename = '/u/m/amw58/dos/MPHYS_PROJECT/spitzer_images/'+RR+'/'+RR+'__e'+epoch+'/'+RR+'__e'+epoch+'_3p6um.fits'

RR43 = fits.open(filename)
image_data = RR43[0].data
image_header = RR43[0].header

th = image_header['CROTA2']*np.pi/180   # angle of Y axis, west of north, in radians
d_RA = RA - image_header['CRVAL1']   # differences between known RA and dec of star and of center pixel
d_dec = dec - image_header['CRVAL2']
dX = -(d_RA*np.cos(th) + d_dec*np.sin(th))*6000   # calculated X and Y distances from center to star
dY = -(d_RA*np.sin(th) - d_dec*np.cos(th))*6000   
Xstar = int(np.floor(dX + image_header['CRPIX1']))   # calculated X and Y positions of star, rounded
Ystar = int(np.floor(dY + image_header['CRPIX2']))
region = image_data[Ystar-100:Ystar+100,Xstar-100:Xstar+100]   # defines the region of interest based on guessed star position

bkg_sigma = mad_std(region,ignore_nan=True)   # the stdev of background noise, to limit detection to real stars
detection_w = 5.   # the width threshold of stars which can be detected
aperture_w = 10.   # the aperture width to be used in aperture photometry

daofind = DAOStarFinder(fwhm=detection_w,threshold=50.*bkg_sigma)   # defines a star finding command
sources = daofind(region)   # applies this command to the region of interest
for col in sources.colnames:   # does something with formatting, admittedly I can't quite remember why
    sources[col].info.format='%.8g'
    
#the aperture photometry bit. Defines an aperture with radius (not width, oops, I'll sort that) of star_w, and
#an annulus around the star to pick up background noise. 
    
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=aperture_w/2)
annuli = CircularAnnulus(positions,r_in=10.,r_out=12.)
totals = [apertures,annuli]
phot_table = aperture_photometry(region, totals)
for col in phot_table.colnames:
    phot_table[col].info.format='%.8g'
    
#calculates the mean background from the annulus, then corrects the
#aperture brightness for this

bkg_mean = (phot_table['aperture_sum_1'])/annuli.area()
bkg_sum = bkg_mean*apertures.area()
final_sum = phot_table['aperture_sum_0']-bkg_sum
phot_table['final_sum']=final_sum
print phot_table['final_sum']

#plots it so I know what's going on (plus it's pretty, inferno is a great colourmap)
    
plt.figure(figsize=(10,10),dpi=80)
plt.imshow(region,cmap='inferno',origin='left',aspect='equal',norm=LogNorm())
apertures.plot(color='green',lw=1.5,alpha=1)
annuli.plot(color='white',lw=.5,alpha=1)
plt.show()