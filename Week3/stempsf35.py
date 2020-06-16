#  stempsf.py
#
#
from numpy import *
from pylab import *
from stemh import *

print ('Plot STEM probe intensity')

#  remember; in 3.5 input alwats returns a str (no numbers)
kev  = float64( input( 'Type electron energy in keV : ') )
Cs3  = float64( input( 'Type spherical aberration Cs3 in mm : ') )
Cs5  = float64( input( 'Type spherical aberration Cs5 in mm : ') )
df   = float64( input( 'Type defocus df in Angstroms : ') )
amax = float64( input( 'Type obj. apert. semiangle in mrad : ') )
ddf  = float64( input( 'Type defocus spread FWHM in Angstroms : ') )
p = [kev, Cs3, Cs5, df, amax, ddf ]
#
wav = wavelen(kev)  # electron wavelength
Cs = abs(Cs3)
if Cs < 0.1 : Cs = 0.1
rmax = sqrt( sqrt( Cs*1.0e7*wav*wav*wav ))
npts = 300;   #  number of points in curve
r = linspace( 0, rmax, npts)
psf = stemhrCc( r, p )
rfwhm = prbsize( r, psf )
print ("FWHM-II= ", rfwhm, " Ang.")
plot( r, psf )
xlabel( 'radius in Angstroms');
ylabel( 'PSF' );
s1 = 'Cs3= %gmm, Cs5= %gmm, ' % (Cs3, Cs5 )
s2 = 'df= %gA, ' % df
s3 = 'E= %gkeV, OA= %gmrad, ' % (kev, amax)
s4 = 'ddf= %gA' % (ddf)
title( s1+s2+s3+s4);
#savefig( 'stempsf.eps' )  # select file format
savefig( 'stempsf.pdf' )
show()
