#
# ------- stemh.py -------------
#
#
#  functions:
#
#   wavelen() return electron wavelength
#   stemhr()  calculate ADF-STEM probe PSF
#   stemhk()  calculate ADF-STEM probe MTF
#   lenr,i()  aberration function to integrate
#
from pylab import *   #  gets sqrt()
from scipy.special import *     # for j0()
from scipy.integrate import *   # for quad

#-------- wavelength() ----------------
#
# electron wavelength for energy kev
def wavelen(kev):
    return( 12.3986/sqrt((2*511.0+kev)*kev) )

#
#-------- stemhr() ----------------
#
# function stemhr to calculate 
#            STEM probe profile vs. r
#    input array r has the radial positions (in Angs.)
#    input variable params has the optical parameters
#    output array contains the psf 
#
#  params = [kev, Cs3, Cs5, df, amax ]
#
#  kev  = electron energy (in keV)
#  Cs3  = third order spherical aberration (in mm)
#  Cs5  = fifth order spherical aberration (in mm)
#  df   = defocus (in Angstroms)
#  amax = objective aperture  (in mrad)
#
def stemhr( r, params ):
    global w2, w4, w6, intr  # constants for lensr,i()
    kev = params[0]
    Cs3 = params[1]*1.0e7
    Cs5 = params[2]*1.0e7
    df = params[3]
    amax = params[4]*0.001
    wav = wavelen(kev)  # electron wavelength
    kmax = amax/wav
    w2 = wav*pi*df
    w4 = 0.5*pi*Cs3*wav*wav*wav
    w6 = pi*Cs5*wav*wav*wav*wav*wav /3.0
    nr = len( r )
    psf = empty(nr) # make array to fill
    for ir in range(0, nr, 1):
        intr = 2*pi*r[ir]
        # use adaptive quadrature because integrand
        #     not well behaved
        hrr= quad( lensr, 0, kmax )
        hri= quad( lensi, 0, kmax )
        psf[ir] = hrr[0]*hrr[0] + hri[0]*hri[0]

    a = max(psf)
    psf = psf/a  # norm. probe intensity to a max. of 1
    return psf

#
#-------- stemhrCc() ----------------
#
# function stemhrCc to calculate 
#       STEM probe profile vs. r including chromatic aber.
#    input array r has the radial positions (in Angs.)
#    input variable params has the optical parameters
#    output array contains the psf 
#
#  params = [kev, Cs3, Cs5, df, amax, ddf ]
#
#  kev  = electron energy (in keV)
#  Cs3  = third order spherical aberration (in mm)
#  Cs5  = fifth order spherical aberration (in mm)
#  df   = defocus (in Angstroms)
#  amax = objective aperture  (in mrad)
#  ddf  = FWHM of defocus spread (in Angstroms)
#
# Gauss-Hermite Quadrature 
#       with exp(-x*x) weighting of integrand 
#       from Abramowitz and Stegun, and Numerical Recipes
#
def stemhrCc( r, params ):
    NGH = 9  # number of Gauss-Hermete coeff. to use
    # absiccas and weights for Gauss-Hermite Quadrature 
    xGH= array([ 3.190993201781528, 2.266580584531843, 1.468553289216668,
        0.723551018752838, 0.000000000000000, -0.723551018752838,
        -1.468553289216668,-2.266580584531843,-3.190993201781528] )
    wGH= array( [3.960697726326e-005, 4.943624275537e-003 ,8.847452739438e-002,
        4.326515590026e-001, 7.202352156061e-001, 4.326515590026e-001,
        8.847452739438e-002, 4.943624275537e-003, 3.960697726326e-005] )
    
    df0 = params[ 3 ]   #  defocus mean
    ddf = params[ 5 ]   # defocus spread FWHM in Angst.

    #  no defocus integration with small df spread
    if( ddf < 1.0 ):
        psf = stemhr( r, params )
        a = max(psf)
        psf = psf/a  # norm. probe intensity to a max. of 1
        return psf

    nr = len( r )
    psf = empty(nr) # make array to fill
    psf = 0

    #  integrate over defocus spread
    ndf = NGH
    ddf2 = sqrt(log(2.0)/(ddf*ddf/4.0))  # convert from FWHM 
    for idf in range(0, ndf, 1):
        df = df0 + xGH[idf]/ddf2
        weight =  wGH[idf]
        #print "df step ", idf, " df= ", df  # debug
        params[3] = df
        psf1 = stemhr( r, params )
        a = sum( psf1 * r )  # normalize to total current
        psf = psf + weight*psf1/a

    params[ 3 ] = df0    #  restore original value
    a = max(psf)
    psf = psf/a  # norm. probe intensity to a max. of 1
    return psf

#
#-------- prbsize() ----------------
#
# function prbsize() to calculate 
#       FWHM-II size from results stemhr()
#       II = integrated intensity meaning diameter of half current
#    input array r has the radial positions (in Angs.)
#    input array psf has probe intensity from stemhr()
#    output is the size
#
def prbsize( r, psf ):
    nr = len( r )
    psf2 = copy(psf)  # make new array to work with
    psf2[0] = psf2[0] * r[0];
    for ir in range(1, nr, 1):  # integrated intensity
        psf2[ir] = psf2[ir-1] + psf2[ir]*r[ir]

    amax = 0.5 * psf2[nr-1]  # should be half max. value
    for ir in range(nr-1, 0, -1):
        j = ir
        if( psf2[ir] < amax ) : break

    if( j <= 1 ): return 2.0*r[1]  # avoid divide by zero below

    #  interpolate to get a little more accurate
    d = abs( (r[j+1]-r[j])*(amax-psf2[j])/(psf2[j+1]-psf2[j]) )
    return 2.0*( r[j] + d);


#
#-------- stemhk() ----------------
#
#   function to calculate STEM mtf vs. k
#    input array k has the spatial freq. (in inv. Angs.)
#    input variable params has the optical parameters
#          [Cs, df, kev, amax] as elements
#    output array contains the transfer function
#
#  params = [kev, Cs3, Cs5, df, amax, ddf ]
#
#  kev  = electron energy (in keV)
#  Cs3  = third order spherical aberration (in mm)
#  Cs5  = fifth order spherical aberration (in mm)
#  df   = defocus (in Angstroms)
#  amax = objective aperture  (in mrad)
#  ddf  = FWHM of defocus spread (in Angstroms)
#
def stemhk( k, params ):
    kev = params[0]
    Cs3 = params[1]*1.0e7
    #  first calculate the psf using stemhr()
    nr = 500  # number of points in integral over r
    wav = wavelen(kev)  # electron wavelength
    Cs = abs(Cs3)
    if Cs < 0.1e7 : Cs = 0.1e7
    rmax = 2.0*sqrt( sqrt( Cs*wav*wav*wav) )
    r = linspace(0, rmax, nr)
    psf = stemhrCc( r, params )
    # next invert psf to get mtf
    nk = len( k )
    mtf = empty( nk )
    for ik in range(0, nk, 1 ):
        h = psf * j0(2*pi*r*k[ik] ) *r
        mtf[ik] = sum(h)

    a = mtf[0]
    mtf = mtf/a  # normalize to mtf(0)=1
    return mtf

#-------- lens(k) ----------------
#
#  dummy function to integrate (used by stempsf)
#  to calculate complex aberr. function
#    input k (in 1/Angs.), wav = electron wavelength
#
#  chi = pi*wav*k^2*[ 0.5*Cs3*wav^2*k^2 
#                  + (1/3)*Cs5*wav^4*k^4 - df ]
#  return exp( -i*chi )
#
# globals:
#   w2 = pi*defocus*wav
#   w4 = 0.5*pi*Cs3*wav^3
#   w6 = (1.0/3.0)*pi*Cs5*wav^5
#   intr = 2*pi*r
#
#  started 25-jun-2011 E. Kirkland
#
def lensr( k ):    #  real part
    global w2, w4, w6, intr  # constants
    k2 = k*k
    w = ( (w6*k2 + w4) *k2 - w2 )*k2
    return cos(w) * j0( intr*k )*k

def lensi( k ):    #  imag part
    global w2, w4, w6, intr  # constants
    k2 = k*k
    w = ( (w6*k2 + w4) *k2 - w2 )*k2
    return -sin(w) * j0( intr*k )*k

