
import numpy as np

geo = 4.*np.pi*(100.*10.*3.0867*10**16)**2 # factor relating the L to M in cm^2



def flux_to_m(flux):

    return -2.5*np.log10(flux/1E9) + 8.9 # -- assumes flux in nJy

def m_to_flux(m):

    return 1E9 * 10**(-0.4*(m - 8.9)) # -- flux returned nJy

def flux_to_L(flux, cosmo, z):

    return flux*(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2)/(1E23 * (1.+z))

def lum_to_flux(lum, cosmo, z):

    return 1E23 * lum * (1.+ z)/(4.*np.pi*cosmo.luminosity_distance(z).to('cm').value**2) 


def lum_to_M(lum):

    return -2.5*np.log10(lum/geo)-48.6