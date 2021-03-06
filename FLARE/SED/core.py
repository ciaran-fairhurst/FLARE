
import numpy as np

from . import IGM

from ..photom import *


def default_cosmo():

    from astropy.cosmology import WMAP9 as cosmo

    return cosmo







class sed():

    def __init__(self, lam, description = False):
    
        self.description = description

        self.lam = lam # \AA
        self.lnu = np.zeros(self.lam.shape) # luminosity ers/s/Hz
        
    def get_l(self): # luminosity  erg/s
      
        nu = physics.constants.c / (self.lam * 1E-6)
        
        return self.Lnu * nu
         
         
    def get_Lnu(self, F): # broad band luminosity/erg/s/Hz
      
        self.Lnu = {f: np.trapz(self.lnu * F[f].T, self.lam) / np.trapz(F[f].T, self.lam) for f in F['filters']}
         
    def get_fnu(self, cosmo, z, include_IGM = True): # flux nJy, depends on redshift and cosmology 

        self.lamz = self.lam * (1. + z)

        self.fnu = 1E23 * 1E9 * self.lnu * (1.+z) / (4 * np.pi * cosmo.luminosity_distance(z).to('cm').value**2) # nJy
        
        if include_IGM:
        
            self.fnu *= IGM.madau(self.lamz, z)
        
    def get_Fnu(self, F): # broad band flux/nJy
                        
        self.Fnu = {f: np.trapz(self.fnu * F[f].T, self.lamz) / np.trapz(F[f].T, self.lamz) for f in F['filters']}
             

        

