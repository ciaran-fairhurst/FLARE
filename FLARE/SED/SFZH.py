import numpy as np


class empty: pass

## rewritten to rely on an SPS dictionary of the form:
## dict_keys(['lam', 'log10Z', 'log10Z_edges', 'log10age', 'log10age_edges', 'fraction_remaining', 'stellar', 'nebular' ])



def simple(SPS, p):

    # --- returns the SFZH for a single age/metallicity

    SFZH = np.zeros(SPS.grid['log10Q'].shape)

    iZ = (np.abs(SPS.grid['log10Z'] - p['log10Z'])).argmin()

    a = np.abs(SPS.grid['log10age'] - p['log10age'])

    if np.min(a)==0.0:

        ia = a.argmin()

        SFZH[ia, iZ] = 1.0

    else:

        # --- determine bracketing indicies

        ias = a.argsort()[:2]
        SFZH[ias, iZ] = 1./a[ias]

    SFZH /= np.sum(SFZH)

    return SFZH


def constant(SPS, p):


    # --- returns the SFZH for a constant SFH and single metallicity. At the moment this just chooses a single metallicity.



    a0 = 10**p['log10_duration']

    log10ages = SPS.grid['log10age']

    SFZH = np.zeros(SPS.grid['log10Q'].shape)

    # --- determine the metallicity bin in which to place all the SF

    iZ = (np.abs(SPS.grid['log10Z'] - p['log10Z'])).argmin()

    prev = 0.0

    cont = True

    for ia, log10age in enumerate(log10ages[:-1]):

        c = 10**np.mean([log10ages[ia+1],log10age])

        if cont:

            if a0<c:

                w = a0 - prev # --- determine age width of bin
                cont = False

            else:

                w = c - prev # --- determine age width of bin

            SFZH[ia, iZ] = w

        else:

            w = 0

        prev = c

        # print(ia,log10age,c,w)

    SFZH /= np.sum(SFZH)

    SFZH *= 10**p['log10M*']

    SFR = 10**p['log10M*']/10**p['log10_duration']

    return SFZH, SFR




def F19(SPS,p):
    '''
    SPS.grid DICT: dict_keys(['lam', 'log10Z', 'log10Z_edges', 'log10age',
                         'log10age_edges', 'fraction_remaining', 'stellar', 'nebular' ])

    p DICT : dict_keys([ 'log10M*', 'tau_age', 'dlog10Z/dt', 'log10Z0', 'sigma_log10Z'])
    '''

    SFZH_at_edges =  np.exp( -(SPS.grid['log10Z_edges'][None,:] - p['dlog10Z/dt']*(10**SPS.grid['log10age_edges'][:,None]) - p['log10Z0'])**2/(2*p['sigma_log10Z']**2) \
    - 10**(SPS.grid['log10age_edges'][:,None]/p['tau_age'] ) )

    SFZH      =  np.mean(  np.array( [grid_edges[:-1,:-1],grid_edges[1:,:-1],grid_edges[:-1,1:],grid_edges[1:,1:]] ), axis=0 )

    SFHZ =  SFHZ * np.diff(SPS.grid['log10Z_edges'][None,:]) * np.diff(10**10**(SPS.grid['log10age_edges'][:,None])

    return SFHZ / np.sum(SFHZ) * 10**p['log10M*']



def exp(SPS,p):

    '''
    SPS.grid DICT: dict_keys(['lam', 'log10Z', 'log10Z_edges', 'log10age',
                         'log10age_edges', 'fraction_remaining', 'stellar', 'nebular' ])

    p DICT : dict_keys([ 'log10M*', 'tau_age', 'log10Z'])
    '''

    iZ = (np.abs(SPS.grid['log10Z'] - p['log10Z'])).argmin()

    SFZH = np.zeros(SPS.grid['log10Q'].shape)


    _SFH_edges = np.exp( -SPS.grid['log10age_edges']/p['tau_age']   ) #SFH density in 1d at the edges
    _SFH       = np.mean( np.array([_SFH_edges[1:],_SFH_edges[:-1]), axis=0) #SFH density in 1d at the edges


    SFHZ[:,iZ] = _SFH

    return SFHZ/np.sum(SFHZ) * 10**p['log10M*']
