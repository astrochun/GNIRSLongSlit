def gaussian(x, mu, sig):
    # + on 04/07/2017
    return 1./(np.sqrt(2.*pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
#enddef

def gaussian_R(x_arr, lambda_cen, R_spec):
    '''
    Generate an array consisting of a Gaussian profile given the
    spectral resolution

    Parameters
    ----------
    x_arr : array
      An array of wavelengths

    x_lambda : float
      The central wavelength of the Gaussian line

    R_spec : float or double
      Spectral resolution to consider width of emission lines (e.g., R = 3000)

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 4 July 2017
     - Copied from Zcalbase_gal.observing.locate_em_lines
    '''

    t_FWHM = lambda_cen / R_spec # FWHM based on the wavelength of interest
    temp   = gaussian(x_arr, lambda_cen, t_FWHM/(2 * np.sqrt(2*np.log(2))))
    return temp
#enddef
