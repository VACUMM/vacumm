# -*- coding: utf8 -*-
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

import numpy as N
from vacumm.misc import kwfilter
from vacumm.misc.grid import resol
import matplotlib.pyplot as P
import cdms2

def get_kxky(shape, dx, dy, verbose=False):
    """Generate 2D wavenumbers"""
    # to work with dimensional axes
    if hasattr(shape, 'shape'):
        shape = shape.shape
    nlat, nlon = shape
    Lx = nlon*dx
    Ly = nlat*dy
    # to work with non-dimensional axes
    # Lx = N.pi
    # Ly = Lx * nlat / nlon * dy / dx
    coefx = Lx / (2*N.pi)
    coefy = Ly / (2*N.pi)
    if verbose: print "Lx : %s Ly : %s coefx : %s coefy : %s  \n" %(Lx,Ly,coefx,coefy)
    kx = N.hstack( ( N.arange(0.,nlon/2), N.array([0]), -N.arange(1.,nlon/2)[::-1] ) ) / coefx
    ky = N.hstack( ( N.arange(0.,nlat/2), N.array([0]), -N.arange(1.,nlat/2)[::-1] ) ) / coefy
    [kkx,kky] = N.meshgrid(kx,ky)
    kk = N.sqrt( kkx**2 + kky**2 )
    return kx,ky,kkx,kky,kk,Lx,Ly


def get_spec(var1, var2=None, dx=None, dy=None, verbose=False, fft=False, **kwargs):
    """Get the spectrum of a 2D variable

    :Return: specvar,nbwave,dk
    """
    # Get the resolution in meters
    kwresol = kwfilter(kwargs, 'resol_')
    if None in [dx, dy]:
        if not cdms2.isVariable(var1):
            raise TypeError('var1 must be a MV2 array to compute dx and dy')
        kwresol.setdefault('proj', True)
        ldx, ldy = resol(var1, **kwresol)
        dx = dx or ldx
        dy = dy or ldy
    if N.ma.isMA(var1): var1= var1.filled(0.)
    if N.ma.isMA(var2): var2 = var2.filled(0.)


    # Part 1 - estimate of the wavenumbers
    [kx,ky,kkx,kky,kk,Lx,Ly] = get_kxky(var1, dx, dy, verbose=verbose)
    if verbose:
        print "dx = %s, fy = %s " %(dx,dy)
        print "kx = ",kx[0:3]
        print "ky = ",ky[0:3]
        print "kkx[0:3,2] = ",kkx[0:3,2]
        print "kky[0:3,1] = ",kky[0:3,1]
        print "kk[0:3,3] = ",kk[0:3,3]
        print "shape",kx.shape,ky.shape,kkx.shape,kky.shape,kk.shape

    # Part 2 - estimate of the spectrum
    # - fast fourier transform
    if fft:
        hat_phi1=N.fft.fft2(var1)
        hat_phi1=hat_phi1*hat_phi1.conj()
        hat_phi1=hat_phi1.real.copy()  # useless
    else:
        hat_phi1 = var1

    if var2 is not None:
        if fft:
            hat_phi2=N.fft.fft2(var2)
            hat_phi2=hat_phi2*hat_phi2.conj()
            hat_phi2=hat_phi2.real.copy()  # useless
        else:
            hat_phi2 = var2
        hat_spec = hat_phi1 + hat_phi2
    else:
        hat_spec = hat_phi1
    # - integration of the spectrum
    #   shift to have values centered on the intervals
    dk = kx[1] - kx[0]
    dk_half = dk/2  # half of the interval
    k_= kx[0:min(var1.shape[1],var1.shape[0])/2] + dk
    specvar = k_*0
    for i in range(len(k_)):
        # get indexes that satisfy two conditions
        # integration over a spectral width
        specvar[i] = ( hat_spec[ (kk<=k_[i]+dk_half) & (kk>k_[i]-dk_half) ] ).sum()

    # Normalisation (Danioux 2011)
    specvar /= (var1.shape[1]*var1.shape[0])**2 * dk

    # Dimensionalization to be able to compare with the litterature
    nbwave = k_*1000.0/2.0/N.pi     # from rad/m to km/m
    dk *= 1000.0/2.0/N.pi
    specvar *= 2.0*N.pi/1000.0
    if verbose:
        if var2 is not None:
            print "\n Normalized co-spectrum : \n",specvar
        else:
            print "\n Normalized spectrum : \n",specvar
    return specvar,nbwave,dk

def get_cospec(var1, var2, dx=None, dy=None, verbose=False, **kwargs):
    """Co-sprectum of two variables"""
    return get_spec(var1, var2, dx=dx, dy=dy, verbose=verbose, **kwargs)

def save_spec(filename,var,absc,time):
    if os.path.isfile(filename):
        #tmpabsc=N.column_stack((N.load(filename)['nbwave'],absc))
        tmpvar=N.column_stack((N.load(filename)['spectrum'],var))
        tmptime=N.append(N.load(filename)['time'],time)
    else:
        tmpvar = var
        tmptime = time
    N.savez(filename,spectrum=tmpvar,nbwave=absc,time=tmptime)


def plot_loglog_kspec(data,abscisse,title=None,subtitle=None,xlabel=None,ylabel=None,slope=-3,savefig=None):
    # bound the range of the abscisses
    abs_min = round( abscisse.min() , min([v for v in range(10) if round(abscisse.min(),v) != 0]) )
    abs_min = abs_min * 0.9  # 90 percent of the value
    abs_max = round( abscisse.max() , min([v for v in range(10) if round(abscisse.max(),v) != 0]) )
    abs_max = abs_max * 1.1  # 110 percent of the value
    ord_min = round( data.min() , min([v for v in range(10) if round(data.min(),v) != 0]) )
    ord_min = ord_min * 0.9  # 90 percent of the value
    ord_max = round( data.max() , min([v for v in range(10) if round(data.max(),v) != 0]) )
    ord_max = ord_max * 1.1  # 110 percent of the value
    # straight line with a slope of k**slope
    a=data[len(abscisse)/2] / abscisse[len(abscisse)/2]**slope
    line=a*abscisse**slope
    # figure
    P.figure()
    P.loglog(abscisse, data, basex=10)
    P.loglog(abscisse[1:-1], line[1:-1], basex=10)
    text='k'+str(slope)
    P.text(abscisse[len(abscisse)/10], data[len(abscisse)/10]*5,text,color='g')
    P.xlim(abs_min,abs_max)
    P.grid(True,which="major",ls="-")
    P.grid(True,which="minor")
    P.suptitle(title)
    P.title(subtitle)
    P.xlabel(xlabel)
    P.ylabel(ylabel)
    P.tight_layout()
    P.subplots_adjust(top=.9)
    P.show()
    if savefig: P.savefig(savefig+'.png')
    P.close()


def plot_semilogx_kspec(data,abscisse,title=None,subtitle=None,xlabel=None,ylabel=None,savefig=None):
    # bound the range of the abscisses
    abs_min = round( abscisse.min() , min([v for v in range(10) if round(abscisse.min(),v) != 0]) )
    abs_min = abs_min * 0.9  # 90 percent of the value
    abs_max = round( abscisse.max() , min([v for v in range(10) if round(abscisse.max(),v) != 0]) )
    abs_max = abs_max * 1.1  # 110 percent of the value
    ord_min = round( data.min() , min([v for v in range(10) if round(data.min(),v) != 0]) )
    ord_min = ord_min * 0.9  # 90 percent of the value
    ord_max = round( data.max() , min([v for v in range(10) if round(data.max(),v) != 0]) )
    ord_max = ord_max * 1.1  # 110 percent of the value
    # figure
    P.figure()
    P.semilogx(abscisse, data, basex=10)
    P.axhline(linewidth=3, color='r')
    P.xlim(abs_min,abs_max)
    P.grid(True,which="major",ls="-")
    P.grid(True,which="minor")
    P.suptitle(title)
    P.title(subtitle)
    P.xlabel(xlabel)
    P.ylabel(ylabel)
    P.tight_layout()
    P.subplots_adjust(top=.9)
    P.show()
    if savefig: P.savefig(savefig+'.png')
    P.close()


def energy_spectrum(var, dx=None, dy=None, ctime=None, dispfig=False, latmean=None, verbose=False):
    """Compute the energy spectrum of a 2D variable"""
    # Guess dx and/or dy
    if None in [dx, dy]:
        if not cdms2.isVariable(var):
            raise TypeError('var must be a MV2 array to compute dx and dy')
        kwresol.setdefault('proj', True)
        ldx, ldy = resol(var, **kwresol)
        dx = dx or ldx
        dy = dy or ldy

    if latmean is None:
        if not cdms2.isVariable(var):
            raise Exception('You must provide explicitly latmean')
        latmean = var.getLatitude().getValue().mean()

    # Numpy arrays
    sshbox = var.filled(0.) if N.ma.isMA(var) else var

    #####################################################
    # estimate of the spatial resolution and the latitude
    #####################################################

    grav = 9.81
    #latmean = int(float(latmean))
    #dx = int(grid[0])*np.cos(np.pi*latmean/180.)
    f0 = 2*N.pi/(24.*3600.)*2.0*N.sin(N.pi*latmean/180.)
    gof0 = grav/f0
    if verbose:
        print "the spatial resolution of the grid is constant \n dx = %s m, dy = %s m" %(dx,dy)
        print "the averaged latitude is %s degrees\n" %(latmean)


    ###############################################
    #  1ST PART ANALYZE FROM THE SEA SURFACE HEIGHT
    #  --------------------------------------------
    ###############################################

    # zonal average and its removing from ssh
    # mean on i (we erase 1 dimension)
    # and duplication along i (depends on sshbox.mean(1).shape)
    # then we need to transpose to have constant value along i and variations along j
    #revoir sshbox = sshbox - N.tile( sshbox.mean(1),(sshbox.shape[1],1) ).T
    # plot of ssh (box domain)
    #if options.dispfig: plot.contour_xy(sshbox,sshbox.shape[0],sshbox.shape[1],'ssh_box_minus_iaverage')

    # estimate of the rms to further check
    sshbox_rms = N.sqrt( (sshbox**2).mean() )
    if verbose: print "rms of the ssh over the region of interest (from physical field) : %s \n" %sshbox_rms

    ################################
    # create a double periodic field
    ################################

    ssh2per = N.tile( sshbox, (2,2) )
    ssh2per[ssh2per.shape[0]/2:ssh2per.shape[0],0:ssh2per.shape[1]/2] = sshbox[::-1,::]
    ssh2per[::,ssh2per.shape[1]/2:] = N.fliplr(ssh2per[::,0:ssh2per.shape[1]/2])
    #if options.dispfig: plot.contour_xy(ssh2per,ssh2per.shape[0],ssh2per.shape[1],'ssh_double_periodic')
    del sshbox

    # estimate of the rms to further check
    ssh2per_rms = N.sqrt( (ssh2per**2).mean() )
    if verbose: print "rms of the double periodic ssh over the region of interest (from physical field) : %s \n" %ssh2per_rms

    ################################
    # wavenumber spectrum of the ssh (double periodic field)
    ################################

    # estimate of the spectrum
    ssh2per_spec, ssh2per_k, dk = get_spec(ssh2per, dx=dx, dy=dy, fft='x2k', verbose=verbose)  # sshbox in physical space

    # make sure the FFT is correct : we must have the same rms either from physics or from wavenumbers
    ssh2per_rmsk = N.sqrt( (ssh2per_spec*dk).sum() )
    if verbose: print "rms of the ssh over the region of interest (from wave number fields) : %s \n" %ssh2per_rmsk

    # figure
    if dispfig: plot_loglog_kspec(ssh2per_spec, ssh2per_k, savefig='ssh2per_spec',
        title='SPECTRUM OF SSH (double periodic)', subtitle=str(ctime),
        xlabel='Wavenumber (cycles/km)', ylabel='SSH [m^2/(cycle/km)]',
        slope=-5)


    ###########################################
    # the ssh pectrum after doubling the domain
    # is much better than the previous one
    ###########################################

    ####################################################
    #  2ND PART ANALYZIS FROM THE GEOSTROPHIC VELOCITIES
    #  -------------------------------------------------
    ####################################################

    ###########################################
    # estimate of the geostrophic velocities (from the spectrum of the periodic field of ssh)
    ###########################################

    # estimate of the wavenumbers
    kx,ky,kkx,kky,kk,Lx,Ly = get_kxky(ssh2per, dx, dy, verbose=verbose)
    # estimate of the velocities from the spectrum of the ssh
    ussh = -gof0*N.fft.ifft2( (kky*N.fft.fft2(ssh2per))*N.complex(0,1) ).real
    vssh =  gof0*N.fft.ifft2( (kkx*N.fft.fft2(ssh2per))*N.complex(0,1) ).real
    #if options.dispfig: plot.contour_xy(ussh,ussh.shape[0],ussh.shape[1],'ussh (from double periodic field spectral method) ')
    #if options.dispfig: plot.contour_xy(vssh,vssh.shape[0],vssh.shape[1],'vssh (from double periodic field and spectral method) ')
    ssh2per_shape = ssh2per.shape
    del ssh2per

    #######################################################
    # wavenumber co-spectrum of the geostrophic velocities (eke)
    #######################################################

    # estimate of the co-spectrum
    uv2per_cospec,eke_k,dk = get_cospec(ussh, vssh, dx=dx, dy=dy, verbose=verbose, fft=True)

    # estimate of the eke to further check
    uv2per_rmsk = N.sqrt( ((uv2per_cospec*dk)/2.).sum() )
    if verbose: print "eke over the region of interest (from physical field) : %s \n" %uv2per_rmsk
    # make sure the FFT is correct : we must have the same rms either from physics or from wavenumbers
    #ssh2per_rmsk = N.sqrt( (ssh2per_spec*dk).mean() )
    #if verbose: print "eke over the region of interest (from physical field) : %s \n" %sshbox_rms
    #if verbose: print "rms of the ssh over the region of interest (from wave number fields) : %s \n" %ssh2per_rmsk
    #if verbose: print "rms difference (wavenumber - physics) : %s \n" %(ssh2per_rmsk-ssh2per_rms)

    # figure
    if dispfig:
        plot_loglog_kspec(uv2per_cospec, eke_k, savefig='eke2per_spec',
            title='SPECTRUM OF EKE', subtitle=str(ctime),
            xlabel='Wavenumber (cycles/km)', ylabel='EKE [(m/s)^2/(cycle/km)]')


    ####################################################
    #  3RD PART ESTIMATE OF THE FLUXES OF ENERGY
    #  -------------------------------------------------
    ####################################################

    ###########################################
    # estimate of the production of eke
    ###########################################

    # geostrophic velocities in the spectral space are non zero in the bottom-left corner
    ug = ussh.copy() ; del ussh
    ug[ssh2per_shape[0]/2:ssh2per_shape[0],::] = 0.
    ug[::,ssh2per_shape[1]/2:] = 0.
    vg = vssh.copy() ; del vssh
    vg[ssh2per_shape[0]/2:ssh2per_shape[0],::] = 0.
    vg[::,ssh2per_shape[1]/2:] = 0.
    #if options.dispfig: plot.contour_xy(ug,ug.shape[0],ug.shape[1],'ug (from the spectrum of the periodic ssh)')
    #if options.dispfig: plot.contour_xy(vg,vg.shape[0],vg.shape[1],'vg (from the spectrum of the periodic ssh)')


    u = N.fft.fft2(ug)
    v = N.fft.fft2(vg)

    ugx = N.fft.ifft2( (u*kkx)*N.complex(0,1) ).real
    ugy = N.fft.ifft2( (u*kky)*1j ).real
    vgx = N.fft.ifft2( (v*kkx)*1j ).real
    vgy = N.fft.ifft2( (v*kky)*1j ).real
    term_u = N.fft.fft2(ug*ugx + vg*ugy)
    term_v = N.fft.fft2(ug*vgx + vg*vgy)
    eke_prod = - ( u.conj()*term_u + v.conj()*term_v ).real
    eke_prod_phys1=N.fft.ifft2(eke_prod).real
    #if options.dispfig: plot.contour_xy(eke_prod_phys1,eke_prod_phys1.shape[0],eke_prod_phys1.shape[1],'eke_prod_phys1')

    term_u = 1j*kkx*N.fft.fft2(ug*ug) + 1j*kky*N.fft.fft2(vg*ug)
    term_v = 1j*kkx*N.fft.fft2(ug*vg) + 1j*kky*N.fft.fft2(vg*vg)
    eke_prod = - ( u.conj()*term_u + v.conj()*term_v ).real
    eke_prod_phys2=N.fft.ifft2(eke_prod).real
    #if options.dispfig: plot.contour_xy(eke_prod_phys2,eke_prod_phys2.shape[0],eke_prod_phys2.shape[1],'eke_prod_phys2')
    #if options.dispfig: plot.contour_xy(eke_prod_phys1-eke_prod_phys2,eke_prod_phys2.shape[0],eke_prod_phys2.shape[1],'diff_eke_prod')
    phi = eke_prod_phys1-eke_prod_phys2

    # estimate of the total energy flux of the momentum equation
    gg_spec, gg_k, dk = get_spec(eke_prod, dx=dx, dy=dy, verbose=verbose, fft=False)  # phi in Fourier space

    # choose kmax1 such that it corresponds to gg_k(kmax1) < 1/(3*dx)*1000 [km/m] (2dx Shanon + take into account diffusion)
    kmax1 = min( N.where(gg_k>= 1/(3.0*dx)*1000.)[0][0] - 1 ,gg_spec.shape[0])
    if verbose: print "Shanon criterium kmax1 = %s, gg_k(kmax1) = %s, gg_k(kmax1+1) = %s, 1/(3.0*dx)*1000. = %s" %(kmax1,gg_k[kmax1],gg_k[kmax1+1],1/(3.0*dx)*1000.)
    # integration from large scales to small ones (shift one large scales)
    tpi = N.array( [sum(gg_spec[ i:kmax1 ]) for i in range(kmax1)] )
    tpi[0]=0.

    # figure
    if dispfig:
        plot_semilogx_kspec(tpi, gg_k[:kmax1], savefig='ssh_spec',
        title='TOTAL ENERGY FLUX', subtitle=str(ctime),
        xlabel='Wavenumber (cycles/km)', ylabel='Sum of EKE from small scales to k [(m/s)^2/(cycle/km)]')


    if verbose: print "The end !"

    del grav,f0,gof0,sshbox_rms,ssh2per_rms,dk,kx,ky,kkx,kky,kk,Lx,Ly,uv2per_rmsk,ug,vg,u,v,term_u,term_v,ugx,ugy,vgx,vgy,eke_prod,eke_prod_phys1,eke_prod_phys2,phi,gg_spec
    gc.collect()

    return ssh2per_spec, ssh2per_k, uv2per_cospec, eke_k, tpi, gg_k[:kmax1]
