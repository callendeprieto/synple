#!/usr/bin/env python3
#-*- coding: utf-8 -*-

"""Python wrapper for synspec 

Calculation of synthetic spectra of stars and convolution with a rotational/Gaussian kernel.
Makes the use of synspec simpler, and retains the main functionalities (when used from
python). The command line interface is even simpler but fairly limited. 

For information on
synspec visit http://nova.astro.umd.edu/Synspec43/synspec.html.

Example
-------

To compute the solar spectrum between 6160 and 6164 angstroms, using a model atmosphere in
the file ksun.mod (provided with the distribution), with the output going into the file
sun.syn

   $synple.py ksun.mod 6160. 6164. 

To force a micro of 1.1 km/s, and convolve the spectrum with a Gaussian kernel with a fwhm 
of 0.1 angstroms

   $synple.py ksun.mod 6160. 6164. 1.1  0.1

To perform the calculations above in python and compare the emergent normalized profiles

   >>> from synple import syn
   >>> s = syn('ksun.mod', (6160.,6164.))
   >>> s2 = syn('ksun.mod', (6160.,6164.), vmicro=1.1, fwhm=0.1)

   in plain python
   >>> import matplotlib.pyplot as plt
   >>> plt.ion()
   >>> plt.plot(s[0],s[1]/s[2], s2[0], s2[1]/s2[2])

   or ipython
   In [1]: %pylab
   In [2]: plot(s[0],s[1]/s[2], s2[0], s2[1]/s2[2])


"""
import os
import re
import sys
import stat
import string
import random
import subprocess
import glob
import time
import copy
import gzip
import yaml
import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from itertools import product
from astropy.io import fits
import astropy.table as tbl
import astropy.units as units
import datetime
import platform

#configuration
#synpledir = /home/callende/synple
synpledir = os.path.dirname(os.path.realpath(__file__))

#relative paths
modeldir = os.path.join(synpledir, "models")
modelatomdir = os.path.join(synpledir , "data")
confdir = os.path.join(synpledir, "config")
griddir = os.path.join(synpledir, "grids")
atlasdir = os.path.join(synpledir, "atlases")
linelistdir = os.path.join(synpledir, "linelists")
linelist0 = ['gfATOc.19','gfMOLsun.20','gfTiO.20','H2O-8.20']
bindir = os.path.join(synpledir,"bin")
synspec = os.path.join(bindir, "synspec54")
rotin = os.path.join(bindir, "rotin")

#internal synspec data files
isdf = ['CIA_H2H2.dat',  'CIA_H2H.dat', 'CIA_H2He.dat', 'CIA_HHe.dat', \
        'irwin_bc.dat',  'tremblay.dat', \
        'tsuji.molec_bc2']

#other stuff
clight = 299792.458 # km/s
epsilon = 0.6 #clv coeff.
bolk = 1.38054e-16  # erg/ K
zero = " 0 "
one =  " 1 "
two =  " 2 "



def syn(modelfile, wrange, dw=None, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0,  intensity=False, lineid=False, tag=False,  \
    clean=True, save=False, synfile=None, \
    lte=None, compute=True, tmpdir=None):

  """Computes a synthetic spectrum

  Interface to the fortran codes synspec/rotin that only requires two mandatory inputs: 
  a model atmosphere (modelfile) and the limits of the spectral range (wrange). The code 
  recognizes Kurucz, MARCS and Phoenix LTE model atmospheres. The sampling of the frequency 
  grid is chosen internally, but can also be set by adding a constant wavelength step (dw).
  The abundances and microturbulence velocity can be set through the abu and vmicro 
  parameters, but default values will be taken from the model atmosphere. Rotational and 
  Gaussian broadening can be introduced (vrot and fwhm parameters), as well as radial-tangential
  macroturbulence. The computed spectrum  can be written to a file (save == True). 


  Parameters
  ----------
  modelfile : str
      file with a model atmosphere
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  dw: float, optional
      wavelength step for the output fluxes
      this will trigger interpolation at the end.
      A negative value will lead to interpolation to a uniform step in ln(lambda)
      (default is None for automatic frequency selection)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  vmicro: float
      microturbulence (km/s)  
      a negative value triggers the use of the APOGEE DR14 equation      
      (default is None to adopt the value in the model atmosphere)  
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  intensity: bool
      set to True to return intensities  at
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
      -- no broadening due to macro, rotation or instrumental/fwhm effects are considered
      (default False)
  lineid: bool
      set to True to add line identifications to the ouput. They will take the form
      of a list with three arrays (wavelengths, lineids and predicted EWs) 
      (default False)  
  tag: bool
      set to True to produce a plot showing the line identifications. It implicitly
      sets lineid to True, overiding whatever value is passed to lineid
      (default False)
  clean: bool
      True by the default, set to False to avoid the removal of the synspec
      temporary files/links (default True)
  save: bool
      set to True to save the computed spectrum to a file (default False)
      the root of the model atmosphere file, with an extension ".syn" will be used
      but see the parameter synfile to change that
  synfile: str
      when save is True, this can be used to set the name of the output file
      (default None)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default None)
  compute: bool
      set to False to skip the actual synspec run, triggering clean=False
      (default True)
  tmpdir: string
      a temporary directory with this name will be created to store
      the temporary synspec input/output files, and the synple log file (usually named
      syn.log) will be named as tmpdir_syn.log. When tmp is None a random string is used
      for this folder
      (default None)

  Returns
  -------
  wave: numpy array of floats
      wavelengths (angstroms)
  flux: numpy array of floats
      flux (H_lambda in ergs/s/cm2/A)
  cont: numpy array of floats
      continuum flux (same units as flux)

  ---- if intensity is True
  inte: 2D numpy array of floats
      intensity (I_lambda in ergs/s/cm2/A/strad) for
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
  continte: 2D numpy array of float
      continuum intensity for the same angles


  ---- if lineid (or tag) is True
  lalilo: list with three arrays
        la: numpy array of floats
            wavelenghts of lines (angstroms)
        li: numpy array of str
            line identifidation
        lo: numpy array of floats
            predicted equivalent width (miliangstroms)

  """

  #basic checks on the line list and model atmosphere
  linelist, modelfile = checksynspec(linelist,modelfile)

  #read model atmosphere
  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)
 
  if vmicro is None:
      vmicro = vmicro2
  elif vmicro < 0.0:
      #Holtzman et al. 2015, AJ 150, 148
      vmicro = 2.478 - 0.325 * logg

  #print('teff,logg,vmicro=',teff,logg,vmicro)

  if abu is None: abu = abu2
  #we take a step of 1/3 of the Gaussian (thermal + micro) FWHM at the lowest T and for an atomic mass of 100
  space = np.mean(wrange) / clight * 2.355 / 3. * np.sqrt(0.1289**2 * np.min(atmos['t']) / 100. + vmicro** 2 / 2.) 

  #check input parameters are valid
  imode = checkinput(wrange, vmicro, linelist)

  #find out additional info for Tlusty models
  if atmostype == 'tlusty':
    madaffile, nonstdfile, nonstd, numpar, datadir, inlte, atommode, atominfo = read_tlusty_extras(modelfile)
    if (inlte == -1): nonstd['IBFAC'] = 1
    if (lte == True): inlte = 0
  else:
    nonstd = None
    inlte = 0
    atommode = None
    atominfo = None


  print(modelfile,'is a',atmostype,' model')
  print('teff,logg,vmicro=',teff,logg,vmicro)
  #print ('abu=',abu)
  #print (len(abu))
  #print ('nd=',nd)
  #print ('linelist=',linelist)
  #print ('wrange=',wrange)

  logfile = 'syn.log'
  if tmpdir is None:
    tmpdir = ''.join(random.choices(string.ascii_lowercase + string.digits, k = 16))
  startdir = os.getcwd()
  logfile = os.path.join(startdir,os.path.split(tmpdir)[-1]) + "_" + logfile
  try:
    if tmpdir != '.': os.mkdir(tmpdir)
  except OSError:
    print( "cannot create tmpdir %s " % (tmpdir) )
  try:
    os.chdir(tmpdir)
  except OSError:
    print("cannot enter tmpdir %s " % (tmpdir) )


  cleanup_fort()
  if os.path.islink('data'): os.unlink('data')
  if os.path.isfile('tas'): os.remove('tas')

  assert (not os.path.isdir('data')), 'A subdirectory *data* exists in this folder, and that prevents the creation of a link to the data directory for synple'

  #link data folder for used-provided model atoms
  dd = ''
  if atmostype == 'tlusty':
    # data dir
    hdd, dd = os.path.split(datadir)
    os.symlink(datadir,dd)
    if dd == 'data':
      for entry in isdf:
        assert (os.path.isfile(os.path.join(dd,entry))), 'Cannot find the data file:'+dd+'/'+entry     
    #dd2 = ''
    #hdd2, dd2 = os.path.split(madaffile)
    #if hdd2 == '': hdd2 = '..'
    #os.symlink(os.path.join(hdd2,madaffile),'./fort.5')      
  #else:

  write5(teff,logg,abu,atom,inlte=inlte,
           atommode=atommode,atominfo=atominfo) #abundance/opacity file

  write8(teff,logg,nd,atmos,atmostype)          #model atmosphere


  #link the data folder (synspec data + default tlusty model atoms)
  if dd != 'data': os.symlink(modelatomdir,'./data')          #data directory  


  #set iprin to 2 to ouput lineids from synspec
  if tag: lineid = True
  iprin = 0
  if lineid: iprin = 2

  cutoff0=250.
  if logg > 3.0:
    if teff < 4000.: cutoff0=500.
    if teff < 3500.: cutoff0=1000.
    if teff < 3000.: cutoff0=1500. 
  write55(wrange,dw=space,imode=imode,iprin=iprin, inlte=inlte, hydprf=2, \
          cutoff0=cutoff0,strength=strength,vmicro=vmicro,   \
          linelist=linelist,atmostype=atmostype,intensity=intensity)
  #synspec control file
  writetas('tas',nd,linelist,nonstd=nonstd)               #non-std param. file
  create_links(linelist)                                  #auxiliary data

  if compute == False:

    wave = None
    flux = None  
    cont = None
    if intensity: 
      inte = None
      continte = None

  else:

    synin = open('fort.5')
    synout = open(logfile,'w')

    start = time.time()
    p = subprocess.Popen([synspec], stdin=synin, stdout = synout, stderr= synout, shell=True)
    p.wait()

    synout.flush()
    synout.close()
    synin.close()

    assert (os.path.isfile('fort.7')), 'Error: I cannot read the file *fort.7* in '+tmpdir+' -- looks like synspec has crashed, please look at '+logfile

    assert (os.path.isfile('fort.17')), 'Error: I cannot read the file *fort.17* in '+tmpdir+' -- looks like synspec has crashed, please look at '+logfile


    wave, flux = np.loadtxt('fort.7', unpack=True)
    if np.any(np.diff(wave) <= 0.0):
      wave, win = np.unique(wave,return_index=True)
      flux = flux[win] 
    wave2, flux2 = np.loadtxt('fort.17', unpack=True)
    if np.any(np.diff(wave2) <= 0.0):
      wave2,win = np.unique(wave2,return_index=True)
      flux2 = flux2[win]
    if intensity:
      nmu = 10
      assert (os.path.isfile('fort.10')), 'Error: I cannot read the file *fort.10* in '+tmpdir
      iwave, inte = read10('fort.10')
      contiwave, continte1 = read10('fort.18')
      if np.any(np.diff(iwave) <= 0.0):
        iwave, win = np.unique(iwave,return_index=True)
        inte = inte[win,:]

      clight_cgs = clight * 1e5 # km to cm
      iwave_cgs = iwave * 1e-8  # AA to cm
      # I_nu to I_lambda (cgs)
      inte = (inte.transpose() * clight_cgs / iwave_cgs**2 * 1e-8).transpose()

      continte = np.zeros((len(iwave),nmu))
      #mapping continte onto the same wavelength grid as the actual intensity
      for entry in range(nmu):
        continte[:,entry] = np.interp(iwave,contiwave,continte1[:,entry]) 
      # I_nu to I_lambda
      continte = (continte.transpose() * clight_cgs / iwave_cgs**2 * 1e-8).transpose()

      assert (np.max(iwave-wave) < 1e-7), 'Error: the wavelengths of the intensity (fort.10) and flux arrays  (fort.7) are not the same'
      assert (fwhm < 1e-7), 'Error: computing the intensity at various angles is not compatible with the fwhm keyword'
      assert (vmacro < 1e-7), 'Error: computing the intensity at various angles is not compatible with the vmacro keyword'
      assert (vrot < 1e-7), 'Error: computing the intensity at various angles is not compatible with the vrot keyword'

    if dw == None and fwhm <= 0. and vrot <= 0.: cont = np.interp(wave, wave2, flux2)
    end = time.time()
    print('syn ellapsed time ',end - start, 'seconds')

    if fwhm > 0. or vrot > 0. or vmacro > 0.:
      start = time.time()
      print( vrot, fwhm, vmacro, space, steprot, stepfwhm)
      wave, flux = call_rotin (wave, flux, vrot, fwhm, vmacro, \
        space, steprot, stepfwhm, clean=False, \
        reuseinputfiles=True, logfile=logfile)
      if dw == None: cont = np.interp(wave, wave2, flux2)
      end = time.time()
      print('convol ellapsed time ',end - start, 'seconds')

    if (dw != None): 
      if (dw < 0.):
        ldw=np.abs(dw)/np.mean(wrange)
        nsamples = int((np.log(wrange[1]) - np.log(wrange[0]))/ldw) + 1
        wave3 = np.exp(np.arange(nsamples)*ldw + np.log(wrange[0]))
      else:
        nsamples = int((wrange[1] - wrange[0])/dw) + 1
        wave3 = np.arange(nsamples)*dw + wrange[0]
      cont = np.interp(wave3, wave2, flux2)
      flux = np.interp(wave3, wave, flux)
      if intensity:
        nmu = 10
        inte2 = np.zeros((nsamples,nmu))
        continte2 = np.zeros((nsamples,nmu))
        for entry in range(nmu): 
          inte2[:,entry] = np.interp(wave3, wave, inte[:,entry])
          continte2[:,entry] = np.interp(wave3, wave, continte[:,entry])
        inte = inte2
        continte = continte2
      #flux = interp_spl(wave3, wave, flux)      
      wave = wave3

    if lineid == True:
      assert (os.path.isfile('fort.12')), 'Error: I cannot read the file *fort.12* in '+tmpdir+' -- looks like synspec has crashed, please look at '+logfile

      d = np.loadtxt('fort.12', usecols=(0,1,2,3,4,5,6), dtype=str)
      la = []
      li = []
      lo = []
      for i in range(len(d[:,0])): 
        la.append(d[i,0])
        li.append(d[i,1]+d[i,2])
        lo.append(d[i,6])
      la = np.array(la, dtype=float)
      li = np.array(li, dtype=str)
      lo = np.array(lo, dtype=float)
      if (os.path.isfile('fort.15')): 
        d = np.loadtxt('fort.15', usecols=(0,1,2,3,4,5), dtype=str)
        la2 = []
        li2 = []
        lo2 = []
        for i in range(len(d[:,0])): 
          la2.append(d[i,0])
          li2.append(d[i,1])
          lo2.append(d[i,5])
        la2 = np.array(la2, dtype=float)
        li2 = np.array(li2, dtype=str)
        lo2 = np.array(lo2, dtype=float)
        la = np.concatenate((la, la2))
        li = np.concatenate((li, li2))
        lo = np.concatenate((lo, lo2))

    if clean == True: 
      cleanup_fort()
      if os.path.islink('data'): os.unlink('data')
      if os.path.isfile('tas'): os.remove('tas')
      if atmostype == 'tlusty':
        if os.path.islink(dd): os.unlink(dd)


    try:
      os.chdir(startdir)
    except OSError:
      print("cannot change directory from tmpdir %s to startdir %s"  % (tmpdir,startdir) ) 
    if clean == True:
      try:
        os.rmdir(tmpdir)
      except OSError:
        print("cannot remove directory tmpdir %s" % (tmpdir) )
     

    if save == True:

      out = ['MODEL   = '+modelfile+'\n']
      out.append('TEFF    = '+str(teff)+'\n')
      out.append('LOGG    = '+str(logg)+'\n')
      out.append('VMICRO  = '+str(vmicro)+'\n')
      out.append('WRANGE  = '+' '.join(map(str,wrange))+'\n')
      out.append('STRENGTH= '+str(strength)+'\n')
      out.append('LINELIST= '+' '.join(linelist)+'\n')
      out.append('ATOM    = '+atom+'\n')
      out.append('VROT    = '+str(vrot)+'\n')
      out.append('FWHM    = '+str(fwhm)+'\n')
      out.append('VMACRO    = '+str(vmacro)+'\n')
      out.append('STEPROT = '+str(steprot)+'\n')
      out.append('STEPFWHM= '+str(stepfwhm)+'\n')
      out.append('LTE     = '+str(lte)+'\n')
      out.append('ABU     = '+' '.join(map(str,abu))+'\n')

      header = ''.join(out)

      if synfile == None: 
        tmpstr = os.path.split(modelfile)[-1]
        synfile = tmpstr[:tmpstr.rfind('.')]+'.syn'
      np.savetxt(synfile,(wave,flux,cont),header=header)

  if lineid: 
    if intensity:
      s = wave, flux, cont, inte, continte, [la,li,lo]
    else:
      s = wave, flux, cont, [la,li,lo]
  else:
    if intensity:
      s = wave, flux, cont, inte, continte
    else:
      s = wave, flux, cont

  if tag: tags(s)

  return(s)


def mpsyn(modelfile, wrange, dw=None, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0, intensity=False, lineid=False, tag=False,  \
    clean=True, save=False, synfile=None, \
    lte=False, compute=True, nthreads=0):

  """Computes a synthetic spectrum, splitting the spectral range in nthreads parallel calculations

  Wrapper for syn, using multiprocessing, to speed-up the calculation of a broad spectral range

  Parameters
  ----------
  modelfile : str
      file with a model atmosphere
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  dw: float, optional
      wavelength step for the output fluxes
      this will trigger interpolation at the end
      A negative value will lead to interpolation to a uniform step in ln(lambda)
      (default is None for automatic selection)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  vmicro: float
      microturbulence (km/s)
      a negative value triggers the use from the APOGEE DR14 equation        
      (default is 'model' meaning taken from the model atmosphere)  
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  intensity: bool
      set to True to return intensities  at
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
      -- no broadening due to macro, rotation or instrumental/fwhm effects are considered
      (default False)
  lineid: bool
      set to True to add line identifications to the ouput. They will take the form
      of a list with three arrays (wavelengths, lineids and predicted EWs) 
      (default False)  
  tag: bool
      set to True to produce a plot showing the line identifications. It implicitly
      sets lineid to True, overiding whatever value is passed to lineid
      (default False)
  clean: bool
      True by the default, set to False to avoid the removal of the synspec
      temporary files/links (default True)
  save: bool
      set to True to save the computed spectrum to a file (default False)
      the root of the model atmosphere file, with an extension ".syn" will be used
      but see the parameter synfile to change that
  synfile: str
      when save is True, this can be used to set the name of the output file
      (default None)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)
  compute: bool
      set to False to skip the actual synspec run, triggering clean=False
      (default True)
  nthreads: int
      choose the number of cores to use in the calculation
      (default 0, which has the meaning that the code should take all the cores available except 1)

  Returns
  -------
  wave: numpy array of floats
      wavelengths (angstroms)
  flux: numpy array of floats
      flux (H_lambda in ergs/s/cm2/A)
  cont: numpy array of floats
      continuum flux (same units as flux)

  ---- if intensity is True
  inte: 2D numpy array of floats
      intensity (I_lambda in ergs/s/cm2/A/strad) for
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
  continte: 2D numpy array of float
      continuum intensity for the same angles


  ---- if lineid (or tag) is True 
  lalilo: list with three arrays
        la: numpy array of floats
            wavelenghts of lines (angstroms)
        li: numpy array of str
            line identifidation
        lo: numpy array of floats
            predicted equivalent width (miliangstroms)

  """

  from multiprocessing import Pool,cpu_count

  #basic checks on the line list and model atmosphere
  linelist, modelfile = checksynspec(linelist,modelfile)

  #read model atmosphere
  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)

  if vmicro is None: 
      vmicro = vmicro2
  elif vmicro < 0.0:
      #Holtzman et al. 2015, AJ 150, 148
      vmicro = 2.478 - 0.325 * logg


  if abu is None: abu = abu2

  if nthreads == 0: 
    nthreads = int(cpu_count() - 1)

  #override lineid when tag is True
  if tag: lineid = True

  tmpdir = ''.join(random.choices(string.ascii_lowercase + string.digits, k = 16))

  #delta = (wrange[1]-wrange[0])/nthreads #linear
  power=1 # best load balancing splitting in wavelength as (log(lambda))**(1/power)
  delta = ( (np.log10(wrange[1]))**(1./power) - (np.log10(wrange[0]))**(1./power) )/nthreads
  pars = []
  for i in range(nthreads):

    #wrange1 = (wrange[0]+delta*i,wrange[0]+delta*(i+1))
    wrange1 = ( 10.**( ( (np.log10(wrange[0]))**(1./power) +delta*i     )**power ),
                10.**( ( (np.log10(wrange[0]))**(1./power) +delta*(i+1) )**power )   )

    pararr = [modelfile, wrange1, dw, strength, vmicro, abu, \
      linelist, atom, vrot, fwhm, vmacro, \
      steprot, stepfwhm,  lineid, intensity, tag, clean, False, None, lte, \
      compute, tmpdir+'-'+str(i) ]
    pars.append(pararr)

  pool = Pool(nthreads)
  results = pool.starmap(syn,pars)
  pool.close()
  pool.join()

  x = results[0][0]
  y = results[0][1]
  z = results[0][2]
  if lineid: la, li, lo = results[0][3]

  if len(results) > 1:
    for i in range(len(results)-1):
      x = np.concatenate((x, results[i+1][0][1:]) )
      y = np.concatenate((y, results[i+1][1][1:]) )
      z = np.concatenate((z, results[i+1][2][1:]) )
      if lineid: 
        la2, li2, lo2 = results[i+1][3]
        la = np.concatenate((la, la2[1:]) )
        li = np.concatenate((li, li2[1:]) )
        lo = np.concatenate((lo, lo2[1:]) )

  if save == True:

    out = ['MODEL   = '+modelfile+'\n']
    out.append('TEFF    = '+str(teff)+'\n')
    out.append('LOGG    = '+str(logg)+'\n')
    out.append('VMICRO  = '+str(vmicro)+'\n')
    out.append('WRANGE  = '+' '.join(map(str,wrange))+'\n')
    out.append('STRENGTH= '+str(strength)+'\n')
    out.append('LINELIST= '+' '.join(linelist)+'\n')
    out.append('ATOM    = '+atom+'\n')
    out.append('VROT    = '+str(vrot)+'\n')
    out.append('FWHM    = '+str(fwhm)+'\n')
    out.append('VMACRO  = '+str(vmacro)+'\n')
    out.append('STEPROT = '+str(steprot)+'\n')
    out.append('STEPFWHM= '+str(stepfwhm)+'\n')
    out.append('LTE     = '+str(lte)+'\n')
    out.append('ABU     = '+' '.join(map(str,abu))+'\n')

    header = ''.join(out)

    if synfile == None: 
      tmpstr = os.path.split(modelfile)[-1]
      synfile = tmpstr[:tmpstr.rfind('.')]+'.syn'
    np.savetxt(synfile,(x,y,z),header=header)



  if lineid: 
    s = x, y, z, [la,li,lo]
  else:
    s = x, y, z

  if tag: tags(s)

  return(s)


def raysyn(modelfile, wrange, dw=None, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0,  intensity=False, lineid=False, tag=False, \
    clean=True, save=False, synfile=None, \
    lte=False, compute=True, nthreads=0):

  """Computes a synthetic spectrum, splitting the spectral range in nthreads parallel calculations 

  Wrapper for syn, using ray, to speed-up the calculation of a broad spectral range

  Parameters
  ----------
  modelfile : str
      file with a model atmosphere
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  dw: float, optional
      wavelength step for the output fluxes
      this will trigger interpolation at the end
      A negative value will lead to interpolation to a uniform step in ln(lambda)
      (default is None for automatic selection)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  vmicro: float
      microturbulence (km/s)        
      a negative value triggers the use of the APOGEE DR14 equation
      (default is 'model' meaning taken from the model atmosphere)  
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)      
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  intensity: bool
      set to True to return intensities  at 
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
      -- no broadening due to macro, rotation or instrumental/fwhm effects are considered
      (default False)
  lineid: bool
      set to True to add line identifications to the ouput. They will take the form
      of a list with three arrays (wavelengths, lineids and predicted EWs) 
      (default False)  
  tag: bool
      set to True to produce a plot showing the line identifications. It implicitly
      sets lineid to True, overiding whatever value is passed to lineid
      (default False)
  clean: bool
      True by the default, set to False to avoid the removal of the synspec
      temporary files/links (default True)
  save: bool
      set to True to save the computed spectrum to a file (default False)
      the root of the model atmosphere file, with an extension ".syn" will be used
      but see the parameter synfile to change that
  synfile: str
      when save is True, this can be used to set the name of the output file
      (default None)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)
  compute: bool
      set to False to skip the actual synspec run, triggering clean=False
      (default True)
  nthreads: int
      choose the number of cores to use in the calculation
      (default 0, which has the meaning that the code should take all the cores available except 1)

  Returns
  -------
  wave: numpy array of floats
      wavelengths (angstroms)
  flux: numpy array of floats
      flux (H_lambda in ergs/s/cm2/A)
  cont: numpy array of floats
      continuum flux (same units as flux)

  ---- if intensity is True
  inte: 2D numpy array of floats
      intensity (I_lambda in ergs/s/cm2/A/strad) for
      mu=0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0
  continte: 2D numpy array of float
      continuum intensity for the same angles


  ---- if lineid (or tag) is True
  lalilo: list with three arrays
        la: numpy array of floats
            wavelenghts of lines (angstroms)
        li: numpy array of str
            line identifidation
        lo: numpy array of floats
            predicted equivalent width (miliangstroms)

  """

  import psutil
  import ray

  @ray.remote
  def fun(vari,cons):

    wrange,tmpdir = vari

    modelfile,dw,strength,vmicro,abu,linelist, \
    atom,vrot,fwhm,vmacro,steprot,stepfwhm,intensity, \
    lineid,tag,clean,save,synfile,compute = cons

    s = syn(modelfile, wrange, dw, strength, vmicro, abu, \
              linelist, atom, vrot, fwhm, vmacro, \
              steprot, stepfwhm,  intensity,  \
              lineid, tag, clean, save, synfile, \
              lte, compute, tmpdir)

    return(s)


  #basic checks on the line list and model atmosphere
  linelist, modelfile = checksynspec(linelist,modelfile)

  #read model atmosphere
  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)

  if vmicro is None: 
      vmicro = vmicro2
  elif vmicro < 0.0:
      #Holtzman et al. 2015, AJ 150, 148
      vmicro = 2.478 - 0.325 * logg


  if abu is None: abu = abu2

  if nthreads == 0: 
    nthreads = int(psutil.cpu_count(logical=False) - 1)

  if tag: lineid = True

  print('nthreads=',nthreads)

  tmpdir = ''.join(random.choices(string.ascii_lowercase + string.digits, k = 16))

  ray.init(num_cpus=nthreads)

  rest = [ modelfile,dw,strength,vmicro,abu,linelist, \
    atom,vrot,fwhm,vmacro,steprot,stepfwhm,intensity, \
    lineid,tag,clean,False,None,compute ]

  constants = ray.put(rest)

  #delta = (wrange[1]-wrange[0])/nthreads #linear split
  power=1 # best load balancing splitting in wavelength as (log(lambda))**(1/power)
  delta = ( (np.log10(wrange[1]))**(1./power) - (np.log10(wrange[0]))**(1./power) )/nthreads
  pars = []
  for i in range(nthreads):

    #wrange1 = (wrange[0]+delta*i,wrange[0]+delta*(i+1))
    wrange1 = ( 10.**( ( (np.log10(wrange[0]))**(1./power) +delta*i     )**power ),
                10.**( ( (np.log10(wrange[0]))**(1./power) +delta*(i+1) )**power )   )
    folder = tmpdir+'-'+str(i)

    pararr = [wrange1, folder ]
    pars.append(pararr)

  results = ray.get([fun.remote(pars[i],constants) for i in range(nthreads)])

  ray.shutdown()

  x = results[0][0]
  y = results[0][1]
  z = results[0][2]
  if lineid: la, li, lo = results[0][3]

  if len(results) > 1:
    for i in range(len(results)-1):
      x = np.concatenate((x, results[i+1][0][1:]) )
      y = np.concatenate((y, results[i+1][1][1:]) )
      z = np.concatenate((z, results[i+1][2][1:]) )
      if lineid: 
        la2, li2, lo2 = results[i+1][3]
        la = np.concatenate((la, la2[1:]) )
        li = np.concatenate((li, li2[1:]) )
        lo = np.concatenate((lo, lo2[1:]) )

  if save == True:

    out = ['MODEL   = '+modelfile+'\n']
    out.append('TEFF    = '+str(teff)+'\n')
    out.append('LOGG    = '+str(logg)+'\n')
    out.append('VMICRO  = '+str(vmicro)+'\n')
    out.append('WRANGE  = '+' '.join(map(str,wrange))+'\n')
    out.append('STRENGTH= '+str(strength)+'\n')
    out.append('LINELIST= '+' '.join(linelist)+'\n')
    out.append('ATOM    = '+atom+'\n')
    out.append('VROT    = '+str(vrot)+'\n')
    out.append('FWHM    = '+str(fwhm)+'\n')
    out.append('VMACRO  = '+str(vmacro)+'\n')
    out.append('STEPROT = '+str(steprot)+'\n')
    out.append('STEPFWHM= '+str(stepfwhm)+'\n')
    out.append('LTE     = '+str(lte)+'\n')
    out.append('ABU     = '+' '.join(map(str,abu))+'\n')

    header = ''.join(out)

    if synfile == None: 
      tmpstr = os.path.split(modelfile)[-1]
      synfile = tmpstr[:tmpstr.rfind('.')]+'.syn'
    np.savetxt(synfile,(x,y,z),header=header)

  if lineid: 
    s = x, y, z, [la,li,lo]
  else:
    s = x, y, z

  if tag: tags(s)

  return(s)



def multisyn(modelfiles, wrange, dw=None, strength=1e-4, abu=None, \
    vmicro=None, vrot=0.0, fwhm=0.0, vmacro=0.0, nfe=0.0, \
    linelist=linelist0, atom='ap18', \
    steprot=0.0, stepfwhm=0.0, clean=True, save=None, lte=False, nthreads=0):

  """Computes synthetic spectra for a list of files. The values of vmicro, vrot, 
  fwhm, and nfe can be iterables. Whether or not dw is specified the results will be 
  placed on a common wavelength scale by interpolation. When not specified, dw will be 
  chosen as appropriate for the first model in modelfiles.


  Parameters
  ----------
  modelfiles : list of str
      files with model atmospheres
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  dw: float
      wavelength step for the output fluxes.
      A negative value will lead to interpolation to a uniform step in ln(lambda).
      Unlike in 'syn', interpolation to a constant step will always be done
      (default is None for automatic selection based on the first model of the list)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  vmicro: float
      microturbulence (km/s)        
      a negative value triggers the use of the APOGEE DR14 equation
      (default is None to take the value from the model atmosphere)  
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)            
  nfe: float, can be an iterable
      [N/Fe] nitrogen abundance change from the one specified in the array 'abu' (dex)
      (default 0.)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  clean: bool
      True by the default, set to False to avoid the removal of the synspec
      temporary files/links (default True)
  save: bool
      set to True to save the computed spectra to files (default False)
      the root of the model atmosphere file, with an extension ".syn" will be used
      if multiple values of vmicro, vrot, fwhm or nfe are used, their values are
      prepended to the file names 
      (default None)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)
  nthreads: int
      choose the number of cores to use in the calculation
      (default 0, which has the meaning that the code should take all the cores available but one)



  Returns
  -------
  wave: numpy array of floats (1D)
      wavelengths (angstroms)
  flux: numpy array of floats (2D -- as many rows as models input)
      flux (H_lambda in ergs/s/cm2/A)
  cont: numpy array of floats (2D -- as many rows as models input)
      continuum flux (same units as flux)


  """


  #when vmicro, vrot, fwhm or nitrogen are not iterables, we create ones, otherwise we copy them
  try: 
    nvmicro = len(vmicro)
    vmicros = vmicro
  except TypeError:
    nvmicro = 1
    vmicros = [ vmicro ] 
  try: 
    nvrot = len(vrot)
    vrots = vrots
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]   
  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]   
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]
  try: 
    nnfe = len(nfe)
    nnfes = nfe
  except TypeError:
    nnfe = 1
    nfes = [ nfe ] 

  assert (len(modelfiles) > 0), 'multisyn needs at least one model to work with'
  wave = None
  flux = None
  cont = None

  for entry in modelfiles:
    for vmicro1 in vmicros:
      for nfe1 in nfes:

        abu1 = copy.copy(abu)        

        #if need be, adjust nitrogen abundance according to nfe
        if (abs(nfe1) > 1e-7):
          if (abu1 == None):
            linelist, entry = checksynspec(linelist,entry)
            atmostype, teff, logg, vmicro2, abu1, nd, atmos = read_model(entry)
          abu1[6] = abu1[6] * 10.**nfe1

        x, y, z = mpsyn(entry, wrange, dw=None, strength=strength, \
        vmicro=vmicro1, abu=abu1, linelist=linelist, atom=atom, \
        clean=clean, save=save, lte=lte, nthreads=nthreads)

        space = np.mean(np.diff(x))
            
        for vrot1 in vrots:
          for fwhm1 in fwhms:
            for vmacro1 in vmacros:

              if fwhm1> 0. or vrot1 > 0. or vmacro1 > 0.:
                start = time.time()
                print( entry, vmicro1, nfe1, vrot1, fwhm1, vmacro1, space)
                x2, y2 = call_rotin (x, y, vrot1, fwhm1, vmacro1, space, \
                steprot, stepfwhm, clean=False, reuseinputfiles=True)
                z2 = np.interp(x2, x, z)
                end = time.time()
                print('convol ellapsed time ',end - start, 'seconds')
              else:
                x2, y2, z2 = x, y, z

              if (entry == modelfiles[0] and vmicro1 == vmicros[0] and vrot1 == vrots[0] 
                  and fwhm1 == fwhms[0] and vmacro1 == vmacros[0] and nfe1 == nfes[0] ):
                if dw == None: dw = np.median(np.diff(x2))
                if (dw < 0.):
                  ldw=np.abs(dw)/np.mean(wrange)
                  nsamples = int((np.log(wrange[1]) - np.log(wrange[0]))/ldw) + 1
                  wave = np.exp(np.arange(nsamples)*ldw + np.log(wrange[0]))
                else:
                  nsamples = int((wrange[1] - wrange[0])/dw) + 1
                  wave = np.arange(nsamples)*dw + wrange[0]

                flux = np.interp(wave, x2, y2)
                #flux = interp_spl(wave, x2, y2)
                cont = np.interp(wave, x2, z2)
              else:
                flux = np.vstack ( (flux, np.interp(wave, x, y) ) )
                #flux = np.vstack ( (flux, interp_spl(wave, x, y) ) )
                cont = np.vstack ( (cont, np.interp(wave, x, z) ) )


  return(wave, flux, cont)


def polydelta(modelfile, wrange, elem, enhance=0.2, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0,  lte=False):

  """Sets a a dir tree to compute synthetic spectra for an input model, and then as
many spectra as elements are input in the elem array (symbols), increasing their 
abundances for one at a time.


  Parameters
  ----------
  modelfile : str
      file with a model atmosphere
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  elem: array of str
      symbols for the elements for which the abundance will be enhanced
  enhance: float, optional
      abundance enhancement (dex)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  vmicro: float
      microturbulence (km/s)        
      a negative value triggers the use of the APOGEE DR14 equation
      (default is None to take it from the model atmosphere)  
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)      
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)

  Returns
  -------
  A directory tree set up to perform calculations (hyd0000001/2/3...). Each
  folder contains a *.job script. The output spectra go to files named 
  0000001fort.7 in each .
  """

  #synspec does not currently run in parallel
  nthreads = 1

  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)

  if vmicro is None: 
      vmicro = vmicro2
  elif vmicro < 0.0:
      #Holtzman et al. 2015, AJ 150, 148
      vmicro = 2.478 - 0.325 * logg


  if abu is None: abu = abu2


  #change chemical symbols into indices
  symbol, mass, sol = elements()
  index = []
  for entry in elem: index.append(symbol.index(entry))

  idir = 0
  for j in range(len(elem)+1):

      idir = idir + 1
      dir = ( "hyd%07d" % (idir) )
      try:
        if dir != '.': os.mkdir(dir)
      except OSError:
        print( "cannot create dir hyd%07d" % (idir) )
      try:
        os.chdir(dir)
      except OSError:
        print( "cannot change dir to hyd%07d" % (idir) )

      if modelfile == 'missing':
        pass
      else:
        #setup the slurm script
        sfile = dir+".job"
        now=time.strftime("%c")
        s = open(sfile ,"w")
        s.write("#!/bin/bash \n")
        s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
        s.write("#This script was written by synple on "+now+" \n")          
        s.write("#SBATCH  -J "+dir+" \n")
        s.write("#SBATCH  -o "+dir+"_%j.out"+" \n")          
        s.write("#SBATCH  -e "+dir+"_%j.err"+" \n")
        #s.write("#SBATCH  -n "+str(nthreads)+" \n")
        s.write("#SBATCH  --ntasks-per-node="+str(1)+" \n")
        s.write("#SBATCH  --cpus-per-task="+str(1)+" \n")
        s.write("#SBATCH  -t 00:15:00"+" \n") #hh:mm:ss
        s.write("#SBATCH  -D "+os.path.abspath(os.curdir)+" \n")
        s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n\n\n")

      abu3 = abu[:]
      if j > 0: abu3[index[j-1]] = abu3[index[j-1]]*10.**(enhance)

      x, y, z = syn(modelfile, wrange, dw=None, strength=1e-4, vmicro=vmicro, abu=abu3, \
          linelist=linelist, atom=atom, vrot=vrot, fwhm=fwhm, vmacro=vmacro, \
          steprot=steprot, stepfwhm=stepfwhm,  clean=False, lte=lte, \
          compute=False, tmpdir='.')

      s.write(synspec+" < "+"fort.5"+"\n")

      si = open("fort.55",'r')
      for i in range(6): line = si.readline()
      entries = line.split()
      space = float(entries[5])
      si.close()
            
      iconv = 1
      inconv = ("%07dfort.5" % (iconv) )
      outconv = ("'%07dfort.7'" % (iconv) )
      if fwhm> 0. or vrot > 0.:
        f = open(inconv,'w')
        f.write( ' %s %s %s \n' % ("'fort.7'", "'fort.17'", outconv) )
        f.write( ' %f %f %f \n' % (vrot, space, steprot) )
        f.write( ' %f %f %f \n' % (fwhm, stepfwhm, vmacro) )
        print('stepfwhm=',stepfwhm)
        f.write( ' %f %f %i \n' % (wrange[0], wrange[1], 0) )
        f.close()
        s.write(rotin+" < "+inconv+"\n")
      else:
        s.write("cp "+" fort.7 "+outconv[1:-1]+"\n")

      s.close()
      os.chmod(sfile ,0o755)

      try:
          os.chdir('..')
      except OSError:
          print( "cannot exit dir hyd%07d" % (idir) )


  return(None,None,None)

def collectdelta(modelfile, wrange, elem, enhance=0.2, 
    strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0,  lte=False):

  """Collects the spectra, after computed, in a dir tree created with polydelta, and writes them out to an output file (modelfile.dlt)


  Parameters
  ----------
  modelfile : str
      file with a model atmosphere
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  elem: array of str
      symbols for the elements for which the abundance will be enhanced
  enhance: float, optional
      abundance enhancement (dex)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  vmicro: float
      microturbulence (km/s)        
      a negative value triggers the use of the APOGEE DR14 equation
      (default is None to take it from the model atmosphere)  
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)

  Returns
  -------
  No data are return, but a file (modelfile.dlt) is produced with a header, and then
  the wavelength array, the flux for the input abundances, and the perturbed flux with
  enhanced abundance for each of the elements in elem.
  """


  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)

  if vmicro is None: 
      vmicro = vmicro2
  elif vmicro < 0.0:
      #Holtzman et al. 2015, AJ 150, 148
      vmicro = 2.478 - 0.325 * logg


  if abu is None: abu = abu2

  if save:
      out = open(modelfile+'.dlt','w')
      out.write('MODEL   = '+modelfile+'\n')
      out.write('TEFF    = '+str(teff)+'\n')
      out.write('LOGG    = '+str(logg)+'\n')
      out.write('VMICRO  = '+str(vmicro)+'\n')
      out.write('WRANGE  = '+' '.join(map(str,wrange))+'\n')
      out.write('ELEM    = '+' '.join(elem)+'\n')
      out.write('ENHANCE = '+str(enhance)+'\n')
      out.write('STRENGTH= '+str(strength)+'\n')
      out.write('LINELIST= '+' '.join(linelist)+'\n')
      out.write('ATOM    = '+atom+'\n')
      out.write('VROT    = '+str(vrot)+'\n')
      out.write('FWHM    = '+str(fwhm)+'\n')
      out.write('VMACRO  = '+str(vmacro)+'\n')
      out.write('STEPROT = '+str(steprot)+'\n')
      out.write('STEPFWHM= '+str(stepfwhm)+'\n')
      out.write('LTE     = '+str(lte)+'\n')
      out.write('ABU     = '+' '.join(map(str,abu))+'\n')

  idir = 0
  for j in range(len(elem)+1):

      idir = idir + 1
      dir = ( "hyd%07d" % (idir) )

      #print('j,idir=',j,idir)

      try:
        os.chdir(dir)
      except OSError:
        print( "cannot change dir to hyd%07d" % (idir) )

      x, y  = np.loadtxt('0000001fort.7', unpack=True)

      if j == 0:
          xx = x
          if save: 
              np.savetxt(out,[xx], fmt='%12.5e')
              np.savetxt(out,[y], fmt='%12.5e')
      else:
          yy2 = np.interp(xx, x, y)
          #yy2 = interp_spl(xx, x, y)
          if save: np.savetxt(out,[yy2], fmt='%12.5e')


      try:
          os.chdir('..')
      except OSError:
          print( "cannot exit dir hyd%07d" % (idir) )

  out.close()

  return None

def mkflt(dltfile,wavelengths,blocks=[],fwhm=0.0,unit='km/s',outdir='.'):

  """produces FERRE filters from a dlt file (output from collectdelta)


  Parameters
  ----------
  dltfile : str
      file with the fluxes for abundance perturbations in various elements, produced by 
      polydelta+collectdelta
  wavelengths: float array
      wavelength array for which we want the filters to be resampled (angstroms). These 
      should match those of the FERRE grid with which the filters will be used
  blocks: array of 2-element arrays/tuples defining wavelength intervals to be blocked in the filters (set to zero values, i.e. not to be used by FERRE). Should have the same units as the wavelengths array (angstroms).
      (default is [], so no region is to be blocked)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms or km/s)
      (default 0.)
  unit: str
      Units for the FWHM of the Gaussian kernel ('km/s' or 'A')
       (default 'km/s')

  outdir: str
      Folder to place the output flt files
       (default is './')

  Returns
  -------
  No data are return, but a FERRE filter file is written for each of the elements which
  abundance has been perturned in the input dlt file (see ELEM in the header of the file).
  """

  assert (unit == 'km/s' or unit == 'A'),'unit for FWHM must be km/s or A (Angstroms)'

  f = open(dltfile,'r')
  flux = False

  k = 0
  hd = {}
  for line in f:
    #print('line=',line)
    if '=' in line:
      b = line.split('=')
      hd[b[0].strip()] = b[1]
    else:
      if flux:
        y = np.array(line.split(), dtype=float)
        if fwhm > 0.0:
          if unit == 'km/s':
            xc,yc = vgconv(x,y,fwhm)
          else:
            xc,yc = lgconv(x,y,fwhm)
        else:
          xc = x[:]
          yc = y[:]

        y2 = np.interp(wavelengths, xc, yc)
        #y2 = interp_spl(wavelengths, xc, yc)

        if k == 0:
          wrange = list(map(float,hd['WRANGE'].split()))
          assert(np.min(wavelengths) >= wrange[0]),'Attempted to interpolate to wavelengths shorter than the minimum in the input dlt file '+dltfile 
          assert(np.max(wavelengths) <= wrange[1]),'Attempted to interpolate to wavelengths longer than the maximum in the input dlt file '+dltfile 
          yref = y2
        else:
          #f2 = open(elem[k-1]+'.flt','w')
          we = y2/yref/np.median(y2/yref)
          wp = we > 1.0
          we[wp] = 1.0
          if np.min(we) > 0.999:
            we[:] = 1.0
          #else:
            #we = (we - np.min(we))/(np.max(we) - np.min(we))
          we = 1.- we

          #apply blocks
          if len(blocks) > 0:
            for entry in blocks:
              wblock = (wavelengths >= entry[0]) & (wavelengths <= entry[1])
              if len(np.where(wblock)[0]) > 1: we[wblock] = 0.0

          if k == 1:
            yy = we
          else:
            yy = np.vstack ( (yy, we ) )
          #np.savetxt(f2,we, fmt='%12.5e')
          #f2.close()
        k = k + 1
      else:
        x = np.array(line.split(), dtype=float)
        flux = True
        elem = hd['ELEM'].strip().split()

  f.close()

  #print(yy.shape)
 

  #write out the filter files
  k = 0
  total = np.sum(yy,0)
  for el in elem:
    f = open(os.path.join(outdir,el+'.flt'),'w')
    res = 2.* yy[k,:] - total 
    wn = res < 0.
    res[wn] = 0.0 
    np.savetxt(f,res, fmt='%12.5e')
    k = k + 1
    f.close()

  return None


def polysyn(modelfiles, wrange, strength=1e-4, abu=None, \
    vmicro=None, vrot=0.0, fwhm=0.0, vmacro=0.0,  \
    linelist=linelist0, atom='ap18', \
    steprot=0.0, stepfwhm=0.0,  clean=True, save=None, lte=True, 
    nchem=1, **kargs):

  """Sets up a directory tree for computing synthetic spectra for a list of files in 
  parallel. The values of vmicro, vrot, fwhm, and nfe can be iterables. 


  Parameters
  ----------
  modelfiles : list of str
      files with model atmospheres
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  abu: array of floats (99 elements), optional
      chemical abundances relative to hydrogen (N(X)/N(H))
      (default taken from input model atmosphere)
  vmicro: float, can be an iterable 
      microturbulence (km/s)       
      a negative value triggers the use of the APOGEE DR14 equation 
      (default is None to take it from the model atmosphere)
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  nfe: float, can be an iterable
      [N/Fe] nitrogen abundance change from the one specified in the array 'abu' (dex)
      (default 0.)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  clean: bool
      True by the default, set to False to avoid the removal of the synspec
      temporary files/links (default True)
  save: bool
      set to True to save the computed spectra to files (default False)
      the root of the model atmosphere file, with an extension ".syn" will be used
      if multiple values of vmicro, vrot, fwhm or nfe are used, their values are
      prepended to the file names 
      (default None)
  lte: bool
      this flag can be set to True to enforce LTE in NLTE models. MARCS, Kurucz, the 
      class of Phoenix models used here are always LTE models. Tlusty models
      can be LTE or NLTE, and this keyword will ignore the populations and compute
      assuming LTE for a input NLTE Tlusty model.
      (default False)
  nchem: int
      number of combinations of abundances to include in the synspec 
      calculations when kargs are pairs (irregular grids).
      All these will be performed for every model in modelfiles. 
      Which elements to
      vary and the range of values for each are specified through kargs. Note: 
      calculations for regular grids will become irregular when nchem > 1.
      This parameter is updated automatically in the subroutine for regular 
      grids.
      (default 1)
  kargs:  tuples
       For irregular grids with random abundances as many pairs as necessary, 
       indicating the range for elemental variations [X/Fe]
       e.g. Na=(-0.2,0.2), Al=(-0.5, 0.2), ...
       For regular grids as many triplets as necessary, giving the number of 
       steps, the lower limits and the stepsize for elemental variations [X/Fe]
       e.g. Na=(9,-0.2,0.05), ...
       All entries must be pairs or triplets.


  Returns
  -------
  builds a directory tree ready to perform synspec calculations

  """

  #synspec does not currently run in parallel
  nthreads = 1


  #when vmicro, vrot, fwhm or nitrogen are not iterables, we create ones, otherwise we copy them
  try: 
    nvmicro = len(vmicro)
    vmicros = vmicro
  except TypeError:
    nvmicro = 1
    vmicros = [ vmicro ] 
  try: 
    nvrot = len(vrot)
    vrots = vrot
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]   
  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ] 
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]  
  chems = dict() # abundance variations [X/Fe]
  symbols = [] #elemental symbols X
  for entry in list(map(str,kargs.keys())): symbols.append(entry)
  iel = 0
  for entry in kargs.values():
    print(entry)
    if iel == 0:
      n_p = []
      symbol, mass, sol = elements()
      zatom = dict()
      for i in range(len(symbol)):
        zatom[symbol[i]] = i + 1
    assert(len(entry) == 2 or len(entry) == 3),'kargs entries must have 2 (irregular grids with random values) or 3 (regular grids) entries'
    if len(entry) == 2:
      chems[symbols[iel]] = np.random.random_sample(nchem)*(entry[1]-entry[0])+entry[0]
    else:
      chems[symbols[iel]] = np.arange(entry[0])*entry[2]+entry[1]
      n_p.append(entry[0])
    iel += 1
   
  lenen = 2 
  if iel > 0 and len(entry) == 3:
    lenen = 3
    kems = np.array(list(product(*chems.values())))
    nchem = len(kems[:,0])
    iel = 0
    for entry in kargs.keys():
      chems[entry] = kems[:,iel]
      iel += 1



  idir = 0
  ichem = -1
  dirfile = open('dirtree.txt','w')
  for entry in modelfiles:
    for vmicro1 in vmicros:
      for ichem in range(nchem):

        idir = idir + 1
        dir = ( "hyd%07d" % (idir) )
        cadena = str(idir)+' folder='+dir+' model='+entry+' vmicro='+str(vmicro1)
        iel = 0
        for el in symbols:
           cadena = cadena + ' ['+el+'/Fe]='+str(chems[el][ichem]) + ' '
           iel += 1
        cadena = cadena + '\n'
        dirfile.write(cadena)
        try:
          if dir != '.': os.mkdir(dir)
        except OSError:
          print( "cannot create dir hyd%07d" % (idir) )
        try:
          os.chdir(dir)
        except OSError:
          print( "cannot change dir to hyd%07d" % (idir) )

        if entry == 'missing' or os.path.getsize(entry) == 0:
          pass
        else:
          atmostype, teff, logg, vmicro2, abu1, nd, atmos = read_model(entry)
          if teff is None: 
            pass
          else:
            #setup the slurm script
            sfile = dir+".job"
            now=time.strftime("%c")
            s = open(sfile ,"w")
            s.write("#!/bin/bash \n")
            s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
            s.write("#This script was written by synple on "+now+" \n") 
            s.write("#SBATCH  -J "+dir+" \n")
            s.write("#SBATCH  -o "+dir+"_%j.out"+" \n")
            s.write("#SBATCH  -e "+dir+"_%j.err"+" \n")
            #s.write("#SBATCH  -n "+str(nthreads)+" \n")
            s.write("#SBATCH  --ntasks-per-node="+str(1)+" \n")
            s.write("#SBATCH  --cpus-per-task="+str(1)+" \n")
            s.write("#SBATCH  -t 04:00:00"+" \n") #hh:mm:ss
            s.write("#SBATCH  -D "+os.path.abspath(os.curdir)+" \n")
            s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n\n\n")


            abu1 = copy.copy(abu)

            #if need be, adjust nitrogen abundance according to nfe
            if (ichem > -1):
              if (abu1 == None):
                linelist, entry = checksynspec(linelist,entry)
                atmostype, teff, logg, vmicro2, abu1, nd, atmos = read_model(entry)
              iel = 0 
              for el in symbols:
                abu1[zatom[el]-1] = abu1[zatom[el]-1] * 10.**chems[el][ichem]
                iel += 1


            x, y, z = syn(entry, wrange, dw=None, strength=strength, vmicro=vmicro1, \
            abu=abu1, linelist=linelist, atom=atom, lte=lte, clean=False, compute=False, tmpdir=".")

            s.write("cd "+os.path.abspath(os.curdir)+" \n")
            s.write(synspec+" < "+"fort.5"+"\n")

            si = open("fort.55",'r')
            for i in range(6): line = si.readline()
            entries = line.split()
            space = float(entries[5])
            si.close()
            
            iconv = 0
            for vrot1 in vrots:
              for fwhm1 in fwhms:
                for vmacro1 in vmacros:

                  print('iconv=',iconv)

                  dirfile.write(' -- '+str(iconv+1)+' vrot='+str(vrot1)+' fwhm='+str(fwhm1)+'\n')
                  iconv = iconv + 1
                  inconv = ("%07dfort.5" % (iconv) )
                  outconv = ("'%07dfort.7'" % (iconv) )
                  if vrot1 > 0.0 or fwhm1 > 0. or vmacro1 > 0.:
                    f = open(inconv,'w')
                    f.write( ' %s %s %s \n' % ("'fort.7'", "'fort.17'", outconv) )
                    f.write( ' %f %f %f \n' % (vrot1, space, steprot) )
                    f.write( ' %f %f %f \n' % (fwhm1, stepfwhm, vmacro1) )
                    print('stepfwhm=',stepfwhm)
                    f.write( ' %f %f %i \n' % (wrange[0], wrange[1], 0) )
                    f.close()
                    s.write(rotin+" < "+inconv+"\n")
                  else:
                    s.write("cp "+" fort.7 "+outconv[1:-1]+"\n")

            s.close()
            os.chmod(sfile ,0o755)

        try:
          os.chdir('..')
        except OSError:
          print( "cannot exit dir hyd%07d" % (idir) )


  return()


def build_exomol(isosum_folder='isosum_data'):

#   code to build linelists from exomol data using Ivan Hubeny's scripts

    isosum_files = sorted(glob.glob(os.path.join(isosum_folder,'*.5')))

    for entry in isosum_files:
      print(entry,'build_exomol_diatomic(entry)')
      fh = open(entry,'r')
      zs = fh.readline()
      print('first line of ',entry,'=',zs)
      fh.close()
      data = np.loadtxt(entry, skiprows=1, dtype=str)
      if data.ndim == 1: 
        niso = 1
      elif data.ndim == 2:
        niso = data.shape[0]
      else:
        print('error!! the isosum data in ',entry,' does not conform to expectations')
        continue 
 
      
      tmp = os.path.split(entry)
      molecule = tmp[-1][:-2]

      #create folder
      try:
        folder = os.path.join('linelists',molecule)
        os.mkdir(folder)
        print('created folder ',folder)
      except OSError:
        print('cannot create folder ',folder)

      #write the .5 file
      fh = open(os.path.join(folder,molecule+'.5'),'w')
      fh.write(zs)
      if niso == 1:
          fh.write('  '+' '.join(data)+'\n')
      else:
          for i in range(niso):
              fh.write('  '+' '.join(data[i,:])+'\n')
      fh.close()

      #write the download.list file
      fh = open(os.path.join(folder,'download.list'),'w')
      for i in range(niso):
          if niso == 1:
              cad = data[2][1:-4]
          else:
              cad = data[i,2][1:-4]
          parts = cad.split('__')
          iso = parts[0]
          tag = parts[1]
          if i == 0:
              url = os.path.join('https://www.exomol.com/db',molecule,tag+'_README.txt')
              fh.write(url+'\n')
          url = os.path.join('https://www.exomol.com/db',molecule,iso,tag,cad+'.states.bz2')
          fh.write(url+'\n')
          url = os.path.join('https://www.exomol.com/db',molecule,iso,tag,cad+'.trans.bz2')
          fh.write(url+'\n')
          url = os.path.join('https://www.exomol.com/db',molecule,iso,tag,cad+'.pf')
          fh.write(url+'\n')
      fh.close()

      #write the cc?.5 files
      zzo = np.array(np.sort(np.array(zs.split(),dtype=int)),dtype=str)
      if int(zzo[1]) < 10: zzo[1] = '0'+zzo[1]
      if molecule[-1] == '+':
          zzs = zzo[0]+zzo[1]+'.01'
      else:
          zzs = zzo[0]+zzo[1]+'.00'
      zz = zs.split()
      for i in range(niso):
          fh = open(os.path.join(folder,'cc'+str(i+1)+'.5'),'w')
          fh.write(str(zzs)+" '"+molecule+"'"+'\n')
          if niso == 1:
              fh.write(' '.join([zz[0],data[0],zz[1],data[1]])+'\n')
              fh.write(data[2][:-4]+"'"+'\n')
          else:
              fh.write(' '.join([zz[0],data[i,0],zz[1],data[i,1]])+'\n')
              fh.write(data[i,2][:-4]+"'"+'\n')
          fh.write('100  100000  -9.0\n')
          fh.close()

      #write the R0 file
      fh = open(os.path.join(folder,'R0'),'w')
      fh.write('awk '+"'"+'{print "wget "$0}'+"'"+' download.list  |sh \n')
      fh.write('bunzip2 *bz2 \n')
      fh.write('ln -s -f ../../xprog/isotops . \n')
      fh.write('../../xprog/isosum.exe < '+molecule+'.5 > '+molecule+'.log \n')
      fh.write('cp fort.10 '+molecule+'.pf \n')
      fh.close()
      os.chmod(os.path.join(folder,'R0') ,0o755)

      #write the R1 file
      fh = open(os.path.join(folder,'R1'),'w')
      for i in range(niso):
          fh.write('echo building cc'+str(i+1)+'.list \n')
          fh.write('../../xprog/list.exe <cc'+str(i+1)+'.5 >cc'+str(i+1)+'.log \n')
          fh.write('../../xprog/reverse.exe <fort.10 >cc'+str(i+1)+'.list \n')
          if i > 0:
              fh.write('echo  merging cc'+str(i)+'.list and cc'+str(i+1)+'.list \n')
              fh.write('ln -s -f cc'+str(i)+'.list fort.10 \n')
              fh.write('../../xprog/merge.exe <cc'+str(i+1)+'.list >cc'+str(i)+str(i+1)+'.log \n')
              fh.write('mv fort.11 cc'+str(i+1)+'.list \n')

      fh.write('mv cc'+str(i+1)+'.list '+molecule+'.list \n') 
      fh.close()
      os.chmod(os.path.join(folder,'R1'),0o755)

    return()



def polyopt(wrange=(9.e2,1.e5), dlw=2.1e-5, binary=False, strength=1e-4, inttab=1, \
    abu=None, linelist=linelist0, \
    tlt = (20,3.08,0.068), tlrho = (20,-14.0,0.59), \
    tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), \
    tofe=(1,0.0,0.0), trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), tvmicro=(1,1.0,0.0), \
    zexclude=None, atom='ap18'):

  """Sets up a directory tree for computing opacity tables for TLUSTY. The table collection forms 
  a regular grid defined by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in (log10 of) temperature (tlt) and (log10 of) density
  (tlrho) are mandatory. Triads in [Fe/H] (tfeh), [alpha/Fe] (tafe), [C/Fe] (tcfe), 
  [N/Fe] (tnfe), [O/Fe] (tofe), [r/Fe] (rfe), and [s/Fe] (sfe) are optional since 
  arrays with just one 0.0 are included by default.

  Parameters
  ----------
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  dlw: float
      step in log10(lambda) for the output opacity table.
      Unlike in 'syn', interpolation to a constant step will always be done
      to build the opacity table. If one wishes to resolve lines, it is advisable
      to use dlw ~ 1e-6 or smaller. However, for tables used to build model atmospheres
      dlw ~ 1e-5 suffices in more cases. NOTE that the actual calculation uses a different
      step, computed internally to ensure that lines are resolved.
      (default value is 2.1e-5)
  binary: boolean
      when true, the output table is written in binary format to speed up
      reading it from tlusty
      (default is False)
  strength: float, optional
      threshold in the line-to-continuum opacity ratio for 
      selecting lines (default is 1e-4)
  inttab: int
      a switch for determining the mode of transformation of opacities from
      the actual calculation to the opacity table:
      1 indicates that the opacities in the table are simply interpolated
      any other integer will cause that the stored opacities are averaged 
      over bins to approximately preserve the integral over wavelength. 
      This latter option can be useful for model construction, where the integral 
      is what matters, and be sufficient to work with values dlw>1e-5.
      (default is 1)
  abu: array of floats (99 elements), optional
    chemical abundances relative to hydrogen (N(X)/N(H))
    (default is solar -- see the function 'elements' )
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines. Give an empty array for a continuum-only (+ H and HeII lines) table
      (default is in array linelist0)
  tlt: tuple
    log10(T) triad (n, llimit, step) for opacity grid
    (default values  chosen for grid lt = np.arange(20)*0.068 + 3.08,
     to cover the range in the DR16 APOGEE MARCS grids)
  tlrho: tuple
    log10(rho) triad (n, llimit, step) for opacity grid
    (default values  chosen for grid lrho = np.arange(20)*0.59 -14.0,
     to cover the range in the DR16 APOGEE MARCS grids)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tnfe: tuple
    [N/Fe] triad
  tofe: tuple
    [O/Fe] triad
  trfeh: tuple
    [r/Fe] triad (r-elements abundance ratio)
  tsfeh: tuple
    [s.Fe] triad (s-elements abundance ratio)
  tvmidro: tuple
    vmicro triad (km/s)
  zexclude: list
    atomic numbers of the elements whose opacity is NOT to be
    included in the table
    (default None)
  atom: str
      'ap18' -- generic opacities used in Allende Prieto+ 2018
      'yo19' -- restricted set for NLTE calculations for APOGEE 2019 (Osorio+ 2019)
      'hhm' -- continuum opacity is simplified to H and H-
      (default 'ap18')
  """

  #synspec does not currently run in parallel
  nthreads = 1

  #expanding the triads t* into iterables
  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    ncfe = len(tcfe)
    assert (ncfe == 3), 'Error: cfe triad must have three elements (n, llimit, step)'
    cfes = np.arange(tcfe[0])*tcfe[2] + tcfe[1]
  except TypeError:
    print('Error: cfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nnfe = len(tnfe)
    assert (nnfe == 3), 'Error: nfe triad must have three elements (n, llimit, step)'
    nfes = np.arange(tnfe[0])*tnfe[2] + tnfe[1]
  except TypeError:
    print('Error: nfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nofe = len(tofe)
    assert (nofe == 3), 'Error: ofe triad must have three elements (n, llimit, step)'
    ofes = np.arange(tofe[0])*tofe[2] + tofe[1]
  except TypeError:
    print('Error: ofe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nrfe = len(trfe)
    assert (nrfe == 3), 'Error: rfe triad must have three elements (n, llimit, step)'
    rfes = np.arange(trfe[0])*trfe[2] + trfe[1]
  except TypeError:
    print('Error: rfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nsfe = len(tsfe)
    assert (nsfe == 3), 'Error: sfe triad must have three elements (n, llimit, step)'
    sfes = np.arange(tsfe[0])*tsfe[2] + tsfe[1]
  except TypeError:
    print('Error: sfe triad must have three elements (n, llimit, step)')
    return ()
  
  try: 
    nvmicro = len(tvmicro)
    assert (nvmicro == 3), 'Error: vmicro triad must have three elements (n, llimit, step)'
    vmicros = np.arange(tvmicro[0])*tvmicro[2] + tvmicro[1]
  except TypeError:
    print('Error: vmicro triad must have three elements (n, llimit, step)')
    return ()
  

  #ranges for the opacity table
  try: 
    nlt = len(tlt)
    assert (nlt == 3), 'Error: lt triad must have three elements (n, llimit, step)'
    lt = np.arange(tlt[0])*tlt[2] + tlt[1]  #log10(T)
  except TypeError:
    print('Error: tlt triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlrho = len(tlrho)
    assert (nlrho == 3), 'Error: lrho triad must have three elements (n, llimit, step)'
    lrho = np.arange(tlrho[0])*tlrho[2] + tlrho[1] #log10(density)
  except TypeError:
    print('Error: tlrho triad must have three elements (n, llimit, step)')
    return ()

  linelist = checklinelistpath(linelist)

  space = np.mean(wrange) / clight * 2.355 / 3. * np.sqrt(0.1289**2 * np.min(10.**lt) / 100. + np.min(vmicros)** 2 / 2.)

  cutoff0=250.
  if 10.**min(lt) < 4000.: cutoff0=500.
  if 10.**min(lt) < 3500.: cutoff0=1000.
  if 10.**min(lt) < 3000.: cutoff0=1500.    

  
  symbol, mass, sol = elements()
  z_metals = np.arange(97,dtype=int) + 3
  #Ar usually included among alphas in MARCS and not in Kurucz/Meszaros
  z_alphas = np.array([8,10,12,14,16,18,20,22],dtype=int) 
  # rs increases: notes and data below from comments in the MARCS code 
  # (provided by B.Edvardsson) 
  # Fractional r-process abundance for Ga-Bi (r+s simply assumed == 100%) | Date 2000-01-18
  # (Note: Ga-Sr (31-38) was just copied from Kaeppeler et al. 1989, below)
  # s-process from Stellar models: Arlandini C., Kaeppeler F., Wisshak K.,
  # Gallino R., Busso M., Straniero O., 1999, Astrophys J. 525, 886-900
  #   Fractions corrected to the revised meteoritic abundances
  #   of Grevesse N., Sauval A.J. 1998, Space Science Review 85, 161-174  
  # -0.99 is assigned to unstable elements
  z_rs = np.arange(62,dtype=int) + 31
  rfrac= np.array([.43, .47, .81, .85, .39, .47, 
                   .41, .11, .08, .17, .15, .50,-.99, .68, .86, 
                   .54, .80, .48, .65, .35, .75, .83, .80, .80, 
                   .85, .19, .38, .23, .51, .44,-.99, .71, .93, 
                   .85, .93, .85, .92, .83, .87, .67, .80, .44, 
                   .59, .44, .91, .91, .99, .95, .94, .41, .24, 
                   .54, .95,-.99,-.99,-.99,-.99,-.99,-.99, 1.0, 
                   -.99, 1.0], dtype=float)                     



  idir = 0
  for feh in fehs:
    for afe in afes:
      for cfe in cfes:
        for nfe in nfes:
          for ofe in ofes:
            for rfe in rfes:
              for sfe in sfes: 
                for vmicro in vmicros:
                
                  print(feh,afe,cfe,nfe,ofe,rfe,sfe)

                  idir = idir + 1
                  dir = ( "hyd%07d" % (idir) )
                  try:
                    if dir != '.': os.mkdir(dir)
                  except OSError:
                    print( "cannot create dir hyd%07d" % (idir) )
                  try:
                    os.chdir(dir)
                  except OSError:
                    print( "cannot change dir to hyd%07d" % (idir) )

                  #check input parameters are valid
                  imode = checkinput(wrange, vmicro, linelist)

                  #setup the slurm script
                  sfile = dir+".job"
                  now=time.strftime("%c")
                  s = open(sfile ,"w")
                  s.write("#!/bin/bash \n")
                  s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
                  s.write("#This script was written by synple on "+now+" \n") 
                  s.write("#SBATCH  -J "+dir+" \n")
                  s.write("#SBATCH  -o "+dir+"_%j.out"+" \n")
                  s.write("#SBATCH  -e "+dir+"_%j.err"+" \n")
                  #s.write("#SBATCH  -n "+str(nthreads)+" \n")
                  #s.write("#SBATCH  --ntasks-per-node "+str(4)+" \n")
                  s.write("#SBATCH  --ntasks-per-node="+str(1)+" \n")
                  s.write("#SBATCH  --cpus-per-task="+str(1)+" \n")
                  s.write("#SBATCH  -t 48:00:00"+" \n") #hh:mm:ss
                  s.write("#SBATCH  -D "+os.path.abspath(os.curdir)+" \n")
                  s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n\n\n")

                
                  if abu is None: 
                    abu2 = copy.copy(sol)
                  else:
                    abu2 = copy.copy(abu)

                  if (abs(feh) > 1e-7): 
                    for i in range(len(z_metals)): 
                      abu2[z_metals[i] - 1] = abu2[z_metals[i] - 1] * 10.**feh
                  if (abs(afe) > 1e-7): 
                    for i in range(len(z_alphas)):
                      abu2[z_alphas[i] - 1] = abu2[z_alphas[i] - 1] * 10.**afe
                  if (abs(cfe) > 1e-7): abu2[5] = abu2[5] * 10.**cfe
                  if (abs(nfe) > 1e-7): abu2[6] = abu2[6] * 10.**nfe
                  if (abs(ofe) > 1e-7): abu2[7] = abu2[7] * 10.**ofe
                  if (abs(rfe) > 1e-7): 
                      for i in range(len(z_rs)): 
                        if rfrac[i] > 0.0: abu2[z_rs[i] - 1] = abu2[z_rs[i] - 1] * rfrac[i] * 10.**rfe
                  if (abs(sfe) > 1e-7): 
                      for i in range(len(z_rs)): 
                        if rfrac[i] > 0.0: abu2[z_rs[i] - 1] = abu2[z_rs[i] - 1] * (1.0 - rfrac[i]) * 10.**sfe

                  if (len(linelist) == 0): 
                      imode = -4 
                  else: imode = -3



                  write55(wrange,dw=space,imode=imode,iprin=0,inlte=0,hydprf=0,      \
                          cutoff0=cutoff0, strength=strength, vmicro=vmicro, \
                          linelist=linelist)

                  write5(9999.,9.9,abu2,atom)
                  
                  writetas('tas',1,linelist)


                  write2(lt,lrho,wrange, filename='opt.data', \
                         dlw=dlw, binary=binary,strength=strength,inttab=inttab)

                  if zexclude != None: 
                    write3(zexclude)
                    
                  create_links(linelist)
 
                  assert (not os.path.isdir('data')), 'A subdirectory *data* exists in this folder, and that prevents the creation of a link to the data directory for synple'

                  #link the data folder (synspec data + default tlusty model atoms)
                  os.symlink(modelatomdir,'./data')          #data directory  

                  
                  s.write('time ' + synspec + " < "+"fort.5"+"\n")
                  s.close()
                  os.chmod(sfile ,0o755)
                  
                  try:
                    os.chdir('..')
                  except OSError:
                    print( "cannot exit dir hyd%07d" % (idir) )		  

  return()
  
def merge_slurm(path='./',ext='slurm',nmerge=2,concurrent=False):
  """identifies all the *.slurm files in the path and merge them 
     in groups of nmerge so that there are fewer/longer jobs. 
     The scripts are named job-*.slurm and written to the current folder

  Parameters
  ----------
  path: str
     path under which the slurm jobs are to be found
     (default is './')
  ext: str
     extension of the slurm jobs to be found
     (default is 'slurm')
  nmerge: int
     size of the groups to be created
     (default is 2)
  concurrent: bool
     whether or not the grouped jobs are to be executed concurrently
     (default is False, so the jobs in a group are to be executed serially)
     
  Returns
  -------
  None
  
  """

  slurms = glob.glob(os.path.join(path,'**','*'+ext), recursive=True)

  nfiles = len(slurms)

  assert nfiles > 0, 'There are no input files ending in '+ext

  k = 0 
  wtime = -1
  print('nfiles=',nfiles)
  for i in range(nfiles):
    f1 = open(slurms[i],'r')
    j = i % nmerge
    if j == 0:
      k = k + 1
      if k > 1: 
        if wtime > -1:
          if concurrent: time = int(time*3.) #factor 3 is a safety margin
          entries = header[wtime].split('=')
          header[wtime] = entries[0]+'='+str(time)+'\n'
        f2.writelines(header)
        if concurrent: body.append("wait\n")
        f2.writelines(body)
        f2.close()
      f2 = open('job-'+"{:04d}".format(k)+'.slurm','w')
      time = 0
      header = []
      body = []
    if concurrent: 
      body.append("(\n")
    for line in f1: 
      if line[0] == "#":
        if j == 0: header.append(line)
        if '--time' in line:
          entries = line.split('=') 
          newtime = int(entries[1])
          if concurrent: 
            if newtime > time: time = newtime
          else:
            time = time + newtime
          if j == 0: wtime = len(header)-1
      else:
        body.append(line)
    if concurrent:
      body.append(") & \n")

  if wtime > -1: 
    if concurrent: time = int(time*3.) #factor 3 is a safety margin
    entries = header[wtime].split('=')
    header[wtime] = entries[0]+'='+str(time)+'\n' 
  f2.writelines(header)
  if concurrent: body.append("wait\n")
  f2.writelines(body)
  f2.close()
    

  print(slurms)

  return None
  
  
def merge_slurm_parallel(path='./',ext='slurm',nmerge=2,ncpu=2):
  """identifies all the *.slurm files in the path and merge them 
     in groups of nmerge so that there are fewer/longer jobs. 
     Inside of the new jobs, gnu parallel will be used to have ncpu
     of the input jobs run simultaneously.
     The new scripts are named job-*.slurm and written to the current 
     folder

  Parameters
  ----------
  path: str
     path under which the slurm jobs are to be found
     (default is './')
  ext: str     extension of the slurm jobs to be found
     (default is 'slurm')
  nmerge: int
     size of the groups to be created
     (default is 2)
     
  ncpu: int
     number of the input jobs to be run simultaneously 
     as part of the output jobs
     
  Returns
  -------
  None
  
  """

  slurms = glob.glob(os.path.join(path,'**','*'+ext), recursive=True)

  nfiles = len(slurms)

  assert nfiles > 0, 'There are no input files ending in '+ext

  k = 0 
  wtime = -1
  print('nfiles=',nfiles)
  for i in range(nfiles):
    f1 = open(slurms[i],'r')
    j = i % nmerge
    if j == 0:
      k = k + 1
      if k > 1: 
        if wtime > -1:
          #if concurrent: time = int(time*3.) #factor 3 is a safety margin
          entries = header[wtime].split('=')
          header[wtime] = entries[0]+'='+str(time)+'\n'
        f2.writelines(header)
        f2.write('module load gnuparallel\n')
        f2.writelines('parallel -j'+str(ncpu)+" :::: "+"input-"+"{:04d}".format(k-1)+".txt\n")
        for entry in infiles: f3.write(entry+'\n')
        f2.close()
        f3.close()
      f2 = open('job-'+"{:04d}".format(k)+'.slurm','w')
      f3 = open('input-'+"{:04d}".format(k)+'.txt','w')
      time = 0
      header = []
      infiles = []

    infiles.append(slurms[i])

    for line in f1: 
      if line[0] == "#":
        if j == 0: header.append(line)
        if '--time' in line:
          entries = line.split('=') 
          newtime = int(entries[1])
          time = time + newtime
          if j == 0: wtime = len(header)-1

  if wtime > -1: 
    time = int(time/ncpu*1.3) #factor 1.3 is a safety margin
    entries = header[wtime].split('=')
    header[wtime] = entries[0]+'='+str(time)+'\n' 

  f2.writelines(header)
  f2.write('module load gnuparallel\n')
  f2.write('module load gnu\n')
  f2.write('parallel -j'+str(ncpu)+" :::: "+"input-"+"{:04d}".format(k)+".txt\n")
  for entry in infiles: f3.write(entry+'\n')
  f2.close()
  f3.close()
  
    

  print(slurms)

  return None
  
  

def grid_builder(config,  modeldir=modeldir):

    conf = load_conf(config)

    wrange = tuple(map(float,conf['wrange'].split()))
    if 'vmicro' in conf: vmicro = float(conf['vmicro'])
    
    for entry in conf['grids']:
       print('grid=',entry)
       os.mkdir(entry)
       os.chdir(entry)
       if 'vmicro' in conf[entry]:  vmicro = float(conf[entry]['vmicro']) 
       if conf[entry]['type'] == 'marcs':
          files = collect_marcs(modeldir=modeldir, 
                   tteff = tuple(map(float,conf[entry]['tteff'].split())),
                   tlogg = tuple(map(float,conf[entry]['tlogg'].split())), 
                   tfeh  = tuple(map(float,conf[entry]['tfeh'].split())),
                   tafe  = tuple(map(float,conf[entry]['tafe'].split())),
                   tcfe  = tuple(map(float,conf[entry]['tcfe'].split())),
                   ignore_missing_models = True,
                   ext = 'mod.gz')                   
       elif conf[entry]['type'] == 'kurucz':
          files = collect_kurucz(modeldir=modeldir, 
                   tteff = tuple(map(float,conf[entry]['tteff'].split())),
                   tlogg = tuple(map(float,conf[entry]['tlogg'].split())), 
                   tfeh  = tuple(map(float,conf[entry]['tfeh'].split())),
                   tafe  = tuple(map(float,conf[entry]['tafe'].split())),
                   tcfe  = tuple(map(float,conf[entry]['tcfe'].split())),
                   ignore_missing_models = True,
                   ext = 'mod')
       else:
          print('only APOGEE marcs or kurucz models are accepted')
          continue
                             
       polysyn(files, wrange = wrange, vmicro = vmicro )
                   
       os.chdir('..')
       
    return(None)
           

#gather config. info
def load_conf(config='desi-n.yaml',confdir='.'):

  try:
    yfile=open(os.path.join(confdir,config),'r')
  except:
    print('ERROR in load_conf: cannot find the file ',config)
    return(None)
  #conf=yaml.full_load(yfile)
  conf=yaml.load(yfile, Loader=yaml.SafeLoader)
  yfile.close()

  return(conf)


def collect_marcs(modeldir=modeldir, tteff=None, tlogg=None, \
                  tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), \
                  tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), \
                  tofe=(1,0.0,0.0), trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), \
                  files_in_folders = True, \
                  tie_afe=False, ignore_missing_models=False, ext='mod.gz'):

  """Collects all the MARCS models in modeldir that are part of a regular grid defined
  by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff (tteff) and logg
  (tlogg) are mandatory. Triads in [Fe/H] (tfeh), [alpha/Fe] (tafe), [C/Fe] (tcfe), 
  [N/Fe] (tnfe), [O/Fe] (tofe), [r/Fe] (rfe), and [s/Fe] (sfe) are optional since 
  arrays with just one 0.0 are included by default.

  Parameters
  ----------
  modeldir: str
    directory where model atmosphere files are
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tnfe: tuple
    [N/Fe] triad
  tofe: tuple
    [O/Fe] triad
  rfeh: tuple
    [r/Fe] triad (r-elements abundance ratio)
  sfeh: tuple
    [s.Fe] triad (s-elements abundance ratio)
  files_in_folders: bool
    True when model files are organized in subfolders by metallicity 
    (e.g. mod_z+0.00 for solar-metallicity models)
  tie_afe: boolean
    if active, when there is no loop in [alpha/Fe] (n in tafe is 1),
    [alpha/Fe] is tied to [Fe/H]:
    [alpha/Fe] is 0.5, 0.25, and 0. for [Fe/H]<=-1.5, -1 and -0.5, and >=0,
    respectively
    (default: False)    
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list
 
  Returns
  -------
  files: list of str
    file names with MARCS models that are in modeldir and match
    the parameters in the requested grid

  """

  #expanding the triads t* into iterables
  try: 
    nteff = len(tteff)
    assert (nteff == 3), 'Error: Teff triad must have three elements (n, llimit, step)'
    teffs = np.arange(tteff[0])*tteff[2] + tteff[1]
  except TypeError:
    print('Error: Teff triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlogg = len(tlogg)
    assert (nlogg == 3), 'Error: logg triad must have three elements (n, llimit, step)'
    loggs = np.arange(tlogg[0])*tlogg[2] + tlogg[1]
  except TypeError:
    print('Error: logg triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    ncfe = len(tcfe)
    assert (ncfe == 3), 'Error: cfe triad must have three elements (n, llimit, step)'
    cfes = np.arange(tcfe[0])*tcfe[2] + tcfe[1]
  except TypeError:
    print('Error: cfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nnfe = len(tnfe)
    assert (nnfe == 3), 'Error: nfe triad must have three elements (n, llimit, step)'
    nfes = np.arange(tnfe[0])*tnfe[2] + tnfe[1]
  except TypeError:
    print('Error: nfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nofe = len(tofe)
    assert (nofe == 3), 'Error: ofe triad must have three elements (n, llimit, step)'
    ofes = np.arange(tofe[0])*tofe[2] + tofe[1]
  except TypeError:
    print('Error: ofe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nrfe = len(trfe)
    assert (nrfe == 3), 'Error: rfe triad must have three elements (n, llimit, step)'
    rfes = np.arange(trfe[0])*trfe[2] + trfe[1]
  except TypeError:
    print('Error: rfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nsfe = len(tsfe)
    assert (nsfe == 3), 'Error: sfe triad must have three elements (n, llimit, step)'
    sfes = np.arange(tsfe[0])*tsfe[2] + tsfe[1]
  except TypeError:
    print('Error: sfe triad must have three elements (n, llimit, step)')
    return ()

  if not os.path.isabs(modeldir): 
        modeldir = os.path.join (os.getcwd(), modeldir)

  files = []

  fi = open('files.txt','w')

  for teff in teffs:
    for logg in loggs:
      for feh in fehs:
        for afe in afes:
          for cfe in cfes:
            for nfe in nfes:
              for ofe in ofes:
                for rfe in rfes:
                  for sfe in sfes: 
					  
                    if tie_afe and len(afes) == 1: 
                        if feh <= -1.5: 
                            afe = 0.5
                        elif feh <= -0.4:
                            afe = 0.25
                        else:
                            afe = 0.0
                            
                    print(teff,logg,feh,afe,cfe,nfe,ofe,rfe,sfe)
                    code = 'm*_t*_x3'

                    if logg >= 3.5: 
                      a1 = 'p' 
                    else: 
                      a1 = 's'

                    #make [O/Fe] = [alpha/Fe]
                    if ((tofe[0] == 1) and (abs(tofe[1]) < 1e-7) and (abs(tofe[2]) < 1e-7) ):
                      ofe = afe 

                    sformat = "%s%4i_g%+.1f_%s_z%+.2f_a%+.2f_c%+.2f_n%+.2f_o%+.2f_r%+.2f_s%+.2f."+ext
                    filename = (sformat % (a1,teff,logg,code,feh,afe,cfe,nfe,ofe,rfe,sfe) )

                    if files_in_folders:
                      folder = ("mod_z%+.2f" % (feh) )
                    else:
                      folder = ''

                    file = glob.glob(os.path.join(modeldir,folder,filename))
                    print(filename,file)

                    assert len(file) < 2, 'Multiple files matching the pattern'
                    if len(file) < 1: 
                      file = 'missing'
                    else:
                      file = file[0]

                    if ignore_missing_models == False:
                      assert file != 'missing'
                      assert os.path.isfile(file), 'Cannot find model '+filename+' in modeldir '+modeldir                   
                    else:
                      if file != 'missing': 
                        if not os.path.isfile(file): file = 'missing'
                      
                    files.append(file) 

                    fi.write( "%s  %4i %+.1f %+.2f %+.2f %+.2f %+.2f %+.2f %+.2f %+.2f\n" % (files[-1],teff,logg,feh,afe,cfe,nfe,ofe,rfe,sfe) )



  fi.close()

  return(files)

def collect_kurucz(modeldir=modeldir, tteff=None, tlogg=None, tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), \
  tcfe=(1,0.0,0.0), tie_afe=False, ignore_missing_models=False, ext='mod'):

  """Collects all the (APOGEE ATLAS9) Kurucz models in modeldir that are part of a regular grid defined
  by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff (tteff) and logg
  (tlogg) are mandatory. Triads in [Fe/H] (tfeh), [alpha/Fe] (tafe), [C/Fe] (tcfe) are optional since 
  arrays with just one 0.0 are included by default.

  Parameters
  ----------
  modeldir: str
    directory where model atmosphere files are
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tie_afe: boolean
    if active, when there is no loop in [alpha/Fe] (n in tafe is 1),
    [alpha/Fe] is tied to [Fe/H]:
    [alpha/Fe] is 0.5, 0.25, and 0. for [Fe/H]<=-1.5, -1 and -0.5, and >=0,
    respectively
    (default: False)
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list
 
  Returns
  -------
  files: list of str
    file names with Kurucz models that are in modeldir and match
    the parameters in the requested grid

  """

  #expanding the triads t* into iterables
  try: 
    nteff = len(tteff)
    assert (nteff == 3), 'Error: Teff triad must have three elements (n, llimit, step)'
    teffs = np.arange(tteff[0])*tteff[2] + tteff[1]
  except TypeError:
    print('Error: Teff triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlogg = len(tlogg)
    assert (nlogg == 3), 'Error: logg triad must have three elements (n, llimit, step)'
    loggs = np.arange(tlogg[0])*tlogg[2] + tlogg[1]
  except TypeError:
    print('Error: logg triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    ncfe = len(tcfe)
    assert (ncfe == 3), 'Error: cfe triad must have three elements (n, llimit, step)'
    cfes = np.arange(tcfe[0])*tcfe[2] + tcfe[1]
  except TypeError:
    print('Error: cfe triad must have three elements (n, llimit, step)')
    return ()

  if not os.path.isabs(modeldir): 
        modeldir = os.path.join (os.getcwd(), modeldir)

  files = []

  fi = open('files.txt','w')

  for teff in teffs:
    for logg in loggs:
      for feh in fehs:
        for afe in afes:
          for cfe in cfes:
                
                    if tie_afe and len(afes) == 1: 
                        if feh <= -1.5: 
                            afe = 0.5
                        elif feh <= -0.4:
                            afe = 0.25
                        else:
                            afe = 0.0
                    print(teff,logg,feh,afe,cfe)
                    mcode = 'm'
                    acode = 'm'
                    ccode = 'm'
                    if afe >= 0.: acode='p'
                    if cfe >= 0.: ccode='p'
                    if feh >= 0.: mcode='p'

                    sformat = "t%05ig%3.1fm%s%02ic%s%02io%s%02i."+ext
                    filename = (sformat % (teff,logg,mcode,ceil(abs(feh)*10.),ccode,cfe*10.,acode,ceil(abs(afe)*10.)))

                    file = glob.glob(os.path.join(modeldir,filename))
                    print(filename,file)
                    assert len(file) < 2, 'Multiple files matching the pattern'
                    if len(file) < 1: 
                      file = 'missing'
                    else:
                      file = file[0]

                    print('file=',file)

                    if ignore_missing_models == False:
                      assert file != 'missing'
                      assert os.path.isfile(file), 'Cannot find model '+filename+' in modeldir '+modeldir                   
                    else:
                      if file != 'missing': 
                        if not os.path.isfile(file): file = 'missing'
                      
                    files.append(file)

                    fi.write( "%s  %4i %+.1f %+.2f %+.2f %+.2f \n" % (files[-1],teff,logg,feh,afe,cfe) )



  fi.close()

  return(files)




def collect_k2odfnew(modeldir=modeldir, tteff=None, tlogg=None, tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), \
    tie_afe=False, ignore_missing_models=False):

  """Collects all the ODFNEW Castelli/Kurucz models in modeldir that are part of a regular grid defined
  by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff (tteff) and logg
  (tlogg) are mandatory. Triads in [Fe/H] (tfeh), and [alpha/Fe] (tafe) are optional since 
  arrays with just one 0.0 are included by default. 

  NOTE: There are ODFNEW models with only afe=[alpha/Fe]=0.0 or 0.4. The latter are used whenever
  afe takes values > 0.0, while the afe=0.0 models are used otherwise.

  Parameters
  ----------
  modeldir: str
    directory where model atmosphere files are
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tie_afe: boolean
    if active, when there is no loop in [alpha/Fe] (n in tafe is 1),
    [alpha/Fe] is tied to [Fe/H]:
    [alpha/Fe] is 0.4 and 0. for [Fe/H]<=-1. and any other case, 
    respectively
    (default: False)    
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list
  ext: str
    extension of the model files, usually 'mod' for MARCS but
    could be '.7' or '.22' for Tlusty NLTE models based on MARCS
    (default 'mod')
 
  Returns
  -------
  files: list of str
    file names with Kurucz ODFNEWS models that are in modeldir and match
    the parameters in the requested grid

  """

  #expanding the triads t* into iterables
  try: 
    nteff = len(tteff)
    assert (nteff == 3), 'Error: Teff triad must have three elements (n, llimit, step)'
    teffs = np.arange(tteff[0])*tteff[2] + tteff[1]
  except TypeError:
    print('Error: Teff triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlogg = len(tlogg)
    assert (nlogg == 3), 'Error: logg triad must have three elements (n, llimit, step)'
    loggs = np.arange(tlogg[0])*tlogg[2] + tlogg[1]
  except TypeError:
    print('Error: logg triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  if not os.path.isabs(modeldir): 
        modeldir = os.path.join (os.getcwd(), modeldir)

  files = []

  fi = open('files.txt','w')

  for teff in teffs:
    for logg in loggs:
      for feh in fehs:
        for afe in afes:
			
                    if tie_afe and len(afes) == 1: 
                        if feh <= -1.: 
                            afe = 0.4
                        else:
                            afe = 0.0
                
                    print(teff,logg,feh,afe)
                    code = 'k2odfnew.dat'

                    if afe > 0.0: 
                      a1 = 'a' 
                    else: 
                      a1 = ''

                    if feh < 0.0:
                      a2 = 'am'
                    else:
                      a2 = 'ap'



                    sformat = "t%05ig%.1f%s%02i%s"
                    filename = (sformat % (teff,logg,a2,int(abs(feh)*10),a1+code) )

                    file = os.path.join(modeldir,filename)

                    if ignore_missing_models == False:
                      assert os.path.isfile(file), 'Cannot find model '+filename+' in modeldir '+modeldir                   
                    else:
                      if not os.path.isfile(file): file = 'missing'
                      
                    files.append(file)

                    fi.write( "%s  %4i %+.1f %+.2f %+.2f \n" % (files[-1],teff,logg,feh,afe) )

  fi.close()

  return(files)
  
  
def mkgrid(synthfile=None, tteff=None, tlogg=None, 
           tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0),  
           tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), tofe=(1,0.0,0.0), 
           trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), 
           vmicro=None, nfe=0.0, vrot=0.0, fwhm=0.0, vmacro=0.0, 
           wrange=None, dw=None, logw=0, ignore_missing_models=False, **elements):



  """Collects the synthetic spectra part of a regular grid defined
  by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff (tteff) 
  and logg (tlogg) are mandatory. Triads in [Fe/H] (tfeh), [alpha/Fe] (tafe), 
  [C/Fe] (tcfe), [N/Fe] (tnfe), [O/Fe] (tofe), [r/Fe] (rfe), and 
  [s/Fe] (sfe) are optional since  arrays with just one 0.0 are included by default. 
  The wavelength sampling can be chosen (the spectral range must be limited 
  to the range of the computations), but the default is to take it from 
  the first model.

  Parameters
  ----------
  synthfile: str
    Name of the output FERRE synth file
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tnfe: tuple
    [N/Fe] triad
  tofe: tuple
    [O/Fe] triad
  rfeh: tuple
    [r/Fe] triad (r-elements abundance ratio)
  sfeh: tuple
    [s.Fe] triad (s-elements abundance ratio)
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is None to adopt whatever was used in the calculation)
  nfe: float, can be an iterable
      [N/Fe] nitrogen abundance change from the one specified in the array     
      'abu' (dex)
      (default 0.)
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.0)
  vmacro: float, can be an iterable
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  wrange: tuple or list of two floats, optional
      initial and ending wavelengths (angstroms)
      (default None -- chosen by the code from the first input spectrum)
  dw: float, optional
      wavelength step for the output fluxes
      (default is None for automatic frequency selection)
  logw: int
      parameter that indicates whether the wavelength scale should be 
      linear (0), log10 (1), or log (2)
      (default 0)
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list    
  elements:  tuples
    as many triads as necessary, for other elemental variations [X/Fe]
    e.g. Na=(3,-0.2,0.2), Al=(9, -0.5, 0.1), ...
    
 
  Returns
  -------
  None

  """


  pars = []
  n_p = []
  steps = []
  llimits = []
  #expanding the triads t* into iterables
  try: 
    nteff = len(tteff)
    assert (nteff == 3), 'Error: Teff triad must have three elements (n, llimit, step)'
    teffs = np.arange(tteff[0])*tteff[2] + tteff[1]
    if len(teffs) > 1: 
      pars.append('teff')
      n_p.append(len(teffs))
      steps.append(tteff[2])
      llimits.append(tteff[1])
  except TypeError:
    print('Error: Teff triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlogg = len(tlogg)
    assert (nlogg == 3), 'Error: logg triad must have three elements (n, llimit, step)'
    loggs = np.arange(tlogg[0])*tlogg[2] + tlogg[1]
    if len(loggs) > 1: 
      pars.append('logg')
      n_p.append(len(loggs))      
      steps.append(tlogg[2])
      llimits.append(tlogg[1])
  except TypeError:
    print('Error: logg triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
    if len(fehs) > 1: 
      pars.append('feh')
      n_p.append(len(fehs))      
      steps.append(tfeh[2])
      llimits.append(tfeh[1])      
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
    if len(afes) > 1: 
      pars.append('afe')
      n_p.append(len(afes))      
      steps.append(tafe[2])
      llimits.append(tafe[1])      
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    ncfe = len(tcfe)
    assert (ncfe == 3), 'Error: cfe triad must have three elements (n, llimit, step)'
    cfes = np.arange(tcfe[0])*tcfe[2] + tcfe[1]
    if len(cfes) > 1: 
      pars.append('cfe')
      n_p.append(len(cfes))      
      steps.append(tcfe[2])
      llimits.append(tcfe[1])      
  except TypeError:
    print('Error: cfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nnfe = len(tnfe)
    assert (nnfe == 3), 'Error: nfe triad must have three elements (n, llimit, step)'
    nfes = np.arange(tnfe[0])*tnfe[2] + tnfe[1]
    if len(nfes) > 1: 
      pars.append('nfe')
      n_p.append(len(nfes))      
      steps.append(tnfe[2])
      llimits.append(tnfe[1])      
  except TypeError:
    print('Error: nfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nofe = len(tofe)
    assert (nofe == 3), 'Error: ofe triad must have three elements (n, llimit, step)'
    ofes = np.arange(tofe[0])*tofe[2] + tofe[1]
    if len(ofes) > 1: 
      pars.append('ofe') 
      n_p.append(len(ofes))      
      steps.append(tofe[2])
      llimits.append(tofe[1])      
  except TypeError:
    print('Error: ofe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nrfe = len(trfe)
    assert (nrfe == 3), 'Error: rfe triad must have three elements (n, llimit, step)'
    rfes = np.arange(trfe[0])*trfe[2] + trfe[1]
    if len(rfes) > 1: 
      pars.append('rfe')
      n_p.append(len(rfes))      
      steps.append(trfe[2])
      llimits.append(trfe[1])      
  except TypeError:
    print('Error: rfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nsfe = len(tsfe)
    assert (nsfe == 3), 'Error: sfe triad must have three elements (n, llimit, step)'
    sfes = np.arange(tsfe[0])*tsfe[2] + tsfe[1]
    if len(sfes) > 1: 
      pars.append('sfe')  
      n_p.append(len(sfes))      
      steps.append(tsfe[2])
      llimits.append(tsfe[1])       
  except TypeError:
    print('Error: sfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nvmicro = len(vmicro)
    vmicros = vmicro
    pars.append('vmicro')
    n_p.append(len(vmicros))    
    steps.append(vmicros[1]-vmicros[0])
    llimits.append(vmicros[0])
  except TypeError:
    nvmicro = 1
    vmicros = [ vmicro ] 
  try: 
    nnfe1 = len(nfe)
    nfes1 = nfe
    pars.append('nfe')
    n_p.append(len(nfes1))
    steps.append(nfes1[1]-nfes1[0])
    llimits.append(nfes1[0])    
  except TypeError:
    nnfe1 = 1
    nfes1 = [ nfe ] 
  try: 
    nvrot = len(vrot)
    vrots = vrot
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]   
  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]   
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]

  for entry in elements.keys():
    triad = elements[entry]
    nentry = len(triad)
    assert (nentry == 3), 'Error: element '+entry+' triad must have three elements (n, llimit, step)'
    vals = np.arange(triad[0])*triad[2] + triad[1]
    if len(vals) > 1: 
      pars.append(entry)
      n_p.append(len(vals))
      steps.append(triad[2])
      llimits.append(triad[1])       
    
    
  
  #define indices for grid loops
  ll = []
  ind_n_p = []
  i = 0
  print('pars=',pars)
  for entry in pars:
    ind_n_p.append(i)
    ll.append(np.arange(n_p[i]))
    i = i + 1
  ind = np.array(list(product(*ll)))
  
  print


  hdr = mkhdr(tteff=tteff, tlogg=tlogg, tfeh=tfeh, tafe=tafe, 
              tcfe=tcfe, tnfe=tnfe, tofe=tofe, trfe=trfe, tsfe=tsfe,
              vmicro=vmicro, nfe=nfe, vrot=vrot, fwhm=fwhm, vmacro=vmacro, **elements)

  if os.path.isfile(synthfile): 
    print('Warning -- the output file ',synthfile,' exists and will be overwritten')
    f = open(synthfile,'w')
    f.close()

  f = open(synthfile,'a')

  #look for the first sucessful calculation and define the wavelength for the grid  and write the header
  nfreq = 0
  break_out = False
  idir = 0

  j = 0
  for i in ind:
                        j = j + 1
                        print('line ',j)
                        print(i,steps,llimits)
                        par = i*steps+llimits
                        print(par)
                    
                        dir = ( "hyd%07d" % (j) )

                        iconv = 1
                        if vrot < 1e-7 and  fwhm < 1e-7 and vmacro < 1e-7:
                          outconv = "fort.7" 
                        else:
                          outconv = ("%07dfort.7" % (iconv) )

                        file = os.path.join(dir,outconv)
 
                        if os.path.isfile(file):
                              print('first successful calculation is for idir=',idir)
                              #assert os.path.isfile(file), 'Cannot find model '+file 

                              try:
                                wave, flux = np.loadtxt(file, unpack=True)
                                if wrange is None: 
                                  minwave = np.min(wave)
                                  maxwave = np.max(wave)
                                else:
                                  minwave = wrange[0]
                                  maxwave = wrange[1]

                                if dw is None:
                                  dw = np.median(np.diff(wave))
                        
                                nfreq = np.floor((maxwave - minwave)/dw + 1)
 
                                if logw == 0:
                                  x = minwave + np.arange(nfreq)*dw
                                elif logw == 1:
                                  minwave = np.log10(minwave)
                                  dw = dw/(np.max(wave)+np.min(wave))*2./np.log(10.)
                                  x = minwave + np.arange(nfreq)*dw
                                  x = 10.**x
                                elif logw == 2:
                                  minwave = np.log(minwave)
                                  dw = dw/(np.max(wave)+np.min(wave))*2.
                                  x = minwave + np.arange(nfreq)*dw
                                  x = np.exp(x)
                                else:
                                  print('Error: logw can only be 0, 1 or 2')
                                  sys.exit()

                                hdr['SYNTHFILE_INTERNAL'] = "'"+synthfile+"'"
                                hdr['ID'] = "'"+synthfile[2:]+"'"
                                hdr['NPIX'] = str(int(nfreq))
                                hdr['WAVE'] = str(minwave) + ' ' + str(dw)
                                hdr['LOGW'] = str(int(logw))
                                if fwhm > 0.0:
                                  hdr['RESOLUTION'] = str(np.min(x)/np.max(fwhm))
                                f.write(' &SYNTH\n')
                                for entry in hdr: f.write(' '+entry + ' = ' + hdr[entry] + '\n')
                                f.write(' /\n')
                                break_out = True
                                break

                              except OSError:
                                if ignore_missing_models == False:
                                  print('Error reading file ',file,' ... aborting!')
                                  sys.exit(1)
                                else:
                                  print('Error reading file:', file,' ... skipping it ...')
                                  continue 
  
                              
  assert nfreq > 0, 'could not find a single successful calculation in this grid'
  
  #now read, interpolate and write out the calculations
  j = 0

  for i in ind:
                        j = j + 1
                        print('line ',j)
                        print(i,steps,llimits)
                        par = i*steps+llimits
                        print(par)
                    
                        dir = ( "hyd%07d" % (j) )

                        iconv = 0
                        for vrot1 in vrots:
                          for fwhm1 in fwhms:
                            for vmacro1 in vmacros:
								
                              iconv = iconv + 1

                              if vrot < 1e-7 and  fwhm < 1e-7 and vmacro < 1e-7:
                                outconv = "fort.7"
                              else:
                                outconv = ("%07dfort.7" % (iconv) )

                              file = os.path.join(dir,outconv)
 
                              if os.path.isfile(file):
                                try:
                                  wave, flux = np.loadtxt(file, unpack=True)
                                except OSError:
                                  if ignore_missing_models == False:
                                    print('Cannot read model ',file,' ... aborting!')
                                    sys.exit(1)
                                  else:
                                    wave, flux = (np.array([np.min(x),np.max(x)]), np.array([0.0, 0.0])) 
                              else:
                                if ignore_missing_models == False:
                                  assert os.path.isfile(file), 'Cannot find model '+file                  
                                else:
                                  wave, flux = (np.array([np.min(x),np.max(x)]), np.array([0.0, 0.0]))
                    
                              print('idir,iconv, dw=',idir,iconv,dw)
                              print(wave.shape,flux.shape)
                              y = np.interp(x, wave, flux)
                              print(x.shape,y.shape)
                              #plt.plot(wave,flux,'b',x,y,'.')
                              #plt.show()
                              np.savetxt(f,[y], fmt='%12.5e')

  f.close()

  return(None)


def mkgrid_old(synthfile=None, tteff=None, tlogg=None, 
           tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0),  
           tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), tofe=(1,0.0,0.0), 
           trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), 
           vmicro=None, nfe=0.0, vrot=0.0, fwhm=0.0, vmacro=0.0, 
           wrange=None, dw=None, logw=0, ignore_missing_models=False):



  """Collects the synthetic spectra part of a regular grid defined
  by triads in various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff (tteff) and logg (tlogg) are mandatory. Triads in [Fe/H] (tfeh), [alpha/Fe] (tafe), [C/Fe] (tcfe), [N/Fe] (tnfe), [O/Fe] (tofe), [r/Fe] (rfe), and [s/Fe] (sfe) are optional since  arrays with just one 0.0 are included by default. The wavelength sampling can be chosen (the spectral range must be limited to the range of the computations), but the default is to take it from the first model.

  Parameters
  ----------
  synthfile: str
    Name of the output FERRE synth file
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tnfe: tuple
    [N/Fe] triad
  tofe: tuple
    [O/Fe] triad
  rfeh: tuple
    [r/Fe] triad (r-elements abundance ratio)
  sfeh: tuple
    [s.Fe] triad (s-elements abundance ratio)
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
  nfe: float, can be an iterable
      [N/Fe] nitrogen abundance change from the one specified in the array     
      'abu' (dex)
      (default 0.)
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.0)
  vmacro: float, can be an iterable
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  wrange: tuple or list of two floats, optional
      initial and ending wavelengths (angstroms)
      (default None -- chosen by the code from the first input spectrum)
  dw: float, optional
      wavelength step for the output fluxes
      (default is None for automatic frequency selection)
  logw: int
      parameter that indicates whether the wavelength scale should be 
      linear (0), log10 (1), or log (2)
      (default 0)
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list
 
  Returns
  -------
  None

  """

  #expanding the triads t* into iterables
  try: 
    nteff = len(tteff)
    assert (nteff == 3), 'Error: Teff triad must have three elements (n, llimit, step)'
    teffs = np.arange(tteff[0])*tteff[2] + tteff[1]
  except TypeError:
    print('Error: Teff triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nlogg = len(tlogg)
    assert (nlogg == 3), 'Error: logg triad must have three elements (n, llimit, step)'
    loggs = np.arange(tlogg[0])*tlogg[2] + tlogg[1]
  except TypeError:
    print('Error: logg triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nfeh = len(tfeh)
    assert (nfeh == 3), 'Error: feh triad must have three elements (n, llimit, step)'
    fehs = np.arange(tfeh[0])*tfeh[2] + tfeh[1]
  except TypeError:
    print('Error: feh triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nafe = len(tafe)
    assert (nafe == 3), 'Error: afe triad must have three elements (n, llimit, step)'
    afes = np.arange(tafe[0])*tafe[2] + tafe[1]
  except TypeError:
    print('Error: afe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    ncfe = len(tcfe)
    assert (ncfe == 3), 'Error: cfe triad must have three elements (n, llimit, step)'
    cfes = np.arange(tcfe[0])*tcfe[2] + tcfe[1]
  except TypeError:
    print('Error: cfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nnfe = len(tnfe)
    assert (nnfe == 3), 'Error: nfe triad must have three elements (n, llimit, step)'
    nfes = np.arange(tnfe[0])*tnfe[2] + tnfe[1]
  except TypeError:
    print('Error: nfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nofe = len(tofe)
    assert (nofe == 3), 'Error: ofe triad must have three elements (n, llimit, step)'
    ofes = np.arange(tofe[0])*tofe[2] + tofe[1]
  except TypeError:
    print('Error: ofe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nrfe = len(trfe)
    assert (nrfe == 3), 'Error: rfe triad must have three elements (n, llimit, step)'
    rfes = np.arange(trfe[0])*trfe[2] + trfe[1]
  except TypeError:
    print('Error: rfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nsfe = len(tsfe)
    assert (nsfe == 3), 'Error: sfe triad must have three elements (n, llimit, step)'
    sfes = np.arange(tsfe[0])*tsfe[2] + tsfe[1]
  except TypeError:
    print('Error: sfe triad must have three elements (n, llimit, step)')
    return ()

  try: 
    nvmicro = len(vmicro)
    vmicros = vmicro
  except TypeError:
    nvmicro = 1
    vmicros = [ vmicro ] 
  try: 
    nnfe1 = len(nfe)
    nfes1 = nfe
  except TypeError:
    nnfe1 = 1
    nfes1 = [ nfe ] 
  try: 
    nvrot = len(vrot)
    vrots = vrot
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]   
  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]   
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]


  hdr = mkhdr(tteff=tteff, tlogg=tlogg, tfeh=tfeh, tafe=tafe, 
              tcfe=tcfe, tnfe=tnfe, tofe=tofe, trfe=trfe, tsfe=tsfe,
              vmicro=vmicro, nfe=nfe, vrot=vrot, fwhm=fwhm, vmacro=vmacro)

  if os.path.isfile(synthfile): 
    print('Warning -- the output file ',synthfile,' exists and will be overwritten')
    f = open(synthfile,'w')
    f.close()

  f = open(synthfile,'a')

  #look for the first sucessful calculation and define the wavelength for the grid  and write the header
  nfreq = 0
  break_out = False
  idir = 0

  for teff in teffs:
    for logg in loggs:
      for feh in fehs:
        for afe in afes:
          for cfe in cfes:
            for nfe in nfes:
              for ofe in ofes:
                for rfe in rfes:
                  for sfe in sfes: 
                    for vmicro1 in vmicros:
                      for nfe1 in nfes1:
                
                        print(teff,logg,feh,afe,cfe,nfe,ofe,rfe,sfe,vmicro1,nfe1)

                        idir = idir + 1
                        dir = ( "hyd%07d" % (idir) )

                        iconv = 1
                        outconv = ("%07dfort.7" % (iconv) )
                        file = os.path.join(dir,outconv)
 
                        if os.path.isfile(file):
                              print('first successful calculation is for idir=',idir)
                              assert os.path.isfile(file), 'Cannot find model '+file 
                              wave, flux = np.loadtxt(file, unpack=True)
                              if wrange is None: 
                                minwave = np.min(wave)
                                maxwave = np.max(wave)
                              else:
                                minwave = wrange[0]
                                maxwave = wrange[1]

                              if dw is None:
                                dw = np.median(np.diff(wave))
                        
                              nfreq = np.floor((maxwave - minwave)/dw + 1)
 
                              if logw == 0:
                                x = minwave + np.arange(nfreq)*dw
                              elif logw == 1:
                                minwave = np.log10(minwave)
                                dw = dw/(np.max(wave)+np.min(wave))*2./np.log(10.)
                                x = minwave + np.arange(nfreq)*dw
                                x = 10.**x
                              elif logw == 2:
                                minwave = np.log(minwave)
                                dw = dw/(np.max(wave)+np.min(wave))*2.
                                x = minwave + np.arange(nfreq)*dw
                                x = np.exp(x)
                              else:
                                print('Error: logw can only be 0, 1 or 2')
                                sys.exit()

                              hdr['SYNTHFILE_INTERNAL'] = "'"+synthfile+"'"
                              hdr['ID'] = "'"+synthfile[2:]+"'"
                              hdr['NPIX'] = str(int(nfreq))
                              hdr['WAVE'] = str(minwave) + ' ' + str(dw)
                              hdr['LOGW'] = str(int(logw))
                              if fwhm > 0.0:
                                hdr['RESOLUTION'] = str(np.min(x)/np.max(fwhm))
                              f.write(' &SYNTH\n')
                              for entry in hdr: f.write(' '+entry + ' = ' + hdr[entry] + '\n')
                              f.write(' /\n')
                              break_out = True
                              break
                              
                      if break_out: break
                    if break_out: break
                  if break_out: break
                if break_out: break
              if break_out: break
            if break_out: break
          if break_out: break
        if break_out: break
      if break_out: break
    if break_out: break
						  
  
                              
  assert nfreq > 0, 'could not find a single successful calculation in this grid'
  
  #now read, interpolate and write out the calculations
  idir = 0

  for teff in teffs:
    for logg in loggs:
      for feh in fehs:
        for afe in afes:
          for cfe in cfes:
            for nfe in nfes:
              for ofe in ofes:
                for rfe in rfes:
                  for sfe in sfes: 
                    for vmicro1 in vmicros:
                      for nfe1 in nfes1:
                
                        print(teff,logg,feh,afe,cfe,nfe,ofe,rfe,sfe,vmicro1,nfe1)

                        idir = idir + 1
                        dir = ( "hyd%07d" % (idir) )

                        iconv = 0
                        for vrot1 in vrots:
                          for fwhm1 in fwhms:
                            for vmacro1 in vmacros:
								
                              iconv = iconv + 1
                              outconv = ("%07dfort.7" % (iconv) )
                              file = os.path.join(dir,outconv)
 
                              if os.path.isfile(file):
                                wave, flux = np.loadtxt(file, unpack=True)
                              else:
                                if ignore_missing_models == False:
                                  assert os.path.isfile(file), 'Cannot find model '+file                  
                                else:
                                  wave, flux = (np.array([np.min(x),np.max(x)]), np.array([0.0, 0.0]))
                    
                              print('idir,iconv, dw=',idir,iconv,dw)
                              print(wave.shape,flux.shape)
                              y = np.interp(x, wave, flux)
                              print(x.shape,y.shape)
                              #plt.plot(wave,flux,'b',x,y,'.')
                              #plt.show()
                              np.savetxt(f,[y], fmt='%12.5e')

  f.close()

  return(None)

def mkhdr(tteff=None, tlogg=None, tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), \
              tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), tofe=(1,0.0,0.0),   \
              trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), \
              vmicro=None, nfe=0.0, vrot=0.0, fwhm=0.0, vmacro=0.0, **elements):	  

  """Returns a dictionary for a FERRE regular grid defined by triads in 
  various parameters. Each triad has three values (n, llimit, step)
  that define an array x = np.range(n)*step + llimit. Triads in teff 
  (tteff) and logg (tlogg) are mandatory. Triads in [Fe/H] (tfeh), 
  [alpha/Fe] (tafe), [C/Fe] (tcfe), [N/Fe] (tnfe), [O/Fe] (tofe), 
  [r/Fe] (rfe), and [s/Fe] (sfe) are meant to be associated to changes 
  in the model atmospheres and are optional since  arrays with 
  just one 0.0 are included by default. The keyword arguments 'elements'
  can accommodate any other variations in chemistry, even those in the
  elements already included in the previous parameters, and in any
  arbitrary order.

  Parameters
  ----------
  synthfile: str
    Name of the output FERRE synth file
  tteff: tuple
    Teff triad (n, llimit, step)
  tlogg: tuple
    logg triad (n, llimit, step)
  tfeh: tuple
    [Fe/H] triad
  tafe: tuple
    [alpha/Fe] triad  
  tcfe: tuple
    [C/Fe] triad
  tnfe: tuple
    [N/Fe] triad
  tofe: tuple
    [O/Fe] triad
  rfeh: tuple
    [r/Fe] triad (r-elements abundance ratio)
  sfeh: tuple
    [s.Fe] triad (s-elements abundance ratio)
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
  nfe: float, can be an iterable
      [N/Fe] nitrogen abundance change from the one specified in the array     
      'abu' (dex)
      (default 0.)
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.0)
  vmacro: float, can be an iterable
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  wrange: tuple or list of two floats, optional
      initial and ending wavelengths (angstroms)
      (default None -- chosen by the code from the first input spectrum)
  dw: float, optional
      wavelength step for the output fluxes
      (default is None for automatic frequency selection)
  logw: int
      parameter that indicates whether the wavelength scale should be 
      linear (0), log10 (1), or log (2)
      (default 0)
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list    
  elements:  tuples
    as many triads as necessary, for other elemental variations [X/Fe]
    e.g. Na=(3,-0.2,0.2), Al=(9, -0.5, 0.1), ...
  """    
  
  ndim = 0
  n_p = []
  labels = []
  llimits = []
  steps = []
  if tteff[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tteff[0]))
    labels.append('Teff')
    llimits.append(tteff[1])
    steps.append(tteff[2])
  if tlogg[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tlogg[0]))
    labels.append('logg')
    llimits.append(tlogg[1])
    steps.append(tlogg[2])
  if tfeh[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tfeh[0]))
    labels.append('[Fe/H]')
    llimits.append(tfeh[1])
    steps.append(tfeh[2])
  if tafe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tafe[0]))
    labels.append('[a/Fe]')
    llimits.append(tafe[1])
    steps.append(tafe[2])
  if tcfe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tcfe[0]))
    labels.append('[C/Fe]')
    llimits.append(tcfe[1])
    steps.append(tcfe[2])
  if tnfe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tnfe[0]))
    labels.append('[N/Fe]')
    llimits.append(tnfe[1])
    steps.append(tnfe[2])
  if tofe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tofe[0]))
    labels.append('[O/Fe]')
    llimits.append(tofe[1])
    steps.append(tofe[2])
  if trfe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(trfe[0]))
    labels.append('[r/Fe]')
    llimits.append(trfe[1])
    steps.append(trfe[2])
  if tsfe[0] > 1:
    ndim = ndim + 1
    n_p.append(int(tsfe[0]))
    labels.append('[s/Fe]')
    llimits.append(tsfe[1])
    steps.append(tsfe[2])
  if np.abs(np.max(vmicro)) > 1e-7 and not np.isscalar(vmicro):
    ndim = ndim + 1    
    n_p.append(len(vmicro))
    labels.append('vmicro')
    llimits.append(vmicro[0])
    steps.append(vmicro[1]-vmicro[0])
    dvmicro=np.diff(vmicro)
    if np.max(dvmicro) - np.min(dvmicro) > 1.e-7:
      vmicro = np.log10(vmicro)
      dvmicro=np.diff(vmicro)
      print(vmicro)
      print(dvmicro)
      assert np.max(dvmicro) - np.min(dvmicro) < 1.e-7, 'Vmicro values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10vmicro'
      llimits[-1]=vmicro[0]
      steps[-1]=vmicro[1]-vmicro[0]
  if tnfe[0] > 1:
    ndim = ndim + 1
    n_p.append(len(nfe))
    labels.append('[N/Fe]')
    llimits.append(nfe[0])
    steps.append(nfe[1]-nfe[0])
    dnfe=np.diff(nfe)
    assert np.max(dnfe) - np.min(dnfe) < 1.e-7, '[N/Fe] values are not linearly spaced!'
    

  for entry in elements.keys():
    print('entry=',entry)
    triad = elements[entry]
    nentry = len(triad)
    assert (nentry == 3), 'Error: element '+entry+' triad must have three elements (n, llimit, step)'
    vals = np.arange(triad[0])*triad[2] + triad[1]
    print('vals=',vals)
    if len(vals) > 1: 
      ndim = ndim + 1
      n_p.append(len(vals))
      labels.append(entry)
      llimits.append(triad[1])       
      steps.append(triad[2])
      dvals=np.diff(vals)
      assert np.max(dvals) - np.min(dvals) < 1.e-7, 'values for element '+entry+' are not linearly spaced!'


  try:
    nvrot = len(vrot)
    vrots = vrot
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]	    
  if np.abs(np.max(vrots)) > 1e-7 and nvrot > 1:
    ndim = ndim + 1
    n_p.append(len(vrots))
    labels.append('vrot')
    llimits.append(vrots[0])
    steps.append(vrots[1]-vrots[0])
    dvrot=np.diff(vrots)
    if np.max(dvrot) - np.min(dvrot) > 1.e-7:
      vrots = np.log10(vrots)
      dvrot = np.diff(vrots)
      assert np.max(dvrot) - np.min(dvrot) < 1.e-7, 'Vrot values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10vrot'
      llimits[-1]=vrots[0]
      steps[-1]=vrots[1]-vrots[0]
      
  try:
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]	          
  if np.abs(np.max(fwhms)) > 1e-7 and len(fwhms) > 1:
    ndim = ndim + 1
    n_p.append(len(fwhms))
    labels.append('FWHM')
    llimits.append(fwhms[0])
    steps.append(fwhms[1]-fwhms[0])
    dfwhm=np.diff(fwhms)
    if np.max(dfwhm) - np.min(dfwhm) > 1.e-7:
      fwhms = np.log10(fwhms)
      dfwhm=np.diff(fwhms)
      assert np.max(dfwhm) - np.min(dfwhm) < 1.e-7, 'FWHM values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10FWHM'
      llimits[-1]=fwhms[0]
      steps[-1]=fwhms[1]-fwhms[0]
      
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]	          
  if np.abs(np.min(vmacros)) > 1e-7 and len(vmacros) > 1:
    ndim = ndim + 1
    n_p.append(len(vmacros))
    labels.append('VMACRO')
    llimits.append(vmacros[0])
    steps.append(vmacros[1]-vmacros[0])
    dvmacro=np.diff(vmacros)
    if np.max(dvmacro) - np.min(dvmacro) > 1.e-7:
      vmacros = np.log10(vmacros)
      dvmacro=np.diff(vmacros)
      assert np.max(dvmacro) - np.min(dvmacro) < 1.e-7, 'vmacro values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10vmacro'
      llimits[-1]=vmacros[0]
      steps[-1]=vmacros[1]-vmacros[0]


  pwd=os.path.abspath(os.curdir)
  nowtime=time.ctime(time.time())
  osinfo=os.uname()

                    
  hdr = {}
  hdr['DATE'] = "'"+nowtime+"'"
  hdr['N_OF_DIM'] = str(ndim)
  hdr['N_P'] = '  '.join(map(str,n_p))
  for i in range(ndim): hdr['LABEL('+str(i+1)+")"] = "'"+labels[i]+"'"
  hdr['LLIMITS'] = '  '.join(map(str,llimits))
  hdr['STEPS'] = '  '.join(map(str,steps))
  hdr['COMMENTS1'] = "'mixed and computed with synple-synspec'"
  hdr['COMMENTS2'] = "'"+osinfo[0]+' '+osinfo[2]+'.'+osinfo[4]+' running on '+osinfo[1]+"'"
  hdr['COMMENTS3'] = "'pwd is "+pwd+"'"

  return(hdr)

  
def mkgrid_irregular(synthfile=None, teff=True, logg=True, feh=True,   
           vmicro=None, vrot=0.0, fwhm=0.0, vmacro=0.0, 
           wrange=None, dw=None, logw=0, ignore_missing_models=False,**elem):



  """Collects the synthetic spectra part of an irregular grid. 
   To track changes in Teff, logg and [Fe/H] one can activate the booleans teff,
   logg and feh. To track changes in other elements additional booleans (e.g. Ca=True)
   can be made active at the end of the parameter list (**element). 
   The wavelength sampling can be chosen (the spectral range must be limited 
   to the range of the computations), but the default is to take it from the first model.

  Parameters
  ----------
  synthfile: str
    Name of the output FERRE synth file
  teff: boolean
    Activate to track this parameter
  logg: boolean
    Activate to track this parameter
  feh: boolean
    Activate to track this parameter
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere: 'model')
  vrot: float, can be an iterable
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float, can be an iterable
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.0)
  vmacro: float, can be an iterable
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  wrange: tuple or list of two floats, optional
      initial and ending wavelengths (angstroms)
      (default None -- chosen by the code from the first input spectrum)
  dw: float, optional
      wavelength step for the output fluxes
      (default is None for automatic frequency selection)
  logw: int
      parameter that indicates whether the wavelength scale should be 
      linear (0), log10 (1), or log (2)
      (default 0)
  ignore_missing_models: bool
    set to True to avoid stopping when a model is missing,
    in which case a None is entered in the returning list
    
  **elem: booleans
    Activate to track additional chemical elements
    e.g. Na=True, Ca=True
 
  Returns
  -------
  None

  """

  try: 
    nvmicro = len(vmicro)
    vmicros = vmicro
  except TypeError:
    nvmicro = 1
    vmicros = [ vmicro ]  
  try: 
    nvrot = len(vrot)
    vrots = vrot
  except TypeError:
    nvrot = 1
    vrots = [ vrot ]   
  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]   
  try:
    nvmacro = len(vmacro)
    vmacros = vmacro
  except TypeError:
    nvmacro = 1
    vmacros = [ vmacro ]

  pars = []
  if teff: pars.append('teff')
  if logg: pars.append('logg')
  if feh: pars.append('feh')
  if nvmicro > 1: pars.append('vmicro')
  for entry in elem.keys():
    pars.append('['+entry+'/H]')
  if nvrot > 1: pars.append('vrot')
  if nfwhm > 1: pars.append('fwhm')
  if nvmacro > 1: pars.append('vmacro')

  #track solar reference abundances to set the scale
  symbol, mass, sol = elements()
  solabu = dict()
  for i in range(len(symbol)):
    solabu[symbol[i]] = sol[i]


  hdr = mkhdr_irregular(pars)

  if os.path.isfile(synthfile): 
    print('Warning -- the output file ',synthfile,' exists and will be overwritten')
    f = open(synthfile,'w')
    f.close()

  f = open(synthfile,'a')

  #look for the first sucessful calculation and define the wavelength for the grid  and write the header
  nfreq = 0
  break_out = False
  idir = 0
  
  folders = sorted(glob.glob('hyd*'))
  
  ntot = 0
  for entry in folders:
                    idir = idir + 1
                    dir = ( "hyd%07d" % (idir) )
	         
                    madaffile = os.path.join(entry,'fort.5')
                    if ignore_missing_models == False:
                      assert os.path.isfile(madaffile), 'Cannot find madaf file '+madaffile                  
                    else:
                      if not os.path.isfile(madaffile): continue

                    teff2,logg2,vmicro2,abu = read_madaf(madaffile,startdir=entry)
                    imode, iprin, inmod, inlte, hydprf, wrange, cutoff,  \
                         strength, dw, molls, vmicro1 = read55(os.path.join(entry,'fort.55'))
                    feh2 = np.log10(abu['Fe'])+12-7.50
	                                
                    print(teff2,logg2,feh2,vmicro1,vmicro2)
                    
                    ntot = ntot + 1

                    iconv = 1
                    if vrot < 1e-7 and  fwhm < 1e-7 and vmacro < 1e-7:
                      outconv = "fort.7"
                    else:
                      outconv = ("%07dfort.7" % (iconv) )

                    file = os.path.join(dir,outconv)
                     
                    if break_out == False and os.path.isfile(file) and os.path.getsize(file) > 0:
                      fh = open(file,'r')
                      fdata = fh.read()
                      fh.close()
                      if 'NaN' not in fdata: 
                          print('first successful calculation is for idir=',idir)
                          wave, flux = np.loadtxt(file, unpack=True)
                          if wrange is None: 
                            minwave = np.min(wave)
                            maxwave = np.max(wave)
                          else:
                            minwave = wrange[0]
                            maxwave = wrange[1]

                          if dw is None:
                            dw = np.median(np.diff(wave))
                        
                          nfreq = np.floor((maxwave - minwave)/dw + 1)
 
                          if logw == 0:
                            x = minwave + np.arange(nfreq)*dw
                          elif logw == 1:
                            minwave = np.log10(minwave)
                            dw = dw/(np.max(wave)+np.min(wave))*2./np.log(10.)
                            x = minwave + np.arange(nfreq)*dw
                            x = 10.**x
                          elif logw == 2:
                            minwave = np.log(minwave)
                            dw = dw/(np.max(wave)+np.min(wave))*2.
                            x = minwave + np.arange(nfreq)*dw
                            x = np.exp(x)
                          else:
                            print('Error: logw can only be 0, 1 or 2')
                            sys.exit()

                          hdr['SYNTHFILE_INTERNAL'] = "'"+synthfile+"'"
                          hdr['ID'] = "'"+synthfile[2:]+"'"
                          hdr['NPIX'] = str(int(nfreq))
                          hdr['WAVE'] = str(minwave) + ' ' + str(dw)
                          hdr['LOGW'] = str(int(logw))
                          if fwhm > 0.:
                            hdr['RESOLUTION'] = str(np.min(x)/np.max(fwhm))
  
                          break_out = True
                             
                    #if break_out: pass						  
  
                              
  assert nfreq > 0, 'could not find a single successful calculation in this grid'
  ntot = ntot * nvrot * nfwhm * nvmacro
  hdr['NTOT'] = str(ntot)
  
  #write header
  f.write(' &SYNTH\n')
  for entry in hdr: f.write(' '+entry + ' = ' + hdr[entry] + '\n')
  f.write(' /\n')
  
  
  #now read, interpolate and write out the calculations
  idir = 0
  for entry in folders:
                    idir = idir + 1
                    dir = ( "hyd%07d" % (idir) )
	                
                    madaffile = os.path.join(entry,'fort.5')
                    if ignore_missing_models == False:
                      assert os.path.isfile(madaffile), 'Cannot find madaf file '+madaffile                  
                    else:
                      if not os.path.isfile(madaffile): continue


                    teff2,logg2,vmicro2,abu = read_madaf(madaffile,startdir=entry)

                    imode, iprin, inmod, inlte, hydprf, wrange, cutoff, \
                         strength, dw, molls, vmicro1 = read55(os.path.join(entry,'fort.55'))
                    feh2 = np.log10(abu['Fe']) - np.log10(solabu['Fe'])
	                  
                    print(teff2,logg2,feh2,vmicro1,vmicro2)

                    pars = []
                    if teff: pars.append(teff2)
                    if logg: pars.append(logg2)
                    if feh: pars.append(feh2)
                    if nvmicro > 1: pars.append(vmicro1)
                    for el in elem.keys():
                      pars.append(np.log10(abu[el]) - np.log10(solabu[el]) )
                    iconv = 0
                    for vrot1 in vrots:
                      for fwhm1 in fwhms:
                        for vmacro1 in vmacros:
								
                          iconv = iconv + 1
                          if vrot < 1e-7 and  fwhm < 1e-7 and vmacro < 1e-7:
                            outconv = "fort.7"
                          else:
                            outconv = ("%07dfort.7" % (iconv) )

                          file = os.path.join(dir,outconv)
 
                          if nvrot > 1: pars.append(vrot1)
                          if nfwhm > 1: pars.append(fwhm1)
                          if nvmacro > 1: pars.append(vmacro1)

                          fgood = False
                          if os.path.isfile(file) and os.path.getsize(file) > 0:
                            fh = open(file,'r')
                            fdata = fh.read()
                            fh.close()
                            if 'NaN' not in fdata:
                              fgood = True

                          if fgood == True:
                            wave, flux = np.loadtxt(file, unpack=True)
                          else:
                            if ignore_missing_models == False:
                              sys.exit('Cannot find model '+file+' or it contains no data or NaNs')
                            else:
                              wave, flux = (np.array([np.min(x),np.max(x)]), np.array([0.0, 0.0]))
                 
                          print('idir,iconv, dw=',idir,iconv,dw)
                          print(wave.shape,flux.shape)
                          y = np.interp(x, wave, flux)
                          print(x.shape,y.shape)
                          #plt.plot(wave,flux,'b',x,y,'.')
                          #plt.show()
                          np.savetxt(f,[pars+list(y)], fmt='%12.5e')

  f.close()

  return(None)

def mkhdr_irregular(pars):	  
  
  ndim = len(pars)

  pwd=os.path.abspath(os.curdir)
  nowtime=time.ctime(time.time())
  osinfo=os.uname()
                    
  hdr = {}
  hdr['DATE'] = "'"+nowtime+"'"
  hdr['N_OF_DIM'] = str(ndim)
  for i in range(ndim): hdr['LABEL('+str(i+1)+")"] = "'"+pars[i]+"'"
  hdr['TYPE'] = "'irregular'"
  hdr['COMMENTS1'] = "'mixed and computed with synple-synspec'"
  hdr['COMMENTS2'] = "'"+osinfo[0]+' '+osinfo[2]+'.'+osinfo[4]+' running on '+osinfo[1]+"'"
  hdr['COMMENTS3'] = "'pwd is "+pwd+"'"

  return(hdr)
  
  
#create a regular grid of Kurucz model atmospheres
def create_regular_kurucz(tteff=None, tlogg =None, \
                          tfeh = (1,0.0,0.0), tmicro = (1, 1.0, 0.0), \
                          tie_afe=False, **kargs):
							  
    """Creates scripts to compute a regular grid of Kurucz models using Sbordone's version 
    of ATLAS9. The model grid is defined by triads of various parameters.  Each triad has 
    three values (n, llimit, step) that define an array x = np.range(n)*step + llimit. 
    Triads in teff (tteff) and logg (tlogg) are mandatory. Triads in [Fe/H] (tfeh) and 
    microturbulence (tmicro) are optional since arrays with just one 0.0 are included by 
    default. Any other chemical element can be added with additional triads, e.g. to 
    vary sodium with 3 values [Na/Fe] = -0.2, 0.0 and +0.2 one would add a parameter
    Na=(3,-0.2,0.2). The perl script mkk and Kurucz's atlas9 needs to be installed.
    
    Parameters
    ----------
    tteff: tuple
      Teff triad (n, llimit, step)
    tlogg: tuple
      logg triad (n, llimit, step)
    tfeh: tuple
      [Fe/H] triad
    tmicro: tuple
       microturbulence triad
    tie_afe: boolean
      if active, when there is no loop in [alpha/Fe] (n in tafe is 1),
      [alpha/Fe] is tied to [Fe/H]:
      [alpha/Fe] is 0.5, for [Fe/H]<=-1.5, 0.0 for [Fe/H] >=0 and changes
      linearly in between
      (default: False)    
    kargs:  tuples
       as many triads as necessary, for other elemental variations [X/Fe]
       e.g. Na=(3,-0.2,0.2), Al=(9, -0.5, 0.1), ...
    """
    
							  
    n_p = [tteff[0],tlogg[0], tfeh[0], tmicro[0]]
    llimits = [tteff[1], tlogg[1], tfeh[1], tmicro[1]]
    steps  = [tteff[2], tlogg[2] , tfeh[2], tmicro[2]]
    tags = ['teff', 'logg', 'METALS','MICRO'] 

    for entry in list(map(str,kargs.keys())): tags.append(entry)
    for entry in kargs.values(): 
        print(entry)
        n_p.append(entry[0])
        llimits.append(entry[1])
        steps.append(entry[2])
	
    aa = getaa(n_p)	
    
    for i in range(len(aa[:,0])):
       dir = ( "kur%07d" % (i+1) )
       try:
         os.mkdir(dir)
       except OSError:
         print( "cannot create dir kur%07d" % (i+1) )
       
       #setup the slurm script
       sfile = os.path.join(dir,dir+".job")
       now=time.strftime("%c")
       s = open(sfile ,"w")
       s.write("#!/bin/bash \n")
       s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
       s.write("#This script was written by synple on "+now+" \n")
       s.write("#SBATCH  -J "+dir+" \n")
       s.write("#SBATCH  -o "+dir+"_%j.out"+" \n")
       s.write("#SBATCH  -e "+dir+"_%j.err"+" \n")
       #s.write("#SBATCH  -n "+str(nthreads)+" \n")
       s.write("#SBATCH  --ntasks-per-node="+str(1)+" \n")
       s.write("#SBATCH  --cpus-per-task="+str(1)+" \n")
       s.write("#SBATCH  -t 04:00:00"+" \n") #hh:mm:ss
       s.write("#SBATCH  -D "+os.path.abspath(os.curdir)+" \n")
       s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n\n\n")

       comm = 'mkk '
       for j in [0,1]:
         sst = ('%+.3f   ' % (aa[i,j]*steps[j]+llimits[j]) )
         comm = comm + sst
       for j in range(len(tags)-2):
         sst = ('%+.3f   ' % (aa[i,j+2]*steps[j+2]+llimits[j+2]) )
         comm = comm + tags[j+2] + '=' + sst
         if tie_afe and tags[j+2] == "METALS":
           feh = aa[i,j+2]*steps[j+2]+llimits[j+2]
           if feh <= -1.5: 
             afe = 0.5
           elif feh >= 0.:
             afe = 0.0
           else:
             afe = -1./3. * feh
           alphas = ['O', 'Ne','Mg', 'Si', 'S', 'Ca', 'Ti']
           for entry in alphas:
             sst = ('%+.3f   ' % (afe) )
             comm = comm + entry + '=' + sst
 
       print(sst)
       print(comm)
       s.write(comm+'\n')
       st = os.stat(sfile)
       os.chmod(sfile, st.st_mode | stat.S_IEXEC)

    s.close()
    
    return()

#create an irregular grid of Kurucz model atmospheres
def create_irregular_kurucz(n,pteff=None, plogg =None, \
                          pfeh = (0.0,0.0), pmicro = (1.0, 1.0), \
                          tie_afe=False, **kargs):
							  
    """Creates scripts to compute an iregular grid of Kurucz models using Sbordone's version 
    of ATLAS9. The model grid is defined by pairs of various parameters.  Each pair has 
    two values (llimit, ulimit) that define the sampling interval. 
    Pairs in teff (pteff) and logg (plogg) are mandatory. Pairs in [Fe/H] (pfeh) and 
    microturbulence (pmicro) are optional since arrays with just one value (0.0 for [Fe/H]
    and 1.0 for microturbulence) are included by default. Any other chemical element can 
    be added with additional pairs, e.g. to vary sodium in the range -0.2 <= [Na/Fe] <= +0.2 
    one would add a parameter Na=(-0.2,0.2). The perl script mkk and Kurucz's atlas9 needs 
    to be installed.
    
    Parameters
    ----------
    n: int
      Number of models to produce
    pteff: tuple
      Teff pair (llimit, ulimit)
    plogg: tuple
      logg pair (llimit, ulimit)
    pfeh: tuple
      [Fe/H] pair
    pmicro: tuple
       microturbulence pair
    tie_afe: boolean
      if active, when there is no loop in [alpha/Fe] (n in tafe is 1),
      [alpha/Fe] is tied to [Fe/H]:
      [alpha/Fe] is 0.5, for [Fe/H]<=-1.5, 0.0 for [Fe/H] >=0 and changes
      linearly in between
      (default: False)    
    kargs:  tuples
       as many pairs as necessary, for other elemental variations [X/Fe]
       e.g. Na=(-0.2,0.2), Al=(-0.5, 0.2), ...
    """
    
							      
    teff = np.random.random_sample(n)*(pteff[1]-pteff[0])+pteff[0]
    logg = np.random.random_sample(n)*(plogg[1]-plogg[0])+plogg[0]
    feh = np.random.random_sample(n)*(pfeh[1]-pfeh[0])+pfeh[0]
    micro = np.random.random_sample(n)*(pmicro[1]-pmicro[0])+pmicro[0]
    pars = np.vstack ((teff,logg,feh,micro))
    tags = ['teff', 'logg', 'METALS','MICRO'] 
    for entry in list(map(str,kargs.keys())): tags.append(entry)
    for entry in kargs.values(): 
        print(entry)
        newpar = np.random.random_sample(n)*(entry[1]-entry[0])+entry[0]
        pars = np.vstack ((pars, newpar))
	    
    for i in range(len(pars[0,:])):
       dir = ( "kur%07d" % (i+1) )
       try:
         os.mkdir(dir)
       except OSError:
         print( "cannot create dir kur%07d" % (i+1) )
       
       #setup the slurm script
       sfile = os.path.join(dir,dir+".job")
       now=time.strftime("%c")
       s = open(sfile ,"w")
       s.write("#!/bin/bash \n")
       s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
       s.write("#This script was written by synple on "+now+" \n")
       s.write("#SBATCH  -J "+dir+" \n")
       s.write("#SBATCH  -o "+dir+"_%j.out"+" \n")
       s.write("#SBATCH  -e "+dir+"_%j.err"+" \n")
       #s.write("#SBATCH  -n "+str(nthreads)+" \n")
       s.write("#SBATCH  --ntasks-per-node="+str(1)+" \n")
       s.write("#SBATCH  --cpus-per-task="+str(1)+" \n")
       s.write("#SBATCH  -t 04:00:00"+" \n") #hh:mm:ss
       s.write("#SBATCH  -D "+os.path.abspath(os.curdir)+" \n")
       s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n\n\n")

       comm = 'mkk '

       for j in [0,1]:
         sst = ('%+.3f   ' % (pars[j,i]) )
         comm = comm + sst
       
       for j in range(len(tags)-2):
         sst = ('%+.3f   ' % (pars[j+2,i]) )
         comm = comm + tags[j+2] + '=' + sst
         if tie_afe and tags[j+2] == "METALS":
           feh = pars[j+2,i]
           if feh <= -1.5: 
             afe = 0.5
           elif feh >= 0.:
             afe = 0.0
           else:
             afe = -1./3. * feh
           alphas = ['O', 'Ne','Mg', 'Si', 'S', 'Ca', 'Ti']
           for entry in alphas:
             sst = ('%+.3f   ' % (afe) )
             comm = comm + entry + '=' + sst
         

       print(comm)
       s.write(comm+'\n')
       st = os.stat(sfile)
       os.chmod(sfile, st.st_mode | stat.S_IEXEC)

    s.close()
    
    return()

  
#extract the header of a synthfile
def head_synth(synthfile):
    if synthfile[-6:] == "pickle":
        import pickle
        if not os.path.isfile(synthfile):
          sf = os.path.join(griddir,synthfile)
          if os.path.isfile(sf): synthfile = sf
        print('reading header of the grid '+synthfile+'...')
        file=open(synthfile,'rb')
        header, pars, data = pickle.load(file)
        file.close()
    else:
        meta=0
        multi=0
        if not os.path.isfile(synthfile):
          sf = os.path.join(griddir,synthfile)
          if os.path.isfile(sf): synthfile = sf
        print('reading grid '+synthfile+'...')
        file=open(synthfile,'r')
        line=file.readline()
        header={}
        while (1):
            line=file.readline()
            part=line.split('=')
            if (len(part) < 2): 
              meta=meta+1
              if (meta>multi): 
                if multi>0: multi_header.append(header)
                break
              else:
                if (meta > 0): multi_header.append(header)
                header={}
                line=file.readline()
            else:
              k=part[0].strip()
              v=part[1].strip()
              header[k]=v
              if k == 'MULTI': 
                multi=int(v)
                multi_header=[]
        if (multi > 1): header=multi_header

    return header

#extract the wavelength array for a FERRE synth file
def lambda_synth(synthfile):
    multi_header=head_synth(synthfile)
    if np.ndim(multi_header) == 0: multi_header=[multi_header]
    x = None
    xx=[]
    for header in multi_header:
      if 'WAVE' in header:
        tmp=header['WAVE'].split()
        npix=int(header['NPIX'])
        step=float(tmp[1])
        x0=float(tmp[0])
        x=np.arange(npix)*step+x0
        if header['LOGW']:
          if int(header['LOGW']) == 1: x=10.**x
          if int(header['LOGW']) == 2: x=np.exp(x)   
      if 'WAVELENGTHS' in header:
        x = np.array(header['WAVELENGTHS'].split(),dtype=float)
       
      if x is not None:
        if len(xx) == 0:
          xx.append(x)
        else:
          if (np.size(x) != len(xx)):
            xx.append(x)
          else: 
            if ((xx-x).max() > 0.): xx.append(x)

    if len(xx)>1: x=xx[:]

    return x

#read a synthfile
def read_synth(synthfile,nd=False):
    """
  Reads a FERRE spectral grid from disk. It can be pickled.

  Parameters
  ----------
  synthfile: string
   Name of the FERRE synthfile to read. The code will automatically
   look for it in the current directory but if not present will attempt
   to read it from the synple default griddir.
  nd: boolean
   Keyword to reshape the data array from 2D to the actual
   grid dimensions, given in the array N_P, included in the header

  Returns
  -------
  header: dict
   Dictionary containing the header of the grid
  pars: numpy float array
   2D array with the parameters corresponding to the data in the grid
  data: numpy float array
   2D array with the grid data and dimentions product(N_P) x NPIX
   (multidimensional if the grid has more than 1 parameters and 
    the boolean ND is set to True)
  

    """

    if synthfile[-6:] == "pickle":
        import pickle 
        if not os.path.isfile(synthfile):
          sf = os.path.join(griddir,synthfile)
          if os.path.isfile(sf): synthfile = sf
        print('reading grid '+synthfile+'...')
        file=open(synthfile,'rb')
        header, pars, data = pickle.load(file)
        file.close()
    else:
        #header	
        meta=0
        multi=0
        if not os.path.isfile(synthfile):
          sf = os.path.join(griddir,synthfile)
          if os.path.isfile(sf): synthfile = sf
        print('reading grid '+synthfile+'...')
        file=open(synthfile,'r')
        line=file.readline()
        header={}
        nlines=1
        while (1):
            line=file.readline()
            nlines+=1
            part=line.split('=')
            if (len(part) < 2): 
              meta=meta+1
              if (meta>multi): 
                if multi>0: multi_header.append(header)
                break
              else:
                if (meta > 0): multi_header.append(header)
                header={}
                line=file.readline()
                nlines+=1
            else:
              k=part[0].strip()
              v=part[1].strip()
              header[k]=v.strip("'")
              if k == 'MULTI': 
                multi=int(v)
                multi_header=[]
        if (multi > 1): 
      
          header=multi_header
        file.close()

        if type(header) is list:
            header0 = header[1]
        else:
            header0 = header

        #data
        data=np.loadtxt(synthfile, skiprows=nlines, dtype=float)

        #parameters
        ndim = int(header0['N_OF_DIM']) 
        if ('TYPE' in header0 and 'irregular' in header0['TYPE']):
            pars = data[:,0:ndim]
            data = data[:,ndim:]
        else:
            llimits = np.array(header0['LLIMITS'].split(),dtype=float)
            steps = np.array(header0['STEPS'].split(),dtype=float)
            snp=header0['N_P']
            n_p = tuple(np.array(snp.split(),dtype=int)) + (-1,)
            if nd:  data = np.reshape( data, n_p)
            n_p = n_p[:-1]
            pars = getaa(n_p, dtype=float)*steps + llimits



    return header,pars,data
    
def pickle_synth(synthfile,outsynthfile=None):
    """
    Reads a (text) FERRE grid and rewrites it to disk in pickle
    (binary) format with the extension .pickle
    """
    import pickle

    if outsynthfile is None:
      outsynthfile = synthfile[:-3]+'pickle'

    h, p, d = read_synth(synthfile)
    file = open(outsynthfile, 'wb')
    pickle.dump( (h, p, d), file)
    file.close()
	
    return()
    
def write_synth(synthfile,p,d,hdr=None,irregular=False):
    """
    Writes a FERRE spectral grid to disk

    Parameters
    ----------
    synthfile: str
      Name of the output synthfile
    p: numpy array of floats
       2D array with the parameters 
      (nmodel rows x nparam columns)
    d: numpy array
       2D array with the data (fluxes)
      (nmodel rows x nfreqencies columns)
    hdr: dict
      header
      (made up with basic data if None)
    irregular: bool
      force the output grid to be irregular

    Returns
    -------
    writes out the data to a synthfile


    """
	
    #ndim = d.ndim-1
    #n_p = d.shape[:-1]
    #npix = d.shape[-1]
    
    ndim = len(p[0,:])
    npix = len(d[0,:])
    ntot = len(p[:,0])
    
    pwd=os.path.abspath(os.curdir)
    nowtime=time.ctime(time.time())
    osinfo=os.uname()
   
    if hdr is None:
        #minimal header        
        hdr = {}
        hdr['DATE'] = "'"+nowtime+"'"
        hdr['N_OF_DIM'] = str(ndim)
        #hdr['N_P'] = '  '.join(map(str,n_p))
        for i in range(ndim): hdr['LABEL('+str(i+1)+")"] = "'"+"unknown"+"'"
        hdr['LLIMITS'] = '  '.join(map(str,np.zeros(ndim)))
        hdr['STEPS'] = '  '.join(map(str,np.ones(ndim)))
        hdr['COMMENTS1'] = "'created by write_synth, without axis information'"
        hdr['COMMENTS2'] = "'"+osinfo[0]+' '+osinfo[2]+'.'+osinfo[4]+' running on '+osinfo[1]+"'"
        hdr['COMMENTS3'] = "'pwd is "+pwd+"'"

    else:

        if type(hdr) is list:
          hdr0 = hdr[1]
        else:
          hdr0 = hdr.copy()
          hdr = [hdr]

        ndim = int(hdr0['N_OF_DIM']) 
        #n_p = list(map(int,hdr['N_P'].split()))
        if 'TYPE' in hdr0: 
          if 'irregular' in hdr0['TYPE']:
            irregular = True

    if irregular:
       #drop models with zeros and add the model parameteres to the data
       ds = np.sum(d,1)
       wi = np.where(ds + 1e-31 > 1e-30)[0]
       d = d[wi,:]
       p = p[wi,:]
       d = np.hstack((p,d))
       ntot = len(p[:,0])
       for block in hdr:
         block['TYPE'] = "'irregular'"
         block['NTOT'] = str(ntot)


    fout = open(synthfile,'w')
    for block in hdr:
      fout.write(' &SYNTH\n')
      for entry in block: 
        value = block[entry]
        if any(char.isalpha() for char in value):
          fout.write(' '+entry + ' = ' + "'" + str(value) + "'" + '\n')
        else:
          fout.write(' '+entry + ' = ' + str(value) + '\n')
      fout.write(' /\n')

    #now the data	
    #if ndim > 1:
    #    dd = np.reshape( d, (np.product(n_p), npix) )
    for entry in range(ntot):
        d[entry,:].tofile(fout,sep=" ",format="%0.4e")
        fout.write("\n")

    return(None)			

def paste_synth(synthfile,outsynthfile=None):

    """Pastes various synthfiles with the same parameters but
    different spectral ranges into one synthfile with a multiheader.

    Parameters
    ----------
    synthfile:  iterable of strings
      Names of the synthfiles to merge
    outsynthfile: name of the output synthfile

    Returns
    -------
    hh: dict or list
      header of the merged synthfile
    p: numpy array of floats
      parameters
    dd: numpy array of floats
      data (fluxes)
    """

    k = 0
    synthfile2 = ''
    for entry in synthfile:

        h,p,d = read_synth(entry) 
        synthfile2 += entry
        print('synthfile2=',synthfile2)

        if k == 0:
          if type(h) is list:
            hh = h.copy()
          else:
            hh = [h]
            pp = p.copy()
            dd = d.copy()
        else:
          if type(h) is list:
            for block in h:
              hh.append(block)
          else:
            hh.append(h)
          assert(np.max(pp - p) < 1e-30),'The parameters of the file '+entry+ ' do not match those for preceding synthfiles'
          dd = np.hstack((dd,d))	

        k += 1

    h0 = dict()
    h0['MULTI'] = len(hh)
    h0['ID'] = synthfile2
    h0['COMMENTS1'] = 'merged with synple.paste_synth'
    h0['N_OF_DIM'] = hh[0]['N_OF_DIM']
    if 'NTOT' in hh[0]:
      h0['NTOT'] = hh[0]['NTOT']
    if 'TYPE' in hh[0]: 
      h0['TYPE'] = hh[0]['TYPE']

    hh.insert(0,h0)

    if outsynthfile is None: outsynthfile = synthfile2 + '.dat'
    write_synth(outsynthfile,p,dd,hdr=hh)

    return(hh,p,dd)


def merge_synth(synthfile,outsynthfile=None):

    """Merges various synthfiles with the same parameters and
    the same spectral range/resolution into one irregular grid.

    Parameters
    ----------
    synthfile:  iterable of strings
      Names of the synthfiles to merge
    outsynthfile: name of the output synthfile

    Returns
    -------
    hh: dict
      header of the merged synthfile
    pp: numpy array of floats
      parameters
    dd: numpy array of floats
      data (fluxes)
    """

    k = 0
    synthfile2 = 'merged_'
    for entry in synthfile:

        h,p,d = read_synth(entry) 
        synthfile2 += entry
        print('synthfile2=',synthfile2)

        if k == 0:
          if type(h) is list:
            hh = h.copy()
          else:
            hh = [h]
            pp = p.copy()
            dd = d.copy()
            
          #track LABELS in the first header  
          labels = []
          for key in hh[-1].keys():
            if key[:5] == 'LABEL':
              labels.append(hh[-1][key])
        else:
          if type(h) is list:
            for block in h:
              hh.append(block)
          else:
            hh.append(h)
            
          #check the LABELS are the same in the rest
          labels2 = []
          for key in hh[-1].keys():
            if key[:5] == 'LABEL':
              labels2.append(hh[-1][key])

          assert(np.all(np.array(labels) == np.array(labels2))),'The parameters of the file '+entry+ ' do not match those for preceding synthfiles'           

          
          pp = np.vstack((pp,p))
          dd = np.vstack((dd,d))	

        
        k += 1

    h0 = dict()
    h0['MULTI'] = len(hh)
    h0['ID'] = synthfile2
    h0['COMMENTS1'] = 'merged with synple.merge_synth'
    h0['N_OF_DIM'] = hh[0]['N_OF_DIM']
    if 'NTOT' in hh[0]:
      h0['NTOT'] = hh[0]['NTOT']
    if 'TYPE' in hh[0]: 
      h0['TYPE'] = hh[0]['TYPE']

    hh.insert(0,h0)

    if outsynthfile is None: outsynthfile = synthfile2 + '.dat'
    write_synth(outsynthfile,pp,dd,hdr=hh)

    return(hh,pp,dd)



def fill_synth(d,kernel='thin_plate_spline', neighbors=100):
    """
    Completes data rows with zeros interpolating using RBF from 
    non-zero rows 
    ----------
    d: float
      2D array with data (first dim is number of spectra, 2nd number of frequencies)
    kernel: string
      Type of RBF function (linear, thin_plate_spline, cubic, gaussian ...)
    neighbors: int
      Number of nearest neighbors used to compute the interpolation coefficients
      for each grid point

    Returns
    -------
    a new version of the array d with the rows with zeros replaced by interpolated data

    """
    
    from scipy.interpolate import RBFInterpolator
    
    ndim = d.ndim-1
    n_p = d.shape[:-1]
    nfreq = d.shape[-1]	
    
    ds = np.sum(d, 1)
    wi = np.where(ds + 1e-31 > 1e-30)[0]
    wo = np.where(ds + 1e-31 < 1e-30)[0]

    print('ndim=',ndim)
    print('n_p=',n_p)

    #get loop indices for entries in grid
    iarr = getaa(n_p)
    c = RBFInterpolator(iarr[wi,:], d[wi,:], kernel=kernel, neighbors = neighbors )
    d2 = d.copy()
    d2[wo,:] = c(iarr[wo])
    
    return(d2)
        

def getaa(n_p, dtype=int):
    """
    Generates a matrix with len(n_p) columns and product(n_p)
    rows which can be used to transform ndim nested loops into a single 
    loop

    Parameters
    ----------
    np:	iterable giving the	elements in each dimension

    Returns
    -------
    aa:	- int array  with indices

	"""
		    
    ndim = len(n_p)
    ll = []
    for i in np.arange(ndim):
      ll.append(np.arange(n_p[i], dtype=dtype))
    aa = np.array(list(product(*ll)))
  
    return(aa)

def rbf_get(synthfile, kernel='thin_plate_spline', neighbors=100):
  """Computes RBF coefficients for interpolation in an input FERRE grid
  Parameters
  ----------
  synthfile: string
   Name of the FERRE synthfile to interpolate in
  kernel: string
   Type of RBF function (linear, thin_plate_spline, cubic, gaussian ...)
  neighbors: int
   Number of nearest neighbors used to compute the interpolation coefficients 
   for each grid point

  Returns
  -------
  c: - RBFInterpolator object to interpolate in the grid

  pmin: array with the minimum values for each parameter in the grid 

  ptp: array with the peak-to-peak values for each parameter in the grid (max - min)

  """

  from scipy.interpolate import RBFInterpolator


  h, p, d = read_synth(synthfile)

  #n_p = np.array(h['N_P'].split(),dtype=int)
  #nfreq = int(h['NPIX'])
  #iarr = getaa(n_p)

  pmin = p.min(0)
  ptp  = p.ptp(0)

  p2 = (p - pmin) / ptp

  
  print('deriving interpolation coefficients...')
  #c= RBFInterpolator(iarr, d, kernel=kernel, neighbors = neighbors )
  c= RBFInterpolator(p2, d, kernel=kernel, neighbors = neighbors )

  return(c, pmin, ptp)

def rbf_apply(c,pmin,ptp,par):
  """Interpolates using the RBFInterpolate objects in the array c
   and auxiliary minimum on the parameter range (pmin,ptp) to derive
   fluxes for the parameters in the array par

  Parameters
  ----------
  c: RBFInterpolate object 
   previously
   derived from calling rbf_get
  pmin: array of floats
   minimum values for each of the parameters
  ptp: array of floats
   peak-to-peak values for each of the parameters
  par: array of floats
   array of parameters for interpolation
   rows correspond to different interpolations and columns
   should correspond to the parameters 
  """

  from scipy.interpolate import RBFInterpolator
	
  #map the parameters from physical to indices
  par2 = par.copy()
  for i in range(len(par[0,:])):
    #par2[:,i] = (par[:,i] - llimits[i] ) / steps[i] 
    par2[:,i] = (par[:,i] - pmin[i] ) / ptp[i]


  print('applying coefficients ..')
  res = c(par2)

  return(res)
  

  
def getallt(modelfiles):

  """Collects all the values for temperature, density and electron number density
  in a list of files with model atmospheres

  Parameters
  ----------
  modelfiles : list of str
      files with model atmospheres

  Returns
  -------
  t: list
    list of all temperatures in all the layers of the input model atmospheres    
  rho: list
    list of all values of gas pressure in all the layers of the input model atmospheres
    
  ne: list
    list of all values of electron number density in all the layers of the input model atmospheres

  """

  t = []
  rho = []
  ne = []

  for entry in modelfiles:
    print('reading ',entry)
    teff, logg, vmicro, abu, nd, atmos = read_marcs_model2(entry)
    #atmostype,teff,logg,vmicro,abu,nd,atmos = read_model(entry)
    for value in atmos['t']: t.append(value)
    for value in atmos['rho']: rho.append(value)
    for value in atmos['ne']:  ne.append(value)

  return(t,rho,ne)



def call_rotin(wave=None, flux=None, vrot=0.0, fwhm=0.0, vmacro=0.0, 
  space=1e-2, steprot=0.0, stepfwhm=0.0, 
  clean=True, reuseinputfiles=False, logfile='syn.log'):


  """Convolves a synthetic spectrum with a rotation and/or Gaussian kernel

  Interface to the fortran code rotin.

  Parameters
  ----------
  wave: numpy array of floats
      wavelengths (angstroms)
  flux: numpy array of floats
      flux 
  vrot: float
      projected rotational velocity (km/s)
      (default 0.)
  fwhm: float
      Gaussian broadening: macroturbulence, instrumental, etc. (angstroms)
      (default 0.)
  vmacro: float
      Radial-tangential macroturbulence (km/s)
      (default 0.)
  space: float, optional
      characteristic wavelength scale for variations in the spectrum (angstroms)
      (default is 1e-2)
  steprot: float
      wavelength step for convolution with rotational kernel (angstroms)
      set to 0. for automatic adjustment (default 0.)
  stepfwhm: float
      wavelength step for Gaussian convolution (angstroms)
      set to 0. for automatic adjustment (default 0.)
  clean: bool
      True by the default, set to False to avoid the removal of the rotin
      temporary files (default Tr<ue)
  reuseinputfiles: bool
      set to take the input data from the output synspec file (fort.7) rather than 
      from the input arrays (wave, flux)
  logfile: str
      name of the log file
      (default syn.log)

  Returns
  -------
  wave2: numpy array of floats
      wavelengths (angstroms)
  flux2: numpy array of floats
      flux 


  """
  if reuseinputfiles == False:
    f = open('fort.7','w')
    f2 = open('fort.17','w')
    maxflux = np.max(flux)
    for i in range(len(wave)):
      f.write( ' %f %f \n' % (wave[i], flux[i]) )
      f2.write( ' %f %f \n' % (wave[i], maxflux) )
    f.close()
    f2.close()

  f = open('fort.5','w')
  f.write( ' %s %s %s \n' % ("'fort.7'", "'fort.17'", "'fort.11'") )
  f.write( ' %f %f %f \n' % (vrot, space, steprot) )
  f.write( ' %f %f %f \n' % (fwhm, stepfwhm, vmacro) )
  print('stepfwhm=',stepfwhm)
  f.write( ' %f %f %i \n' % (np.min(wave), np.max(wave), 0) )
  f.close()

  synin = open('fort.5')
  synout = open(logfile,'a')
  p = subprocess.Popen([rotin], stdin=synin, stdout = synout, stderr = synout)
  p.wait()
  synout.flush()
  synout.close()
  synin.close()
  
  assert (os.path.isfile('fort.11')), 'Error: I cannot read the file *fort.11* in '+os.getcwd()+' -- looks like rotin has crashed, please look at syn.log'

  wave2, flux2 = np.loadtxt('fort.11', unpack=True)
  print(len(wave),len(wave2))
  
  if clean == True: cleanup_fort()

  return(wave2, flux2)

def read_model(modelfile):
  
  """Reads a model atmosphere into a structure
  
  Parameters
  ----------  
  modelfile : str
      file with a model atmosphere
      
  Returns
  -------
  atmostype :  str
      type of model atmosphere (kurucz/marcs/phoenix/tlusty)
  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density
  """

  #check
  if not os.path.isfile(modelfile):
    mf = os.path.join(modeldir,modelfile)
    if os.path.isfile(mf): modelfile = mf


  atmostype = identify_atmostype(modelfile)

  if atmostype == 'kurucz':
    teff, logg, vmicro, abu, nd, atmos = read_kurucz_model(modelfile) 
  if atmostype == 'marcs':
    teff, logg, vmicro, abu, nd, atmos = read_marcs_model2(modelfile)
  if atmostype == 'phoenix':
    teff, logg, vmicro, abu, nd, atmos = read_phoenix_model(modelfile)
  if atmostype == 'tlusty':
    teff, logg, vmicro, abu, nd, atmos = read_tlusty_model(modelfile)

  return (atmostype,teff,logg,vmicro,abu,nd,atmos)

def identify_atmostype(modelfile):

  """Idenfies the type of model atmosphere in an input file

  Valid options are kurucz, marcs, tlusty (.7) or phoenix

  Parameters
  ----------
  modelfile: str
      file with a model atmosphere

  Returns
  -------
  atmostype: str
      can take the value 'kurucz', 'marcs', 'tlusty' or 'phoenix' 

  """

  if ('PHOENIX' in modelfile and 'fits' in modelfile): atmostype = 'phoenix'
  else: 
    if modelfile[-3:] == '.gz':
      f = gzip.open(modelfile,'rt')
    else:
      f = open(modelfile,'r')
    line = f.readline()
    #print('modelfile / line=',modelfile,line)
    type(line)
    if ('TEFF' in line): atmostype = 'kurucz'
    else: 
      line = f.readline()
      if ('Teff' in line):
        atmostype = 'marcs'
      else:
        atmostype = 'tlusty'
    f.close()
   
  return(atmostype)

def checksynspec(linelist,modelfile):

  """checking that executables and data are where it should be. Prepend
     default directories to linelist and model atmosphere files when necessary

  Parameters
  ----------
  linelist: array of str
      file names of the line lists to be used. The first string should correspond
      to the atomic line list and is mandatory. The remainder are optional and
      correspond to molecular line lists. All files should be in synspec format.
      (see documentation at http://nova.astro.umd.edu/Synspec43/synspec.html)

  """

  dirs = [synpledir,modelatomdir,linelistdir,bindir]
  for entry in dirs: assert (os.path.isdir(entry)), 'dir '+entry+' missing'

  linelist = checklinelistpath(linelist)

  files = [synspec,rotin]
  for entry in linelist: files.append(entry)
  for entry in files: assert (os.path.isfile(entry)), 'file '+entry+' missing'

  if not os.path.isfile(modelfile):
    mf = os.path.join(modeldir,modelfile)
    if os.path.isfile(mf): modelfile = mf

  #print(modeldir)
  #print(modelfile)
  assert (os.path.isfile(modelfile)),'model atmosphere file '+modelfile+' missing'


  return(linelist,modelfile)

def checklinelistpath(linelist):

  """checking whether the line lists in the array linelists are in the current folder
     or in linelistdir and returning the same array with an absolute path
  """

  i = 0 
  for entry in linelist: 
    if os.path.isfile(entry+".11"): 
      entry = entry+".11" #give preference to the binary file
    elif os.path.isfile(entry):
      pass
    else:
      entry = os.path.join(linelistdir,entry)
      if os.path.isfile(entry+".11"): entry = entry+".11"

    assert(os.path.isfile(entry)), 'The line list '+entry+' is neither accessible in the working directory nor in the linelistdir ('+linelistdir+') folder'
    
    linelist[i] = entry
    if not os.path.isabs(entry): 
      ll = os.path.join (os.getcwd(), entry)
      if os.path.isfile(ll): linelist[i] = ll

    i = i + 1

  return(linelist)

#apply lineid tags to an existing plot 
def tags(s, minew=10., normalized=True):
 
  x = s[0]
  y = s[1]
  z = s[2]

  if normalized:
    y = y/z
    ylabel = 'normalized flux'
  else:
    ylabel = 'flux'

  plt.figure()
  plt.clf()
  plt.ion()
  plt.plot(x,y)
  plt.xlim([np.min(x),np.max(x)])
  plt.ylim([np.min(y)*0.1,np.max(y)*1.1])
  plt.xlabel('wavelength (A)')
  plt.ylabel(ylabel)

  lalilo = s[3] 
  la = lalilo[0]
  li = lalilo[1]
  lo = lalilo[2]

  pli = 0.
  pla = 0.
  pres = 0.0

  for i in range(len(la)):
    if lo[i] >= minew:
      res = min(y[(abs(x-la[i]) < la[i]*1e-5)])
      if li[i] == pli and abs(la[i]-pla) < la[i]*1e-4 and abs(res-pres) < res*0.1:
        pass
      else:
        plt.text(la[i]*(1.-2e-5),res-0.07*max(y),li[i])
        pla=la[i]
        pli=li[i]
        pres=res

  plt.show()

  return (pla,pli)

def checkinput(wrange, vmicro, linelist):

  """checking input parameters from user


  Parameters
  ----------
  wrange: tuple or list of two floats
      initial and ending wavelengths (angstroms)
  vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
  linelist: array of str
      filenames of the line lists, the first one corresponds to 
      the atomic lines and all the following ones (optional) to
      molecular lines
      (default is in array linelist0)

  Returns
  ------
  imode: int
      appropriate value for the variable imode, which specifies whether
      one will use many atomic lines (imode=0), just a few (imode=1),
      or none (H and HeII lines are an exception; imode=2)

  """


  #determine imode
  # imode = 0  is default, atoms and molecules, at least 2 line lists 
  # synple sets IFMOL = 1 in 'tas' when an input molecular line list is used
  # but does not set it when only an atomic line list is given
  # imode = 2 for pure continuum
  # imode = 1 for few-lines mode
  # there are two more values which are possible but not considered in this routine
  # imode = -3 for regular opacity tables (TLUSTY)
  # imode = -4 for continuum-only (+ H and HeII lines) opacity tables


  assert (wrange[1] > wrange[0]),'the ending wavelength for the spectral range must be larger than the starting wavelength'

  if len(linelist) == 0: 
    imode = 2  # no atomic or molecular line list -> pure continuum and no molecules
  else:

    nlines, minlambda, maxlambda = getlinelistrange(linelist[0])

    #check
    if nlines > 10:
      assert (wrange[0] > minlambda and wrange[1] < maxlambda),'wrange exceeds the allow range ('+str(minlambda)+' to '+str(maxlambda)+')'
      imode = 0
    else:
      imode = 1

    assert (vmicro >= 0.0),'vmicro = '+str(vmicro)+' but cannot < 0.'
  
  return(imode)

def getlinelistrange(atomiclinelist):
#finds out min and max wavelengths for a line list

  if atomiclinelist[-3:] == '.11':
    file00 = atomiclinelist[:-3]+'.00'
    assert (os.path.isfile(file00)),'The file '+file00+' reporting linelist statistics is missing'
    f = open(file00,'r')
    lines = f.readlines()
    entries = lines[0].split()
    nlines = np.int64(entries[0])
    entries = lines[1].split()
    minlambda = float(entries[0])*10.
    entries = lines[2].split()
    maxlambda = float(entries[0])*10.
    f.close()
  else:
    f = open(atomiclinelist,'r')
    line = f.readline()
    entries = line.split()
    minlambda = float(entries[0])*10.
    fsize = os.path.getsize(atomiclinelist)
    f.seek(fsize-103)
    line = f.readline()
    f.close()
    entries = line.split()
    maxlambda = float(entries[0])*10.
    nlines = int(0.01 * fsize)


  return(nlines, minlambda,maxlambda)



def writetas(filename,nd,linelist,nonstd=None):
#write non-std input parameters
# input: filename -- str -- name of the non-std. param. file to print
#        nd -- int -- number of layers in the model
#        linelist -- list -- names of the linelist files (atomic first, then one 
#				or more molecular ones
#        nonstd -- dict -- additional entries to add to the tas file 
#
  
  f = open(filename,'w')

  if nonstd is not None:
    for entry in nonstd.keys():
      f.write(entry + "=" + str(nonstd[entry])+"\n")

  f.write("ND= "+str(nd)+" \n")
  if len(linelist) > 1:  f.write("IFMOL= "+one+" \n")
  f.write("TMOLIM= 8000. \n")


  f.close()

  return()

def write3(zexclude):
  
  f = open('fort.3','w')
  for z in zexclude:
    f.write( " %d %10.4e \n" % (z, 0.0) )
  f.close()
  
  return()


def write2(lt,lrho,wrange, filename='opt.data', dlw=2.1e-5, binary=False,strength=1e-4,inttab=1):
#write fort.2 file for creating opacity tables for TLUSTY

  f = open('fort.2','w')
  f.write( " %d %10.4e %10.4e \n" % (len(lt),10.**lt[0],10.**lt[-1]) )
  f.write( " %d \n" % (1) )
  f.write( " %d %10.4e %10.4e \n" % (len(lrho),10.**lrho[0],10.**lrho[-1]) )
  
  nsamples = int( (np.log10(wrange[1]) - np.log10(wrange[0]) )/dlw) + 1 
  f.write( " %d %d %10.4e %10.4e \n" % (nsamples,inttab,wrange[0],wrange[1]) )  
  if binary == True: 
    ibingr = 1
  else:
    ibingr = 0
  filename = "'"+filename+"'"
  f.write( " %s %d \n" % (filename,ibingr) )
  f.close()

  return()

def read55(filename='fort.55'):


  #imode,idst,iprin
  #inmod,zero,ichang,ichemc
  #lyman,zero,zero,zero,zero
  #one,nlte,icontl,inlist,ifhe2
  #ihydpr,ihe1pr,ihe2pr
  #wstart,wend,cutoff,zero,strength,wdist 
  #nmol-linelists unit20 unit21 ...
  #vmicro

  f = open(filename,'r')
  line = f.readline()
  entries = line.split()
  imode = int(entries[0])
  iprin = int(entries[2])
  #f.write(" "+str(imode)+" "+zero+" "+str(iprin)+"\n")
  
  line = f.readline()
  entries = line.split()
  inmod = int(entries[0])
  #f.write(" "+str(inmod)+3*zero+"\n")
  
  line = f.readline()
  #f.write(5*zero+"\n")
  
  line = f.readline()
  entries = line.split()
  inlte = int(entries[0])
  #f.write(one+str(abs(inlte))+zero+str(inlist)+zero+"\n")

  line = f.readline()
  entries = line.split()
  hydprf = int(entries[0])
  #f.write(str(hydprf)+2*zero+"\n")

  line = f.readline()
  entries = line.split()
  wrange = [ float(entries[0]), float(entries[1]) ] 
  cutoff0 = float(entries[2])
  strength = float(entries[4])
  dw = float(entries[5])
  #f.write( ' %f %f %f %i %e %f \n ' % (wrange[0],   wrange[1], cutoff0, 0, strength, dw) )

  line = f.readline()
  molls = map(float,line.split())

  line = f.readline()
  vmicro = float(line)
  f.close()

  return(imode,iprin,inmod,inlte,hydprf,wrange,cutoff0,strength,dw,molls,vmicro)

def write55(wrange,dw=1e-2,imode=0,iprin=0,inlte=0,hydprf=2,cutoff0=200., \
  strength=1e-4,vmicro=1.0, \
  linelist=linelist0, atmostype='kurucz',intensity=False):


  #imode,idst,iprin
  #inmod,zero,ichang,ichemc
  #lyman,zero,zero,zero,zero
  #one,nlte,icontl,inlist,ifhe2
  #ihydpr,ihe1pr,ihe2pr
  #wstart,wend,cutoff,zero,strength,wdist 

  if (atmostype == 'tlusty' or atmostype == 'marcs'): inmod = 1 
  else: inmod = 0

  inlist = 0
  all_inlist = []
  for file in linelist:
    inlist = 0
    if file[-3:] == '.11' : inlist = 1
    all_inlist.append(inlist)
    assert (inlist - all_inlist[0] == 0), 'The line list files must be all either text or binary!'

  f = open('fort.55','w')
  f.write(" "+str(imode)+" "+zero+" "+str(iprin)+"\n")
  f.write(" "+str(inmod)+3*zero+"\n")
  f.write(5*zero+"\n")
  f.write(one+str(abs(inlte))+zero+str(inlist)+zero+"\n")
  f.write(str(hydprf)+2*zero+"\n")
  if imode <= -3:
    f.write( ' %f %f %f %i %e %f \n ' % (wrange[0],  -wrange[1], cutoff0, 0, strength, dw) )
  else:
    f.write( ' %f %f %f %i %e %f \n ' % (wrange[0],   wrange[1], cutoff0, 0, strength, dw) )
  ll = len(linelist)
  if ll < 2: f.write(2*zero)
  else: f.write(str(ll-1) + ' ' + ' '.join(map(str,np.arange(ll-1)+20)))
  f.write("\n")
  f.write( ' %f  \n' % (vmicro) )
  if intensity: 
    f.write( ' %i  %f %i \n' % (-10, 0.0001, 1) )
    f.write( ' %f %f %f %f %f %f %f %f %f %f \n' % (0.0001,0.001,0.01,0.1,0.25,0.4,0.55,0.7,0.85,1.0) )
  f.close()

def write5(teff,logg,abu, atom='ap18', ofile='fort.5', inlte=0, atommode=None, atominfo=None):

  symbol, mass, sol = elements()

  f = open(ofile,'w')
  f.write(' '+str(teff)+" "+str(logg).format('%7.4f')+"       ! TEFF, GRAV \n")
  if inlte == 0:
    f.write(" T  F               ! LTE, GRAY \n")
  else:
    f.write(" F  F               ! LTE, GRAY \n")

  f.write(" 'tas'              ! name of non-standard flags \n")
  f.write(" 50                 ! frequencies \n")

  natom = len(abu)
  f.write(" "+str(natom)+"        ! NATOMS \n")  

  assert (atom == 'hhm' or atom == 'ap18' or atom == 'yo19' or atom == 'test'), 'atom must be one of: hhm/ap18/yo19/test!'
  ex = np.ones(natom)

  #tlusty models provide atommode and atominfo that override the defaults for atom and ex
  if atommode is not None: ex[:len(atommode)] = np.array(atommode)
  if atominfo is not None: atom = 'own' 

  if atom == 'hhm' : 
    zex = [1]  #atomic numbers of elements included explicitly (contributing cont. opacity)
  elif atom == 'yo19':
    zex = [1,11,12,19,20]
  elif atom == 'test': 
    zex = [1,26]
  elif atom == 'ap18': 
    zex = [1,2,6,7,8,11,12,13,14,20,26]
  else:
    zex = []

  for i in zex: ex[i-1] = 2

  for i in range(natom):
    f.write(' %2d %e %i %s\n' %  (np.abs(ex[i]), 
				  abu[i], 0, '  ! ' +symbol[i]) )

  for i in range(3): f.write("* \n")
  
  if atom == 'hhm':  # highly simplified continuum opacities -- just H and H-
    f.write("* ../data_atom for ions  \n")
    f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat' \n" )
    f.write("   0    0     3      0 \n")
    f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
    f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
    f.write("   0    0     0     -1     0     0    '    ' ' '  \n")
    f.write("* \n")
    f.write("* end \n")
  elif atom == "yo19": # set for NLTE calculations for APOGEE (see Osorio+ 2019 A&A paper)
    f.write("* ../data_atom for ions  \n")
    f.write("  1    -1     1      0     0     1    ' H 0' 'data/hm.dat'  \n")
    f.write("  0     0     3      0   \n")
    f.write("  1     0     16     0     0     0    ' H 1' 'data/h1_16lev2.dat'  \n")
    f.write("  1     1     1      1     0     0    ' H 2' ' '  \n")
    f.write("  11    0     42     0     0     0    'Na 1' 'data/NaIkas.tl'  \n")
    f.write("  11    1     1      1     0     0    'Na 2' '' \n")
    f.write("  12    0     96     0     0     0    'Mg 1' 'data/Mg1kas_F_ccc.sy'  \n")
    f.write("  12    1     29     0     0     0    'Mg 2' 'data/Mg2kas_F_ccc.sy'  \n")
    f.write("  12    2     1      1     0     0    'Mg 3' ' '  \n")
    f.write("  19    0     31     0     0     0    'K  1' 'data/KIkas.tl'  \n")
    f.write("  19    1     1      1     0     0    'K  2' ''  \n")
    f.write("  20    0     66     0     0     0    'Ca 1' 'data/Ca1kas_F_zat.sy'  \n")
    f.write("  20    1     24     0     0     0    'Ca 2' 'data/Ca2kas_F_zat.sy'  \n")
    f.write("  20    2     1      1     0     0    'Ca 3' ' '  \n")
    f.write("   0    0     0     -1     0     0    '    ' ' '  \n")
    f.write("* \n")
    f.write("* end \n")
  elif atom == 'test': # generic set used in Allende Prieto+ (2018) A&A paper
    f.write("* ../data for ions  \n")
    f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat'  \n")
    f.write("   0    0     3      0 \n")
    f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
    f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
    f.write("  26    0     49     0     0     0    'Fe 1' 'data/tlusty_fe1_topmod.dat'  \n")
    f.write("  26    1     41     0     0     0    'Fe 2' 'data/tlusty_fe2_topmod.dat'  \n")
    f.write("  26    2     1      1     0     0    'Fe 3' ' '  \n")
    f.write("  0     0     0     -1     0     0    '    ' ' '  \n")
    f.write("* \n")
    f.write("* end \n")
  elif atom == 'ap18': # generic set used in Allende Prieto+ (2018) A&A paper
    f.write("* ../data for ions  \n")
    f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat'  \n")
    f.write("   0    0     3      0 \n")
    f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
    f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
    f.write("   2    0     14     0     0     0    'He 1' 'data/he1.dat'  \n")
    f.write("   2    1     14     0     0     0    'He 2' 'data/he2.dat '  \n")
    f.write("   2    2     1      1     0     0    'He 3' ' '  \n")
    f.write("   6    0     104    0     0     0    ' C 1' 'data/c1.t'  \n")
    f.write("   6    1     40     0     0     0    ' C 2' 'data/c2.t'  \n")
    f.write("   6    2     1      1     0     0    ' C 3' ' '  \n")
    f.write("   7    0     89     0     0     0    ' N 1' 'data/n1.t'  \n")
    f.write("   7    1     51     0     0     0    ' N 2' 'data/n2.t'  \n")
    f.write("   7    2     1      1     0     0    ' N 3' ' '  \n")
    f.write("   8    0     54     0      0    0    ' O 1' 'data/o1c.t'  \n")
    f.write("   8    1     74     0      0    0    ' O 2' 'data/o2.t'  \n")
    f.write("   8    2     1      1      0    0    ' O 3' ' '  \n")
    f.write("  11    0     32     0     0     0    'Na 1' 'data/na1.t'  \n")
    f.write("  11    1     8      0     0     0    'Na 2' 'data/na2.t'  \n")
    f.write("  11    2     1      1     0     0    'Na 3' ' '  \n")
    f.write("  12    0     71     0     0     0    'Mg 1' 'data/mg1.t'  \n")
    f.write("  12    1     31     0     0     0    'Mg 2' 'data/mg2.t'  \n")
    f.write("  12    2     1      1     0     0    'Mg 3' ' '  \n")
    f.write("  13    0     33     0     0     0    'Al 1' 'data/al1.t'  \n")
    f.write("  13    1     81     0     0     0    'Al 2' 'data/al2.t'  \n")
    f.write("  13    2     1      1     0     0    'Al 3' ' '  \n")
    f.write("  14    0     57     0     0     0    'Si 1' 'data/si1.t'  \n")
    f.write("  14    1     46     0     0     0    'Si 2' 'data/si2.t'  \n")
    f.write("  14    2     1      1     0     0    'Si 3' ' '  \n")
    f.write("  20    0     79     0     0     0    'Ca 1' 'data/ca1.t'  \n")
    f.write("  20    1     32     0     0     0    'Ca 2' 'data/ca2.t'  \n")
    f.write("  20    2     1      1     0     0    'Ca 3' ' '  \n")
    f.write("  26    0     49     0     0     0    'Fe 1' 'data/tlusty_fe1_topmod.dat'  \n")
    f.write("  26    1     41     0     0     0    'Fe 2' 'data/tlusty_fe2_topmod.dat'  \n")
    f.write("  26    2     1      1     0     0    'Fe 3' ' '  \n")
    f.write("  0     0     0     -1     0     0    '    ' ' '  \n")
    f.write("* \n")
    f.write("* end \n")
  else:
    for line in atominfo: f.write(line)
  f.close()


def write8(teff, logg, nd, atmos, atmostype, ofile='fort.8', 
           kvmicro=2.0, kabu=None,ktitle='written by synple/write8'):

  """Writes the model atmosphere for synspec

     MARCS models can be passed in 'Tlusty' (default, after read with 
           read_marcs_models2) or 'Kurucz' format
     Phoenix and Kurucz models are passed to synspec formatted as 'kurucz'

     kvmicro (microturbulence in km/s), kabu (abundance array), and ktitle 
     (informative text to be added to the header)
     are optional and only used when atmostype is 'kurucz'
   

  """

  f = open(ofile,'w')

  if atmostype == 'tlusty':

    if ('n' in atmos.dtype.names):  # 4th column is number density n
      if ('pop' in atmos.dtype.names):   # explicit (usually NLTE) populations
        numpop = len(atmos['pop'][0]) 
        sformat = '  %f %e %e %e'
        i = 5
        for entry in atmos['pop'][0]: 
           sformat = sformat + ' %e'
           if i % 6 == 0: sformat = sformat + '  \n'
           i = i + 1
        sformat = sformat + ' \n' 
        f.write(" "+str(nd)+" "+str(-(4+numpop))+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i]]
          for j in range(numpop):
            sdata.append(atmos['pop'][i][j])
          f.write( sformat % tuple(sdata) )                 
      elif ('dep' in atmos.dtype.names): # NLTE departure coefficients
        numpop = len(atmos['dep'][0]) 
        sformat = '  %f %e %e %e'
        i = 5
        for entry in atmos['dep'][0]: 
           sformat = sformat + ' %e'
           if i % 6 == 0: sformat = sformat + '  \n'
           i = i + 1
        sformat = sformat + ' \n' 
        f.write(" "+str(nd)+" "+str(-(4+numpop))+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i]]
          for j in range(numpop):
            sdata.append(atmos['dep'][i][j])
          f.write( sformat % tuple(sdata) )         
      else:                              # LTE
        f.write(" "+str(nd)+" "+str(-4)+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          f.write( '%f %e %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['n'][i] ) )  

    else: # n not included, only t, ne and rho given
      if ('pop' in atmos.dtype.names):   # explicit (usually NLTE) populations
        numpop = len(atmos['pop'][0]) 
        sformat = '  %f %e %e '
        i = 5
        for entry in atmos['pop'][0]: 
           sformat = sformat + ' %e'
           if i % 6 == 0: sformat = sformat + '  \n'
           i = i + 1
        sformat = sformat + ' \n' 
        f.write(" "+str(nd)+" "+str(3+numpop)+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i] ]
          for j in range(numpop):
            sdata.append(atmos['pop'][i][j])
          f.write( sformat % tuple(sdata) )                 
      elif ('dep' in atmos.dtype.names): # NLTE departure coefficients
        numpop = len(atmos['dep'][0]) 
        sformat = '  %f %e %e '
        i = 5
        for entry in atmos['dep'][0]: 
           sformat = sformat + ' %e'
           if i % 6 == 0: sformat = sformat + '  \n'
           i = i + 1
        format = sformat + ' \n' 
        f.write(" "+str(nd)+" "+str(3+numpop)+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          sdata = [atmos['t'][i], atmos['ne'][i], atmos['rho'][i] ]
          for j in range(numpop):
            sdata.append(atmos['dep'][i][j])
          f.write( sformat % tuple(sdata) )         
      else:                              # LTE
        f.write(" "+str(nd)+" "+str(3)+"\n")
        for i in range(nd):
          f.write(' %e ' % atmos['dm'][i])
          if (i+1) % 5 == 0: f.write('\n')
        if (i+1) % 5 != 0: f.write('\n')
        for i in range(nd):
          f.write( '%f %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i] ) )  

  else:

    if atmostype == 'marcs':
      f.write(" "+str(nd)+" "+str(-4)+"\n")
      for i in range(nd):
        f.write(' %e ' % atmos['dm'][i])
      f.write("\n")
      for i in range(nd):
        f.write( '%f %e %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['rho'][i]/atmos['mmw'][i]/1.67333e-24 + atmos['ne'][i] ) )

    else: #Kurucz format
        f.write( 'TEFF %7.0f  GRAVITY %7.5f  LTE \n' % (teff, logg) )
        if kabu is None:
            for i in range(21): f.write('\n')
        else:

            abu2 = np.array(kabu,dtype=float)
            ntotalovernh = np.sum(abu2)
            abu2 = abu2 / ntotalovernh
            abu2[2:] = np.log10(abu2[2:])
   
            f.write( 'TITLE %30s\n' % (ktitle) )
            opflags='1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0' 
            f.write( ' OPACITY IFOP %s\n' % opflags)
            f.write( ' CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00\n')
            f.write ('ABUNDANCE SCALE  %8.5f ABUNDANCE CHANGE 1%8.5f 2 %8.5f\n' % (1.00000,abu2[0],abu2[1]) )
            for i in np.arange(16):
                dd = np.zeros(12)
                sixs = np.arange(6, dtype=int)
                dd[np.array(sixs * 2, dtype=int)] = 3 + 6*i + sixs
                dd[np.array(sixs * 2 + 1, dtype=int)] = np.array(abu2,dtype=float)[np.array( 2 + 6*i + sixs, dtype=int)]
                f.write( ' ABUNDANCE CHANGE %2i %6.2f %2i %6.2f %2i %6.2f %2i %6.2f %2i %6.2f %2i %6.2f\n' % tuple(dd) )
#endfor

            dd = np.zeros(2)
            dd[0] = 3 + 6*16
            dd[1] = abu2[2 + 6*16]
            f.write(' ABUNDANCE CHANGE %2i %6.2f\n' % tuple(dd) )

        if kvmicro is None:
            f.write( 'READ DECK6%3i RHOX,T,P,XNE\n' % nd )
            for i in range(nd):
                f.write( '%e %f %e %e\n' % (atmos['dm'][i], atmos['t'][i], atmos['p'][i], atmos['ne'][i]) )

        else:
            f.write( 'READ DECK6%3i RHOX,T,P,XNE,ZERO,ZERO,VTURB\n' % nd )
            for i in range(nd): 
                f.write( '%e %f %e %e %8.2e %8.2e %8.2e\n' % (atmos['dm'][i], atmos['t'][i], atmos['p'][i], atmos['ne'][i], 0.0, 0.0, kvmicro*1e5) )
      
  f.close()

  return()


def write_innterpol_model(filename,feh,cfe,afe,teff,logg,array,
                          vmicro=2.0,title='iNNterpolated model'): 

   """
   write a model atmosphere from iNNterpol to disk 

   Parameters
   ----------
   filename: str
     name of the file to contail the model atmosphere

   feh: float
     metallicity [Fe/H]    
     
   cfe: float
     carbon to iron abundance [C/Fe]

   afe: float
     alpha-element to iron abundance [alpha/Fe]

   teff: float
     effective temperature (K)

   logg: float
     surface gravity (log10 of the surface gravity in cm/s2)

   array: numpy array of floats
     atmospheric structure with as many rows as layers and 4 columns
     corresponding to the log10 of column mass, temperature, gas pressure and 
     electron density. This array is exactly what iNNterpol returns.
     see https://github.com/cwestend/iNNterpol

   vmicro: float
     microturbulence (km/s)
     (default is 2.)

   title: str
     explanatory text to be added to the header

   Returns
   -------    

   no output, it creates a text file with a Kurucz-formatted model atmosphere

   """

   #we start by repeating the top layer, since it was dropped for 
   #building the iNNterpol NNs and it will again be dropped by Synspec 

   array = np.vstack( (array[0,:], array) )

   nd = len(array[:,0])

   atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')})

   array = 10.**array

   atmos['dm'] = array[:,0]
   atmos['t'] = array[:,1]
   atmos['p'] = array[:,2]
   atmos['ne'] = array[:,3]

   symbol, mass, abu = elements()
   if np.abs(feh) > 0.: abu[2:] = abu[2:] * 10.**feh
   if np.abs(cfe) > 0.: abu[5] = abu[5] * 10.**cfe
   if np.abs(afe) > 0.: 
       #O, Ne, Mg, Si, S, Ca, and Ti.
       for item in [8,10,12,14,16,20,22]:
           abu[item-1] = abu[item-1] * 10.**afe


   write8(teff,logg,nd,atmos,'kurucz',ofile=filename,
          kvmicro=vmicro,kabu=abu,
          ktitle=title)
    
   return()



def read10(file='fort.10'):
#read output synspec file with specific intensities
#returning wavelengths and intensities in the file

  f = open(file,'r')
  x =  []
  y = [] 
  for line in f:
    entries = line.split()
    if len(entries) == 2: 
     x.append(entries[0]) 
    else:
     y.extend(entries)
 

  iwave = np.array(x,dtype=float)
  inte = np.array(y,dtype=float).reshape(len(x),int(len(y)/len(x)))

  return(iwave,inte)
  
  

def create_links(linelist):
#create soft links for line lists

  for i in range(len(linelist)):
    file = linelist[i]
    binaryfile = linelist[i][:-2]+'11'
    if os.path.isfile(binaryfile): file = binaryfile
    if i == 0: os.symlink(file,'fort.19')
    else: os.symlink(file,'fort.'+str(20-1+i))

  return()

def cleanup_fort():
#cleanup all fort* files


  files = os.listdir('.')
  for entry in files: 
    if os.path.islink(entry) and entry.startswith('fort'): os.unlink(entry)
    if os.path.isfile(entry) and entry.startswith('fort'): os.remove(entry)

  return()


def read_kurucz_model(modelfile):
  
  """Reads a Kurucz model atmospheres
  
  Parameters
  ----------
  modelfile: str
      file name  
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
  """

  f = open(modelfile,'r')
  line = f.readline()
  entries = line.split()
  assert (entries[0] == 'TEFF' and entries[2] == 'GRAVITY'), 'Cannot find Teff and logg in the file header'
  teff = float(entries[1])
  logg = float(entries[3])

  while entries[0] != 'ABUNDANCE':  
    line = f.readline()
    entries = line.split()

  abu = []

  if entries[1] == 'SCALE': 
    scale = float(entries[2])
    

  while entries[0] == 'ABUNDANCE':
    i = 0
    for word in entries: 
      if (word == 'CHANGE'): w = i
      i = i + 1 
    for i in range(int((len(entries)-w-1)/2)):
      z = int(entries[w+1+2*i])
      if (z == 1): nhntot = float(entries[w+2+2*i])
      if (z < 3): abu.append(float(entries[w+2+2*i]) / nhntot) 
      else: abu.append(scale*10.**(float(entries[w+2+2*i])) / nhntot)

    line = f.readline()
    entries = line.split() 

  assert (entries[0] == 'READ'), 'I cannot find the header of the atmospheric table in the input Kurucz model'

  nd = int(entries[2]) - 1
  line = f.readline()
  entries = line.split()
  line = f.readline()
  entries = line.split()
  vmicro = float(entries[6])/1e5

  dm = [ float(entries[0]) ]
  t = [ float(entries[1]) ]
  p = [ float(entries[2]) ]
  ne = [ float(entries[3]) ] 
 
  good = True
  for i in range(nd-1):
    line = f.readline()
    entries = line.split()
    for char in entries:
      if char.count('E') > 1:
        good = False

    if not good: break

    dm.append(float(entries[0]))
    t.append(float(entries[1]))
    p.append(float(entries[2]))
    ne.append( float(entries[3]))
 
  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 

  if good:
    atmos['dm'] = dm
    atmos['t'] = t
    atmos['p'] = p
    atmos['ne'] = ne

  else:
    return(None,None,None,None,None,None)

  return (teff,logg,vmicro,abu,nd,atmos)


def read_marcs_model(modelfile):
  
  """Reads a MARCS model atmospheres
  
  Parameters
  ----------
  modelfile: str
      file name. It can be a gzipped (.gz) file
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
  """  

  if modelfile[-3:] == '.gz':
    f = gzip.open(modelfile,'rt')
  else:
    f = open(modelfile,'r')
  line = f.readline()
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Teff'), 'Cannot find Teff in the file header'
  teff = float(entries[0])
  line = f.readline()
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Surface' and entries[2] == 'gravity'), 'Cannot find logg in the file header'
  logg = np.log10(float(entries[0]))
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Microturbulence'), 'Cannot find vmicro in the file header'
  vmicro = float(entries[0])

  while entries[0] != 'Logarithmic':  
    line = f.readline()
    entries = line.split()

  abu = []
  line = f.readline()
  entries = line.split()

  i = 0
  while entries[1] != 'Number':
    for word in entries: 
      abu.append( 10.**(float(word)-12.0) )
      i = i + 1 
    line = f.readline()
    entries = line.split() 

  if i < 99: 
    for j in range(99-i):
      abu.append(1e-111)
      i = i + 1

  nd = int(entries[0])
  line = f.readline()
  entries = line.split()

  assert (entries[0] == 'Model'), 'I cannot find the header of the atmospheric table in the input MARCS model'

  line = f.readline()
  line = f.readline()
  entries = line.split()

  t = [ float(entries[4]) ]
  p = [ float(entries[6]) ]
  ne = [ float(entries[5]) / bolk / float(entries[4]) ] 

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    t.append(  float(entries[4]))
    p.append(  float(entries[6]))
    ne.append( float(entries[5]) / bolk / float(entries[4]))

  line = f.readline()
  line = f.readline()
  entries = line.split()

  dm = [ float(entries[-1]) ]

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    dm.append(  float(entries[-1]))

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 
  atmos['dm'] = dm
  atmos['t'] = t
  atmos['p'] = p
  atmos['ne'] = ne

  return (teff,logg,vmicro,abu,nd,atmos)

def read_marcs_model2(modelfile):
  
  """Reads a MARCS model atmospheres. 
  While read_marcs_model returns T, Pg and Ne in the structure 'atmos'
  read_marcs_model2 returns T, rho, mmw, and Ne.
  
  Parameters
  ----------
  modelfile: str
      file name. It can be a gzipped (.gz) file
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, density, 
      mean molecular weight and electron number density  
  
  """  

  if modelfile[-3:] == '.gz':
    f = gzip.open(modelfile,'rt')
  else:
    f = open(modelfile,'r')
  line = f.readline()
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Teff'), 'Cannot find Teff in the file header'
  teff = float(entries[0])
  line = f.readline()
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Surface' and entries[2] == 'gravity'), 'Cannot find logg in the file header'
  logg = np.log10(float(entries[0]))
  line = f.readline()
  entries = line.split()
  assert (entries[1] == 'Microturbulence'), 'Cannot find vmicro in the file header'
  vmicro = float(entries[0])

  while entries[0] != 'Logarithmic':  
    line = f.readline()
    entries = line.split()

  abu = []
  line = f.readline()
  entries = line.split()

  i = 0
  while entries[1] != 'Number':
    for word in entries: 
      abu.append( 10.**(float(word)-12.0) )
      i = i + 1 
    line = f.readline()
    entries = line.split() 

  if i < 99: 
    for j in range(99-i):
      abu.append(1e-111)
      i = i + 1

  nd = int(entries[0])
  line = f.readline()
  entries = line.split()

  assert (entries[0] == 'Model'), 'I cannot find the header of the atmospheric table in the input MARCS model'

  line = f.readline()
  line = f.readline()
  entries = line.split()

  t = [ float(entries[4]) ]
  p = [ float(entries[6]) ]
  ne = [ float(entries[5]) / bolk / float(entries[4]) ] 

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    t.append(  float(entries[4]))
    p.append(  float(entries[6]))
    ne.append( float(entries[5]) / bolk / float(entries[4]))

  line = f.readline()
  line = f.readline()
  entries = line.split()

  rho = [ float(entries[3]) ]
  dm = [ float(entries[-1]) ]
  mmw = [ float(entries[4]) ]

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    rho.append( float(entries[3]))
    dm.append(  float(entries[-1]))
    mmw.append(  float(entries[4]))

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'rho','mmw','ne'),
                          'formats':('f', 'f', 'f','f','f')}) 
  atmos['dm'] = dm
  atmos['t'] = t
  atmos['rho'] = rho
  atmos['mmw'] = mmw
  atmos['ne'] = ne

  return (teff,logg,vmicro,abu,nd,atmos)

def read_tlusty_model(modelfile,startdir=None):
  
  """Reads a Tlusty model atmosphere. 

  Parameters
  ----------
  modelfile: str
      file name (.7, .8, or .22). It will look for the complementary .5 file to read
      the abundances and the micro (when specified in the non-std. parameter file)

  startdir: str
      directory where the calculations are initiated. The code will look at that
      location to find the tlusty model atom directory and the non-std. parameter
      file when a relative path is provided
      (default is None, indicating it is the current working directory)
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s), by default 1.0 unless set with the parameter
      VTB in the non-std. parameter file specified in the .5 file
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, density
      (other variables that may be included, e.g. populations for NLTE models, 
      are ignored). 

  """  

  assert ((modelfile[-2:] == ".8") | (modelfile[-2:] == ".7") | (modelfile[-3:] == ".22")), 'Tlusty models should end in .7, .8, or .22'
  if modelfile[-2] == ".":
    madaffile = modelfile[:-1]+"5"
  else:
    madaffile = modelfile[:-2]+"5"    
  assert (os.path.isfile(madaffile)),'Tlusty model atmosphere file '+modelfile+' should come with an associated .5 file'

  if startdir is None: startdir = os.getcwd()

  #read the madaf (.5) file
  teff, logg, vmicro, abu = read_madaf(madaffile,startdir=startdir)
  
  #now the structure (.8) file
  f = open(modelfile,'r')
  line = f.readline()
  entries = line.split()
  nd = int(entries[0])
  numpar = int(entries[1])
  if (numpar < 0): 
    numpop = abs(numpar) - 4 
  else:
    numpop = numpar - 3

  assert (len(entries) == 2), 'There are more than two numbers in the first line of the model atmosphere'

  dm = read_multiline_fltarray(f,nd)
  atm = read_multiline_fltarray(f,nd*abs(numpar))
  f.close()

  atm = np.reshape(atm, (nd,abs(numpar)) )

  if (numpar < 0):  # 4th column is number density n
    if (numpop > 0): # explicit (usually NLTE) populations
      if modelfile[-2] == ".":  # NLTE populations or departure coefficients
        tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f'), ('pop', 'f', (numpop))])
      else: 
        tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f'), ('dep', 'f', (numpop))])
    else:
      tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('n','f')])  
  else:
    if (numpop > 0):
      if modelfile[-2] == ".": # NLTE populations or departure coefficients
        tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('pop', 'f', (numpop))])
      else:
        tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f'), ('dep', 'f', (numpop))])
    else:
      tp = np.dtype([('dm', 'f'), ('t','f'), ('ne','f'), ('rho','f') ])

  atmos = np.zeros(nd, dtype=tp)

  atmos['dm'] = dm
  atmos['t'] = atm [:,0]
  atmos['ne'] = atm [:,1]
  atmos['rho'] = atm [:,2]
  if (numpar < 0): atmos['n'] = atm [:,3]
  if (numpop > 0): 
    if modelfile[-2] == ".":
      atmos['pop'] = atm [:,3:]
    else:
      atmos['dep'] = atm [:,3:]

  return (teff,logg,vmicro,list(abu.values()),nd,atmos)

def read_madaf(madaffile,startdir=None):
  
  """Reads a Tlusty MADAF (.5) file with parameters and abundances. 

  Parameters
  ----------
  madaffile: str
      file name (.5) including the abundances and the micro (when specified in 
      the non-std. parameter file)

  startdir: str
      directory where the calculations are initiated. The code will look at that
      location to find the tlusty model atom directory and the non-std. parameter
      file when a relative path is provided
      (default is None, indicating it is the current working directory)
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s), by default 1.0 unless set with the parameter
      VTB in the non-std. parameter file specified in the .5 file
  abu : dictionary 
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es) -- the keys are the elemental symbols

  """  

  if startdir is None: startdir = os.getcwd()
  
  assert (os.path.isfile(madaffile)),'The input madaf file '+madaffile+' is not good'

  #we start reading the .5
  f = open(madaffile,'r')
  line = f.readline()
  entries = line.split()
  teff = float(entries[0])
  logg = float(entries[1])
  line = f.readline()
  line = f.readline()
  entries = line.split()
  nonstdfile = entries[0][1:-1]

  nonstdfile0 = nonstdfile
  if nonstdfile != '':
    if not os.path.isabs(nonstdfile): 
      mf = os.path.join(startdir,nonstdfile)
      if os.path.isfile(mf): 
        nonstdfile = mf
      else:
        mf = os.path.join(modeldir,nonstdfile)
        nonstdfile = mf

    assert (os.path.exists(nonstdfile)), 'The non-std parameter file indicated in the tlusty model, '+nonstdfile0+', is not present' 

  nonstd={}
  if nonstdfile != '':
    assert (os.path.isfile(nonstdfile)),'Tlusty model atmosphere file '+modelfile+' invokes non-std parameter file, '+nonstdfile+' which is not present'


    ns = open(nonstdfile,'r')
    nonstdarr = ns.readlines()
    ns.close()
    for entry in nonstdarr:
      entries = entry.replace(',\n','').split(',')
      for piece in entries:
        sides = piece.split('=')
        nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')

    print('Tlusty nonstd params=',nonstd)

  #the micro might be encoded as VTB in the nonstdfile!!
  #this is a temporary patch, but need to parse that file
  vmicro = 1.0
  if 'VTB' in nonstd: vmicro = float(nonstd['VTB'])

  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  entries = line.split()
  natoms = int(entries[0])

  symbol, mass, sol = elements() 
  abu = dict()

  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  for i in range(natoms):
    entries = line.split()
    #abu.append( float(entries[1]) )
    abu[symbol[i]] = float(entries[1])
    line = f.readline()

  if i < 98: 
    for j in range(98-i):
      #abu.append(1e-111)
      abu[symbol[i+1]] = 1e-111
      i = i + 1

  f.close()

  return(teff, logg, vmicro, abu)

def read_tlusty_extras(modelfile,startdir=None):
  
  """Identifies and reads the non-std parameter file and its content, finds out the 
     number of parameters in the model, whether the file contains populations or departure
     coefficients, and the name of the data directory for Tlusty 
     model atmospheres. 

  Parameters
  ----------
  modelfile: str
      file name (.8, .7 or .22). It will look for the complementary .5 file to read
      the abundances and other information

  startdir: str
      directory where the calculations are initiated. The code will look at that
      location to find the tlusty model atom directory and the non-std. parameter
      file when a relative path is provided
      (default is None, indicating it is the current working directory)
  
  
  Returns
  -------

  madaffile: str
       model atom data and abundance file (.5 Tlusty file)

  nonstdfile: str
       non-std parameter file 

  nonstd: dict
       content of the non-std parameter file

  numpar: int
       number of parameters (can be negative when the model includes number density)

  datadir: str
       name of the model atom directory

  inlte: int
       0 when the populations are to be computed internally by synspec (LTE)
       1 the Tlusty model contains populations
      -1 the Tlusty model contains departure coefficients

  atommode: list
       mode for each of the atoms included. The code indicates
       0= not considered
       1= implicit (no cont. opacity)
       2= explicit  (see synspec man.)
       4= semi-explicit (see synspec man.)
       5= quasi-explicit  (see synspec. man)

  atominfo: list
       all the lines in the file that provide info on the model atoms used
  
  """  

  assert ((modelfile[-2:] == ".8") | (modelfile[-2:] == ".7") | (modelfile[-3:] == ".22")), 'Tlusty models should end in .7, .8, or .22'
  if modelfile[-2] == ".":
    madaffile = modelfile[:-1]+"5"
  else:
    madaffile = modelfile[:-2]+"5"    
  assert (os.path.isfile(madaffile)),'Tlusty model atmosphere file '+modelfile+' should come with an associated .5 file'

  if startdir is None: startdir = os.getcwd()

  #we start reading the .5
  f = open(madaffile,'r')
  line = f.readline()
  line = f.readline()
  line = f.readline()
  entries = line.split()
  nonstdfile = entries[0][1:-1]

  nonstdfile0 = nonstdfile  
  if nonstdfile != '':
    if not os.path.isabs(nonstdfile): 
      md = startdir
      mf = os.path.join(md,nonstdfile)
      if os.path.isfile(mf): 
        nonstdfile = mf
      else:
        md = modeldir
        mf = os.path.join(md,nonstdfile)
        nonstdfile = mf
    else:
      md = os.path.split(nonstdfile)[0]

    assert (os.path.exists(nonstdfile)), 'The non-std parameter file indicated in the tlusty model, '+nonstdfile0+', is not present' 


  nonstd={}
  if nonstdfile != '':
    assert (os.path.isfile(nonstdfile)),'Tlusty model atmosphere file '+modelfile+' invokes non-std parameter file, '+nonstdfile+' which is not present'


    ns = open(nonstdfile,'r')
    nonstdarr = ns.readlines()
    ns.close()
    for entry in nonstdarr:
      entries = entry.replace(',\n','').split(',')
      for piece in entries:
        sides = piece.split('=')
        nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')

  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  entries = line.split()
  natoms = int(entries[0])

  line = f.readline()
  while line[0] == '*':
    line = f.readline()
  atommode = []
  for i in range(natoms):
    entries = line.split()
    atommode.append(int(entries[0]))
    line = f.readline()

  atominfo = []
  #keep reading until you find 'dat' to identify data directory 
  #line = f.readline()
  while True: 
    atominfo.append(line)
    if '.dat' in line: break
    line = f.readline()

  entries = line.split()
  cadena = entries[-1][1:-1]
  datadir, file = os.path.split(cadena)


  datadir0 = datadir
  if datadir != '':
    if not os.path.isabs(datadir): 
      mf = os.path.join(md,datadir)
      if os.path.exists(mf): 
        datadir = mf
      else:
        mf = os.path.join(synpledir,'data',datadir)
        datadir = mf

  #continue reading the rest of the file into atominfo
  line = f.readline()
  while True:
    if line == '': break
    atominfo.append(line)
    line = f.readline()

    assert (os.path.exists(datadir)), 'The datadir indicated in the tlusty model, '+datadir0+', is not present' 


  f.close()

  #now the .8
  f = open(modelfile,'r')
  line = f.readline()
  entries = line.split()
  nd = int(entries[0])
  numpar = int(entries[1])
  if abs(numpar) > 4: 
    inlte = 1 
  else: 
    inlte = 0

  if (modelfile[-3:] == ".22"): inlte = -1

  f.close()

  return (madaffile, nonstdfile, nonstd, numpar, datadir, inlte, atommode, atominfo)


def read_phoenix_model(modelfile):

  """Reads a FITS Phoenix model atmospheres
  
  Parameters
  ----------
  modelfile: str
      file name  
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
  """  

  from astropy.io import fits

  h = fits.open(modelfile)[0].header
  f = fits.open(modelfile)[1].data

  nd = len(f['temp'])

  teff = float(h['PHXTEFF'])
  logg = float(h['PHXLOGG'])
  vmicro = float(h['PHXXI_L'])

  m_h = float(h['PHXM_H'])
  alpha = float(h['PHXALPHA'])
  
  symbol, mass,sol = elements(reference='husser') 
  abu = sol 
  z_metals = np.arange(97,dtype=int) + 3
  z_alphas = np.array([8,10,12,14,16,20,22],dtype=int)
  for i in range(len(z_metals)): abu[z_metals[i] - 1] = abu[z_metals[i] - 1] + m_h
  for i in range(len(z_alphas)): abu[z_alphas[i] - 1] = abu[z_alphas[i] - 1] + alpha
  

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 

  atmos['dm'] = f['pgas'] / 10.**logg
  atmos['t'] = f['temp']
  atmos['p'] = f['pgas']
  atmos['ne'] = f['pe']/ bolk / f['temp']

  return (teff,logg,vmicro,abu,nd,atmos)


def read_phoenix_text_model(modelfile):
  
  
  """Reads a plain-text Phoenix model atmospheres
  
  Parameters
  ----------
  modelfile: str
      file name  
  
  Returns
  -------

  teff : float
      effective temperature (K)
  logg : float
      log10 of the surface gravity (cm s-2)
  vmicro : float
      microturbulence velocity (km/s)
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)
  nd: int
      number of depths (layers) of the model
  atmos: numpy structured array
      array with the run with depth of column mass, temperature, gas pressure 
      and electron density  
  
  """  


  f = open(modelfile,'r')
  line = f.readline()
  while line[0:4] != " no.":
    line = f.readline()
  entries = line.split()
  nd = int(entries[5])
  print('nd=',nd)
  while line[0:14] != " model:   teff":
    line = f.readline()
  entries = line.split()
  teff = float(entries[3])
  print('teff=',teff)
  line = f.readline()
  line = f.readline()
  entries = line.split()
  assert (entries[0] == 'log(g):' and entries[2] == '[cm/s**2]'), 'Cannot find logg in the file header'
  logg = float(entries[1])
  print('logg=',logg)
  line = f.readline()
  while line[0:22] !=  "  Element abundances :":  
    line = f.readline()


  symbol,mass,sol = elements()

  sy = []
  ab = []

  while line[0:29] !=  "  Element abundances relative":  
    line = f.readline()
    #print(line)
    if line[0:9] == ' element:':
      entries = line.split()
      for word in entries[1:]: sy.append(word)
    if line[0:11] == ' abundance:':
      entries = line.split()
      for word in entries[1:]: ab.append(word)

  assert (len(sy) == len(ab)), 'different elements in arrays sy (elemental symbols) and ab (abundances)'

  abu = np.ones(99)*1e-99
  i = 0
  for item in sy:
    try:
      index = symbol.index(item)
      abu[index] =  10.**(float(ab[i])-12.) 
    except ValueError:
      print("the symbol ",item," is not recognized as a valid element")
    i = i + 1

  print('abu=',abu)

  while line[0:72] !=  "   l        tstd temperature        pgas          pe     density      mu":  
    line = f.readline()

  line = f.readline()
  entries = line.split()

  t = [ float(entries[2].replace('D','E')) ]
  p = [ float(entries[3].replace('D','E')) ]
  ne = [ float(entries[4].replace('D','E')) / bolk / float(entries[2].replace('D','E')) ] 
  dm = [ float(entries[3].replace('D','E')) / 10.**logg ] #assuming hydrostatic equil. and negliglible radiation and turb. pressure

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    t.append(  float(entries[2].replace('D','E')))
    p.append(  float(entries[3].replace('D','E')))
    ne.append( float(entries[4].replace('D','E')) / bolk / float(entries[2]))
    dm.append ( float(entries[3].replace('D','E')) / 10.**logg )

  vmicro = 1.0
  while (line[0:6] != " greli"):
    line = f.readline()
    if line == '':
        print('Cannot find a value for vmicro (vturb) in the model atmosphere file ',modelfile)
        break
  
  if line != '':
    entries = line.split()
    vmicro = float(entries[5])

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 
  atmos['dm'] = dm
  atmos['t'] = t
  atmos['p'] = p
  atmos['ne'] = ne

  return (teff,logg,vmicro,abu,nd,atmos)


def read_multiline_fltarray(fhandle,arrlen):

  """Reads a float array that spans one or multiple lines in a file
  """

  ndata = 0
  arr = []
  while ndata < arrlen:
    line = fhandle.readline()
    line = line.replace('D','E')
    line = line.replace('d','e')
    entries = line.split()
    for val in entries: arr.append( float(val) )
    ndata = len(arr)

  return (arr)

def read_opt(filename):
    """Reads a synspec/tlusty opacity table (see routine polyopt)

    Parameters
    ----------
       filename: string
                 name of the file containing the opacity table


   Returns
   -------
       lt: float array
            log10(T) for opacity grid (T in K)
       lrho:  float array
            log10(rho)  for opacity grid (rho in gr/cm3)
       lambda0: numpy float array
            vacuum wavelength (Angstroms)
       lopa: numpy float array
            log10(opacity/rho), where opacity is in cm-1, and opacity/rho in cm2/gr
       abu_eos: array of floats (99 elements) 
            chemical abundances relative to hydrogen (N(X)/N(H)) adopted for the 
            equation of state (elements not explicit in the table header are assigned
            zero abundance)
       abu_opa: array of floats (99 elements)
            chemical abundances relative to hydrogen (N(X)/N(H)) adopted for the 
            opacities in the table  (elements not explicit in the table header are assigned
            zero abundance)

    """

    symbol, mass, sol = elements()
    abu_eos = np.zeros(99)
    abu_opa = np.zeros(99)

    f = open(filename,'r')
    entry = f.readline()
    entry = f.readline()
    entry = f.readline()
    while entry != "\n":
        #print(entry)
        entries = entry.split()
        try:
            index = symbol.index(entries[0])
            abu_eos[index] = float(entries[1])
            abu_opa[index] = float(entries[2])
        except ValueError:
            print("the line ",entry," is not recognized as describing the abundance of an element")
        entry = f.readline()

    while entry[:13] != "number of fre":
        entry = f.readline()

    entry = f.readline()
    entries = entry.split()
    nlambda = np.int64(entries[0])
    lambda0 = np.zeros(nlambda, dtype = float)
    nt = int(entries[1])
    nrho = int(entries[2])
    lopa = np.zeros( (nlambda,nrho,nt), dtype=float)
    print('nrho, nt, nlambda = ', nrho, nt, nlambda)
    

    while entry[:13] != "log temperatu":
        entry = f.readline()
    
    lt = read_multiline_fltarray(f,nt)

    while entry[:13] != "log densities":
        entry = f.readline()
    
    lrho = read_multiline_fltarray(f,nrho)

    entry = f.readline()
    i = 0
    while entry != "":
        entry = f.readline()
        if  entry[1:13] == "*** frequenc":
            entries = entry.split()
            assert (int(entries[-2]) == i+1), 'The frequencies do not appear in order: '+entries[-2]+' != '+str(i+1)
            lambda0[i] = float(entries[-1])
            entry = f.readline()
            arr = read_multiline_fltarray(f,nt*nrho)
            lopa[i,:,:] = np.reshape(arr, (nrho,nt))
            i += 1

    return(lrho,lt,lambda0,lopa,abu_eos,abu_opa)
    
def read_copt(filename,nrho,nt):
    """Reads a synspec/tlusty continuum opacity table (fort.26)

    Parameters
    ----------
       filename: string
                 name of the file containing the opacity table,
                 this text file only contains opacity and the information
                 about the densities/temperatures/abundances needs to
                 be extracted from its associated full-opacity table

       nrho: int
                 number of density points in the table
       nt:   int 
                 number of temperature points in the table

   Returns
   -------
       wave: 3D numpy float array (nwave, nrho, nt)
            wavelenght in Angstroms
            
       lopa: 3D numpy float array (nwave, nrho, nt)
            log10(opacity/rho), where opacity is in cm-1, and opacity/rho in cm2/gr
    """
            
    c = np.loadtxt(filename)
    wave = np.transpose( np.reshape(c[:,0], (nt,nrho,int(c.shape[0]/nt/nrho)) ) )
    lopa = np.transpose( np.reshape(c[:,1], (nt,nrho,int(c.shape[0]/nt/nrho)) ) )

    return(wave,lopa)
    
    
def interp_spl(xout, x, y):

  """Interpolates in 1D using cubic splines

  Parameters
  ----------
     x: numpy array or list
        input abscissae
     y: numpy array or list
        input ordinates 
     xout: numpy array or list
        array of abscissae to interpolate to

   Returns
   -------
     yout: numpy array or list
        array of interpolated values

  """

  tck = interpolate.splrep(x, y, s=0)
  yout = interpolate.splev(xout, tck, der=0)

  return(yout)



def interp_spl2(x0, x, y):
    """
    Interpolate a 1-D function using cubic splines.
      x0 : a float or an 1d-array
      x : (N,) array_like
          A 1-D array of real/complex values.
      y : (N,) array_like
          A 1-D array of real values. The length of y along the
          interpolation axis must be equal to the length of x.

    Implement a trick to generate at first step the cholesky matrice L of
    the tridiagonal matrice A (thus L is a bidiagonal matrice that
    can be solved in two distinct loops).

    additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf 
    
    code from
    https://stackoverflow.com/questions/31543775/how-to-perform-cubic-spline-interpolation-in-python
    """
    x = np.asfarray(x)
    y = np.asfarray(y)

    # remove non finite values
    # indexes = np.isfinite(x)
    # x = x[indexes]
    # y = y[indexes]

    # check if sorted
    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    # allocate buffer matrices
    Li = np.empty(size)
    Li_1 = np.empty(size-1)
    z = np.empty(size)

    # fill diagonals Li and Li-1 and solve [L][y] = [B]
    Li[0] = np.sqrt(2*xdiff[0])
    Li_1[0] = 0.0
    B0 = 0.0 # natural boundary
    z[0] = B0 / Li[0]

    for i in range(1, size-1, 1):
        Li_1[i] = xdiff[i-1] / Li[i-1]
        Li[i] = np.sqrt(2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])
        Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
        z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    i = size - 1
    Li_1[i-1] = xdiff[-1] / Li[i-1]
    Li[i] = np.sqrt(2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])
    Bi = 0.0 # natural boundary
    z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    # solve [L.T][x] = [y]
    i = size-1
    z[i] = z[i] / Li[i]
    for i in range(size-2, -1, -1):
        z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]

    # find index
    index = x.searchsorted(x0)
    np.clip(index, 1, size-1, index)

    xi1, xi0 = x[index], x[index-1]
    yi1, yi0 = y[index], y[index-1]
    zi1, zi0 = z[index], z[index-1]
    hi1 = xi1 - xi0

    # calculate cubic
    f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
         zi1/(6*hi1)*(x0-xi0)**3 + \
         (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
         (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
    return f0



def elements(reference=None):

  
  """Reads the solar elemental abundances
  
  Parameters
  ----------
     reference: string, optional
        set to 'husser; for the abundances adopted for Phoenix models by Huser et al. (2013),
        set to 'basti' for abundances adopted for the BaSTI stellar models (Hidalgo et al. 2018),
        set to 'ags2005' for Asplund et al. (2005) are used -- consistent with
        the MARCS (Gustafsson et al. 2008) models and and Kurucz (Meszaros et al. 2012)
        Kurucz model atmospheres.
        (default, reference=None, will trigger 'ags2005')
        
  Returns
  -------
     symbol: numpy array of str
        element symbols
     mass: numpy array of floats
        atomic masses (elements Z=1-99)
     sol: numpy array of floats
        solar abundances N/N(H)
  
  """

  symbol = [
  'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne', 
  'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca', 
  'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn', 
  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr', 
  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', 
  'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd', 
  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', 
  'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg', 
  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', 
  'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es' ]

  mass = [ 1.00794, 4.00260, 6.941, 9.01218, 10.811, 12.0107, 14.00674, 15.9994,
  18.99840, 20.1797, 22.98977, 24.3050, 26.98154, 28.0855, 30.97376, 
  32.066, 35.4527, 39.948, 39.0983, 40.078, 44.95591, 47.867, 50.9415, 
  51.9961, 54.93805, 55.845, 58.93320, 58.6934, 63.546, 65.39, 69.723, 
  72.61, 74.92160, 78.96, 79.904, 83.80, 85.4678, 87.62, 88.90585, 
  91.224, 92.90638, 95.94, 98., 101.07, 102.90550, 106.42, 107.8682, 
  112.411, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.29, 
  132.90545, 137.327, 138.9055, 140.116, 140.90765, 144.24, 145, 150.36, 
  151.964, 157.25, 158.92534, 162.50, 164.93032, 167.26, 168.93421, 
  173.04, 174.967, 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 
  195.078, 196.96655, 200.59, 204.3833, 207.2, 208.98038, 209., 210., 
  222., 223., 226., 227., 232.0381, 231.03588, 238.0289, 237., 244., 
  243., 247., 247., 251., 252. ]

  if  reference == 'husser':
    #a combination of meteoritic/photospheric abundances from Asplund et al. 2009
    #chosen for the Husser et al. (2013) Phoenix model atmospheres
    sol = [  12.00, 10.93,  3.26,  1.38,  2.79,  8.43,  7.83,  8.69,  4.56,  7.93, 
    6.24,  7.60,  6.45,  7.51,  5.41,  7.12,  5.50,  6.40,  5.08,  6.34, 
    3.15,  4.95,  3.93,  5.64,  5.43,  7.50,  4.99,  6.22,  4.19,  4.56, 
    3.04,  3.65,  2.30,  3.34,  2.54,  3.25,  2.36,  2.87,  2.21,  2.58, 
    1.46,  1.88, -9.99,  1.75,  1.06,  1.65,  1.20,  1.71,  0.76,  2.04, 
    1.01,  2.18,  1.55,  2.24,  1.08,  2.18,  1.10,  1.58,  0.72,  1.42, 
   -9.99,  0.96,  0.52,  1.07,  0.30,  1.10,  0.48,  0.92,  0.10,  0.92, 
    0.10,  0.85, -0.12,  0.65,  0.26,  1.40,  1.38,  1.62,  0.80,  1.17,
    0.77,  2.04,  0.65, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99,  0.06,   
   -9.99, -0.54, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99 ]


  elif reference == 'basti':
  #recommended solar abundances from Lodders (2011) https://ui.adsabs.harvard.edu/abs/2010ASSP...16..379L/abstract
  #except for C, N, O, P, S, K, and Fe for which values are from Caffau et al. (2011)
  #https://ui.adsabs.harvard.edu/abs/2011SoPh..268..255C/abstract
    sol = [ 12.00, 10.925, 3.28, 1.32, 2.81,  8.50,  7.86,  8.76, 4.44, 8.05, 
    6.29,  7.54,  6.46,  7.53,  5.46,  7.16,  5.25,  6.50,  5.11, 6.31, 
    3.07,  4.93,  3.99,  5.65,  5.50,  7.52,  4.90,  6.22,  4.27, 4.65, 
    3.10,  3.59,  2.32,  3.36,  2.56,  3.28,  2.38,  2.90,  2.20, 2.57, 
    1.42,  1.94, -9.99,  1.78,  1.10,  1.67,  1.22,  1.73,  0.78, 2.09, 
    1.03,  2.20,  1.57,  2.27,  1.10,  2.18,  1.19,  1.60,  0.77, 1.47, 
   -9.99,  0.96,  0.53,  1.09,  0.34,  1.14,  0.49,  0.95,  0.14, 0.94, 
    0.11,  0.73, -0.14,  0.67,  0.28,  1.37,  1.36,  1.64,  0.82, 1.19, 
    0.79,  2.06,  0.67, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, 0.08, 
   -9.99, -0.52, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99]

  elif (reference is None or reference == 'ags2005'):
    #Asplund, Grevesse and Sauval (2005), basically the same as 
    #Grevesse N., Asplund M., Sauval A.J. 2007, Space Science Review 130, 205
    sol = [  0.911, 10.93,  1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,  7.84, 
    6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,  5.08,  6.31, 
    3.05,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,  6.23,  4.21,  4.60, 
    2.88,  3.58,  2.29,  3.33,  2.56,  3.28,  2.60,  2.92,  2.21,  2.59, 
    1.42,  1.92, -9.99,  1.84,  1.12,  1.69,  0.94,  1.77,  1.60,  2.00, 
    1.00,  2.19,  1.51,  2.27,  1.07,  2.17,  1.13,  1.58,  0.71,  1.45, 
   -9.99,  1.01,  0.52,  1.12,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08, 
    0.06,  0.88, -0.17,  1.11,  0.23,  1.45,  1.38,  1.64,  1.01,  1.13,
    0.90,  2.00,  0.65, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99,  0.06,   
   -9.99, -0.52, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99, -9.99 ]

  else:
    print('NOT a valid reference for the solar composition (ags2005, husser, basti)')
    nans = np.empty((99,))
    nans[:] = np.nan
    return (symbol, mass, nans )
	      
  sol[0] = 1.0
  for i in range(len(sol)-1): sol[i+1] = 10.**(sol[i+1]-12.0)

  return (symbol,mass,sol)

def xyzmass(abu):

  """computing the mass fractions for hydrogen, helium and metals X,Y,Z
  """

  symbol, mass, sol = elements()

  total = sum(np.array(mass)*np.array(abu))

  x = mass[0]*abu[0]/total
  y = mass[1]*abu[1]/total
  z = np.sum(np.array(mass[2:])*np.array(abu[2:])/total)

  return(x,y,z)


def keepz(abu,zs,dex):

  """finds out the array of abundances that keeps the abundance changes in
     the array dex for the elements with the atomic number zs
     and at the same time keeps the metal mass fraction
     constant by modifying the metal abundances
  """

  assert (len(zs) == len(dex)),'zs and dex must have the same number of elements'

  abu1 = abu.copy()
  i = 0
  for entry in zs:
    abu1[entry-1] = abu[entry-1] + dex[i]
    i = i + 1

  x0, y0, z0 = xyzmass(abu)
  x1, y1, z1 = xyzmass(abu1)
  print(z0,z1)
  diff = np.array(abu1[zs]) - np.array(abu[zs])
  zr = np.abs(z1/z0)
  abu1[2:] = abu[2:]/zr 
  diff1 = np.array(abu1[zs]) - np.array(abu[zs])

  while np.max(diff1 - diff) > 0.0001:
    print(np.max(diff1 - diff))
    x1, y1, z1 = xyzmass(abu1)
    zr = np.abs(z1/z0)
    abu1[2:] = abu[2:]/zr
    diff1 = np.array(abu1) - np.array(abu)

  return()
      


  

def lgconv(xinput, yinput, fwhm, ppr=None):

  """convolution with a Gaussian in linear lambda scale
  for a constant resolution

  Parameters
  ----------
  xinput: numpy float array
      wavelengths 
  yinput: numpy array of floats
      fluxes
  fwhm: float
      FWHM of the Gaussian (same units as for xinput)
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default None, to keep the original sampling)

  Returns
  -------
  x: numpy float array
      wavelengths after convolution, will be a subset of xinput when that is linear, 
      otherwise a subset of the linearly resampled version
  y: numpy array of floats
      fluxes after convolution

  """

  #resampling to a linear lambda wavelength scale if need be
  xx = np.diff(xinput)
  if np.max(xx) - np.min(xx) > 1.e-7: #input not linearly sampled
    nel = len(xinput)
    minx = np.min(xinput)
    maxx = np.max(xinput)
    x = np.linspace(minx,maxx,nel)
    y = np.interp( x, xinput, yinput)
    #y = interp_spl( x, xinput, yinput)
  else:                       #input linearly sampled
    x = xinput
    y = yinput

  step = x[1] - x[0]

  assert (fwhm > 2*step), 'cannot convolve since fwhm is <= 2*step'

  sigma=fwhm/2.0/np.sqrt(-2.0*np.log(0.5))
  npoints = 2*int(3*fwhm/2./step)+1
  half = npoints * step /2.
  xx = np.linspace(-half,half,npoints)
  kernel = np.exp(-(xx-np.mean(xx))**2/2./sigma**2)
  kernel = kernel/np.sum(kernel)

  y = np.convolve(y,kernel,'valid')
  #y = ss.fftconvolve(y,kernel,'valid')
  print(npoints)
  edge = int(npoints/2)
  x = x[edge:-edge]

  print(xinput.size,x.size,y.size)

  if ppr != None:
    fac = int(fwhm / step / ppr)
    assert (fac != 0),'cannot resample since fac is = 0, increase fwhm or reduce ppr'
    subset = np.arange(x.size / fac, dtype=int) * fac 
    x = x[subset]
    y = y[subset]

  return(x,y)

def vgconv(xinput,yinput,fwhm, ppr=None):


  """convolution with a Gaussian in log lambda scale
  for a constant resolving power

  Parameters
  ----------
  xinput: numpy float array
      wavelengths 
  yinput: numpy array of floats
      fluxes
  fwhm: float
      FWHM of the Gaussian (km/s)
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default None, to keep the original sampling)

  Returns
  -------
  x: numpy float array
      wavelengths after convolution, will be a subset of xinput when that is equidistant
      in log lambda, otherwise a subset of the resampled version
  y: numpy array of floats
      fluxes after convolution

  """
  #resampling to ln(lambda) if need be
  xx = np.diff(np.log(xinput))
  if np.max(xx) - np.min(xx) > 1.e-7:  #input not equidist in loglambda
    nel = len(xinput)
    minx = np.log(xinput[0])
    maxx = np.log(xinput[-1])
    x = np.linspace(minx,maxx,nel)
    step = x[1] - x[0]
    x = np.exp(x)
    y = np.interp( x, xinput, yinput)
    #y = interp_spl( x, xinput, yinput)
  else:
    x = xinput
    y = yinput
    step = np.log(xinput[1])-np.log(xinput[0])

  fwhm = fwhm/clight # inverse of the resolving power

  assert (fwhm > 2*step), 'cannot convolve since fwhm is <= 2*step'

  sigma=fwhm/2.0/np.sqrt(-2.0*np.log(0.5))
  npoints = 2*int(3*fwhm/2./step)+1
  half = npoints * step /2.
  xx = np.linspace(-half,half,npoints)
  kernel = np.exp(-(xx-np.mean(xx))**2/2./sigma**2)
  kernel = kernel/np.sum(kernel)

  y = np.convolve(y,kernel,'valid')
  edge = int(npoints/2)
  x = x[edge:-edge]

  #print(xinput.size,x.size,y.size)

  if ppr != None:
    fac = int(fwhm / step / ppr)
    assert (fac != 0),'cannot resample since fac is = 0, increase fwhm or reduce ppr'
    print(fwhm,step,ppr,fac)
    subset = np.arange(x.size / fac, dtype=int) * fac 
    x = x[subset]
    y = y[subset]


  return(x,y)

def rotconv(xinput,yinput,vsini, ppr=None):


  """convolution with a Rotation profile 

  Parameters
  ----------
  xinput: numpy float array
      wavelengths 
  yinput: numpy array of floats
      fluxes
  vsini: float
      projected rotational velocity (km/s)
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default None, to keep the original sampling)

  Returns
  -------
  x: numpy float array
      wavelengths after convolution, will be a subset of xinput when that is equidistant
      in log lambda, otherwise a subset of the resampled version
  y: numpy array of floats
      fluxes after convolution

  """

  #resampling to ln(lambda) if need be
  xx = np.diff(np.log(xinput))
  if np.max(xx) - np.min(xx) > 1.e-7:  #input not equidist in loglambda
    nel = len(xinput)
    minx = np.min(np.log(xinput))
    maxx = np.max(np.log(xinput))
    x = np.linspace(minx,maxx,nel)
    step = x[1] - x[0]
    x = np.exp(x)
    y = np.interp( x, xinput, yinput)
    #y = interp_spl( x, xinput, yinput)
  else:
    x = xinput
    y = yinput

  deltamax=vsini/clight

  assert (deltamax > 2*step), 'cannot convolve since vsini is <= 2*step'

  npoints = 2*int(deltamax/step)+1
  xx = np.linspace(-deltamax,deltamax,npoints)
  c1=2.0*(1.0-epsilon)/np.pi/(1.0-epsilon/3.0)/deltamax
  c2=0.5*epsilon/(1.0-epsilon/3.0)/deltamax
  r2=(xx/deltamax)**2
  kernel = c1*np.sqrt(1.0-r2)+c2*(1.0-r2)
  kernel = kernel/np.sum(kernel)


  y = np.convolve(y,kernel,'valid')
  print(xinput.size,x.size,y.size)
  edge = int(npoints/2)
  x = x[edge:-edge]

  if ppr != None:
    fac = int(deltamax / step / ppr)
    assert (fac != 0),'cannot resample since fac is = 0, reduce fwhm or ppr'
    subset = np.arange(x.size / fac, dtype=int) * fac 
    x = x[subset]
    y = y[subset]

  return(x,y)

def smooth(x,n):    

  """Smooth using a Svitzky-Golay cubic filter


  Parameters
  ----------
  x: arr
    input array to smooth
  n: int
    window size
  """

  x2 = savgol_filter(x, n, 3)

  return(x2)


def gsynth(synthfile,fwhm=0.0,units='km/s',ebv=0.0,r_v=3.1,rv=0.0,
    outsynthfile=None,ppr=5,wrange=None,freeze=None):

  """Convolve and/or redden spectra in a FERRE grid

  Parameters
  ----------
  synthfile: str
      name of the input FERRE synth file 
  fwhm: float, can be an iterable
      FWHM of the Gaussian kernel (in A or km/s) for convolution
      (default 0.0, which means no convolution is performed)      
  units: str
      units for the FWHM ('A' for a constant resolution in Angstroms, 
      'km/s' for a constant resolution in velocity
      (default is 'km/s')
  ebv: float or iterable with floats
      E(B-V) to be applied to the model in the grid
      (default is 0.0)
  r_v: float or iterable with floats
      ratio of total to V-band extinction A_V/E(B-V)

  rv: float or iterable with floats
      Radial velocity to the applied to the model in the grid
      (default is 0.0)
  outsynthfile: str
      name of the output FERRE synth file
      (default is the same as synth file, but starting with 'n')
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default is 5, set to None to keep the original sampling)
  wrange: tuple
      Starting and ending wavelengths (if a smaller range that 
      the input's is desired)
      (default None, to keep the original range)
  freeze: dictionary
      Allows to reduce the dimensionality of the grid. The keys are the labels
      of the dimensions to freeze (as given in in the header of the input grid) 
      with the values that should be adopted for those 'frozen' dimensions. 
      Example: set freeze = {'TEFF': 5000.} to fix that value for the Teff dimension
      in a grid.
      (default None, to retain all the original dimensions)
  Returns
  -------
  writes outsynthfile with the smooth spectra

  """
  
  from synple import vgconv
  import numpy as np
  if np.max(ebv) > 0.0: from extinction import apply,ccm89 
  clight = 299792.458

  if outsynthfile is None: 
    assert synthfile[0] != 'n', 'default output file name starts with n_. Given that the input starts with n too, please choose an anternative outsynthfile'
    outsynthfile='n'+synthfile[1:]
  logw=0

  #read header, update and write out
  fin = open(synthfile,'r')
  fout = open(outsynthfile,'w')
  hd = []
  labels = []
  type = "'regular'"
  line = fin.readline()
  hd.append(line)
  while line[1] != "/":
    line = fin.readline()
    if "N_P" in line: n_p = np.array(line.split()[2:],dtype=int)
    if "STEPS" in line: steps = np.array(line.split()[2:],dtype=float)
    if "LLIMITS" in line: llimits = np.array(line.split()[2:],dtype=float)
    if "LABEL" in line: labels.append(line.split()[-1][1:-1])
    if "NPIX" in line: npix = int(line.split()[2])
    if "N_OF_DIM" in line: ndim = int(line.split()[2])
    if "WAVE" in line: wave = np.array(line.split()[2:],dtype=float)
    if "LOGW" in line: logw = int(line.split()[2]) 
    if "RESOLUTION" in line: resolution = float(line.split()[2])
    if "TYPE" in line: type = str(line.split()[2])
    if "NTOT" in line: ntot = int(line.split()[2])
    hd.append(line)

  if 'irregular' in type:
    assert (len(labels) == ndim), 'The number of LABELS in the header does not agree with the dimension of the grid'
  else:
    assert (len(n_p) == len(steps) & len(n_p) == len(llimits) & len(n_p) == len(labels) & len(n_p) == ndim), 'The dimension of the parameters from the header are inconsistent'

  assert (units == 'km/s' or units == 'A'), 'units must be either km/s or A'

  #update header parameters
  x = np.arange(npix)*wave[1]+wave[0]
  if logw == 1: x=10.**x
  if logw == 2: x=np.exp(x)
  
  newcol = []
  try: 
    nebv = len(ebv)
    ebvs = ebv
    print(ebvs)
    print(ndim,labels)
    #check they are uniformly spaced
    debv = np.diff(ebvs)
    print(np.max(debv),np.min(debv))
    assert np.max(debv) - np.min(debv) < 1.e-7, 'ebv values are not linearly spaced!'
    n_p = np.append(n_p,nebv)
    steps = np.append(steps,ebvs[1]-ebvs[0])
    llimits = np.append(llimits,ebvs[0])
    labels.append('E(B-V)')
    ndim = ndim + 1
    newcol.append(ndim-1)
    #update RESOLUTION?
    print(ndim,labels)
  except TypeError:
    nebv = 1
    ebvs = [ ebv ]

  try:
    nrv = len(rv)
    rvs = rv 
    print(rvs)
    print(ndim,labels)
    #check they are uniformly spaced
    drv = np.diff(rvs)
    print(np.max(drv),np.min(drv))
    assert np.max(drv) - np.min(drv) < 1.e-7, 'rv values are not linearly spaced!'
    n_p = np.append(n_p,nrv)
    steps = np.append(steps,rvs[1]-rvs[0])
    llimits = np.append(llimits,rvs[0])
    labels.append('RV')
    ndim = ndim + 1
    newcol.append(ndim-1)
    print(ndim,labels)
  except TypeError:
    nrv = 1
    ervs = [ rv ]

  try: 
    nfwhm = len(fwhm)
    fwhms = fwhm
    print(fwhms)
    print(ndim,labels)
    #check they are uniformly spaced
    dfwhm = np.diff(fwhms)
    print(np.max(dfwhm),np.min(dfwhm))
    assert np.max(dfwhm) - np.min(dfwhm) < 1.e-7, 'fwhm values are not linearly spaced!'
    n_p = np.append(n_p,nfwhm)
    steps = np.append(steps,fwhms[1]-fwhms[0])
    llimits = np.append(llimits,fwhms[0])
    labels.append('FWHM')
    ndim = ndim + 1
    newcol.append(ndim-1)
    minfwhm=np.min(fwhms)
    #update RESOLUTION?
    print(ndim,labels)
  except TypeError:
    nfwhm = 1
    fwhms = [ fwhm ]

  
  #define indices for grid loops
  if 'irregular' in type:
    ind = range(ntot)
    ind_n_p =  list(range(ndim))
    labels2 = list(labels)
  else:
    ll = []
    ind_n_p = []
    i = 0
    print('labels=',labels)
    labels2 = []
    for entry in labels:
      if freeze is not None:   
        lfkeys = list(freeze.keys())
        if entry not in lfkeys: 
            ind_n_p.append(i)
            labels2.append(entry)
      else:
        ind_n_p.append(i)
        labels2.append(entry)
      ll.append(np.arange(n_p[i]))
      i = i + 1
    ind = np.array(list(product(*ll)))
  
  if wrange is not None:
    assert (len(wrange) == 2), 'Error: wrange must have two elements'
    section1 = np.where( (x >= wrange[0]*(1.-10.*np.max(fwhm)/clight)) & (x <= wrange[1]*(1.+10.*np.max(fwhm)/clight)) )
    x = x[section1]
    npix = len(x)
    
  if np.min(fwhms) > 1.e-7:
    y = np.ones(npix)
    if units == 'km/s':
      print('min(fwhm)=',np.min(fwhms))
      xx,yy = vgconv(x,y,np.min(fwhms),ppr=ppr)
    else:
      xx,yy = lgconv(x,y,np.min(fwhms),ppr=ppr)
  else:
    print('Warning -- fwhm <= 1.e-7, no convolution will be performed, ppr will be ignored')
    xx = x
  
  
  if wrange is not None: 
    section2 = np.where( (xx >= wrange[0]) & (xx <= wrange[1]) ) 
    xx = xx [section2]
    
  #print(x,xx)
  #print(len(x),len(xx))
  #print(len(section1),len(section2))
  
  jlabel = 0
  for line in hd:
    if "N_OF_DIM" in line: line = " N_OF_DIM = "+str(len(ind_n_p))+"\n"    
    if "N_P" in line: line = " N_P = "+' '.join(map(str,n_p[ind_n_p]))+"\n"   
    if "STEPS" in line: line = " STEPS = "+' '.join(map(str,steps[ind_n_p]))+"\n"   
    if "LLIMITS" in line: line = " LLIMITS = "+' '.join(map(str,llimits[ind_n_p]))+"\n"
    if "LABEL" in line: 
      if jlabel == 0:
        for entry in labels2:
          jlabel = jlabel + 1
          ilabel = "'"+entry+"'"
          line = " LABEL("+str(jlabel)+") = "+ilabel+"\n"
          print('line=',line)
          fout.write(line)
        continue
      else:
        continue
    if "NPIX" in line: line = " NPIX = "+str(len(xx))+"\n"
    if "WAVE" in line: 
      if units == 'km/s': 
        line = " WAVE = "+str(np.log10(xx[0]))+" "+str(np.log10(xx[1])-np.log10(xx[0]))+"\n"
      else:
        line = " WAVE = "+str(xx[0])+" "+str(xx[1]-xx[0])+"\n"
    if "LOGW" in line: 
      if units == 'km/s':
        line = " LOGW = 1 \n"
      else:
        line = " LOGW = 0 \n"
    if "RESOLUTION" in line: 
        if units == 'km/s': 
            line = " RESOLUTION = "+str(clight/np.sqrt(clight**2/resolution**2 + np.min(fwhms)**2))+"\n"
        else:
            line = " RESOLUTION = "+str(np.mean(xx)/np.sqrt(np.mean(xx)**2/resolution**2 + np.min(fwhms)**2))+"\n"
    if line[1] != "/": fout.write(line)

  try: resolution
  except NameError: 
        if units == 'km/s': 
             line = " RESOLUTION = "+str(clight/np.min(fwhms))+"\n"
        else:
             line = " RESOLUTION = "+str(np.mean(wrange)/np.min(fwhms))+"\n"
        fout.write(line)

  fout.write(" /\n")

  #smooth and write data
  k = 0 #increases only when a line from the original file is used
  j = 0 #increases as we advance through the array ind
  if 'irregular' in type: 
    pass
  else:
    ntot = np.prod(n_p)
  for i in ind:
    j = j + 1
    print('line ',j,' of ',ntot)
    #print(k,ntot,i)
    if len(newcol) == 0 or all(i[newcol] == 0):
      line = fin.readline()
    if freeze is not None:
      skip = True
      for entry in lfkeys: 
        if (abs(freeze[entry] - par[labels.index(entry)]) < 1e-6): skip = False
      if skip: continue
    y = np.array(line.split(),dtype=float)
    if 'irregular' in type:
      par = y[:ndim]
      y = y[ndim:]
    else:
      #print(i,steps,llimits)
      par = i*steps+llimits
    print('par=',par)
    #print('newcol=',newcol)
    #print('i[newcol]=',i[newcol])
    
    #print('len(y)=',len(y))
    if wrange is not None: y = y [section1]
      
    #apply Gaussian convolution
    if 'FWHM' in labels:
      w = np.where(np.array(labels) == 'FWHM')
      fwhmval = par[w[0][0]]
    else:
      fwhmval = fwhm
    #print('fwhmval=',fwhmval)
    if fwhmval > 1.e-7:
      if units == 'km/s':
        xx,yy = vgconv(x,y,fwhmval,ppr=ppr)
      else:
        xx,yy = lgconv(x,y,fwhmval,ppr=ppr)
    else:
      xx,yy = x, y          
      
    #apply extinction
    if 'E(B-V)' in labels:
      w = np.where(np.array(labels) == 'E(B-V)')
      ebvval = par[w[0][0]]
      #print(ebvval)
      yy = apply(ccm89(xx, ebvval* 3.1, 3.1), yy)

    #apply RV  
    if 'RV' in labels:
      w = np.where(np.array(labels) == 'RV')
      rvval = par[w[0][0]]
      #print(rvval)
      yy = np.interp(xx, xx*(1.+rvval/clight), yy)

		
    if wrange is not None: yy = yy[section2]
    
    if 'irregular' in type: yy = np.insert(yy,0,par)
    yy.tofile(fout,sep=" ",format="%0.4e")
    fout.write("\n")
    k = k + 1

  fin.close()
  fout.close()
  

def gsynth_old(synthfile,fwhm=0.0,units='km/s',outsynthfile=None,ppr=5,wrange=None,freeze=None):

  """Smooth the spectra in a FERRE grid by Gaussian convolution

  Parameters
  ----------
  synthfile: str
      name of the input FERRE synth file 
  fwhm: float
      FWHM of the Gaussian kernel in A or km/s
      (default is 0.0, which means no convolution is performed)
  units: str
      units for the FWHM ('A' for a constant resolution in Angstroms, 
      'km/s' for a constant resolution in velocity
      (default is 'km/s')
  outsynthfile: str
      name of the output FERRE synth file
      (default is the same as synth file, but starting with 'n')
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default is 5, set to None to keep the original sampling)
  wrange: tuple
      Starting and ending wavelengths (if a smaller range that 
      the input's is desired)
      (default None, to keep the original range)
  freeze: dictionary
      Allows to reduce the dimensionality of the grid. The keys are the labels
      of the dimensions to freeze (as given in in the header of the input grid) 
      with the values that should be adopted for those 'frozen' dimensions. 
      Example: set freeze = {'TEFF': 5000.} to fix that value for the Teff dimension
      in a grid.
      (default None, to retain all the original dimensions)
  Returns
  -------
  writes outsynthfile with the smooth spectra

  """

  if outsynthfile is None: 
    assert synthfile[0] != 'n', 'default output file name starts with n_. Given that the input starts with n too, please choose an anternative outsynthfile'
    outsynthfile='n'+synthfile[1:]
  logw=0

  #read header, update and write out
  fin = open(synthfile,'r')
  fout = open(outsynthfile,'w')
  hd = []
  labels = []
  line = fin.readline()
  hd.append(line)
  while line[1] != "/":
    line = fin.readline()
    if "N_P" in line: n_p = np.array(line.split()[2:],dtype=int)
    if "STEPS" in line: steps = np.array(line.split()[2:],dtype=float)
    if "LLIMITS" in line: llimits = np.array(line.split()[2:],dtype=float)
    if "LABEL" in line: labels.append(line.split()[-1][1:-1])
    if "NPIX" in line: npix = int(line.split()[2])
    if "N_OF_DIM" in line: ndim = int(line.split()[2])
    if "WAVE" in line: wave = np.array(line.split()[2:],dtype=float)
    if "LOGW" in line: logw = int(line.split()[2]) 
    if "RESOLUTION" in line: resolution = float(line.split()[2])
    hd.append(line)

  assert (len(n_p) == len(steps) & len(n_p) == len(llimits) & len(n_p) == len(labels) & len(n_p) == ndim), 'The dimension of the parameters from the header are inconsistent'

  assert (units == 'km/s' or units == 'A'), 'units must be either km/s or A'

  #update header parameters
  x = np.arange(npix)*wave[1]+wave[0]
  if logw == 1: x=10.**x
  if logw == 2: x=np.exp(x)
  
  #define indices for grid loops
  ll = []
  ind_n_p = []
  i = 0
  for entry in labels:
    if freeze is not None:   
      lfkeys = list(freeze.keys())
      if entry not in lfkeys: ind_n_p.append(i)
    else:
      ind_n_p.append(i)
    ll.append(np.arange(n_p[i]))
    i = i + 1
  ind = list(product(*ll))
  
  if wrange is not None:
    assert (len(wrange) == 2), 'Error: wrange must have two elements'
    section1 = np.where( (x >= wrange[0]*(1.-10.*fwhm/clight)) & (x <= wrange[1]*(1.+10.*fwhm/clight)) )
    x = x[section1]
    npix = len(x)
    
  if fwhm > 1.e-7:
    y = np.ones(npix)
    if units == 'km/s':
      print('fwhm=',fwhm)
      xx,yy = vgconv(x,y,fwhm,ppr=ppr)
    else:
      xx,yy = lgconv(x,y,fwhm,ppr=ppr)
  else:
    print('Warning -- fwhm <= 1.e-7, no convolution will be performed, ppr will be ignored')
    xx = x
  
  
  if wrange is not None: 
    section2 = np.where( (xx >= wrange[0]) & (xx <= wrange[1]) ) 
    xx = xx [section2]
    
  #print(x,xx)
  #print(len(x),len(xx))
  
  jlabel = 0
  for line in hd:
    if "N_OF_DIM" in line: line = " N_OF_DIM = "+str(len(ind_n_p))+"\n"    
    if "N_P" in line: line = " N_P = "+' '.join(map(str,n_p[ind_n_p]))+"\n"   
    if "STEPS" in line: line = " STEPS = "+' '.join(map(str,steps[ind_n_p]))+"\n"   
    if "LLIMITS" in line: line = " LLIMITS = "+' '.join(map(str,llimits[ind_n_p]))+"\n"
    if freeze is not None:
      if "LABEL" in line:
        ilabel = line.split()[-1][1:-1] #drop starting/ending quotes
        if ilabel in lfkeys:
          continue
        else:
          jlabel = jlabel + 1
          line = " LABEL("+str(jlabel)+") = "+ilabel+"\n"
    if "NPIX" in line: line = " NPIX = "+str(len(xx))+"\n"
    if "WAVE" in line: 
      if units == 'km/s': 
        line = " WAVE = "+str(np.log10(xx[0]))+" "+str(np.log10(xx[1])-np.log10(xx[0]))+"\n"
      else:
        line = " WAVE = "+str(xx[0])+" "+str(xx[1]-xx[0])+"\n"
    if "LOGW" in line: 
      if units == 'km/s':
        line = " LOGW = 1 \n"
      else:
        line = " LOGW = 0 \n"
    if "RESOLUTION" in line: 
        if units == 'km/s': 
            line = " RESOLUTION = "+str(clight/np.sqrt(clight**2/resolution**2 + fwhm**2))+"\n"
        else:
            line = " RESOLUTION = "+str(np.mean(wrange)/np.sqrt(np.mean(wrange)**2/resolution**2 + fwhm**2))+"\n"
    if line[1] != "/": fout.write(line)

  try: resolution
  except NameError: 
        if units == 'km/s': 
             line = " RESOLUTION = "+str(clight/fwhm)+"\n"
        else:
             line = " RESOLUTION = "+str(np.mean(wrange)/fwhm)+"\n"
        fout.write(line)

  fout.write(" /\n")

  #smooth and write data
  k = 0
  j = 0
  ntot = np.prod(n_p)
  for i in ind:
    j = j + 1
    print('line ',j,' of ',ntot)
    #print(k,ntot,i)
    #print(i,steps,llimits)
    par = i*steps+llimits
    line = fin.readline()
    if freeze is not None:
      skip = True
      for entry in lfkeys: 
        if (abs(freeze[entry] - par[labels.index(entry)]) < 1e-6): skip = False
      if skip: continue
    y = np.array(line.split(),dtype=float)
    if wrange is not None: y = y [section1]
    if fwhm > 1.e-7:
      if units == 'km/s':
        xx,yy = vgconv(x,y,fwhm,ppr=ppr)
      else:
        xx,yy = lgconv(x,y,fwhm,ppr=ppr)
    else:
      xx,yy = x, y 
    if wrange is not None: yy = yy[section2]
    yy.tofile(fout,sep=" ",format="%0.4e")
    fout.write("\n")
    k = k + 1

  fin.close()
  fout.close()

def xplsf(synthfile,outsynthfile=None,ppr=5):
	
  """Takes an input FERRE grid with sufficient resolution (R>1000)
  and produces an output one with wavelength-dependent resolution
  approximating the Gaia DR3 XP data
  
  Parameters
  ----------
  synthfile: str
    Name of the Input FERRE synthfile 
    Must have R>1000 and a minimum wavelength coverage between 360-990 nm
  outsynthfile: str
    Name of the Output FERRE synthfile
    (default is None)
  ppr: float
    Points per resolution element for the output grid
    
  
  Returns
  -------
  Creates a synthfile with the nominal resolution of the Gaia XP (DR3) data
  
  """

  #Table describing the resolution R of the externally calibrated spectra (ECS) XP DR3 data in Montegriffo et al. 2023
  x = np.array([350.,370.,390.,410,430.,450.,470.,490.,510.,530.,550.,570.,590.,610.,630.,640.,660.,680.,700.,720.,740.,760.,780.,800.,820.,840.,860.,880.,900.,920.,940.,960.,980.])
  y = np.array([71.5,67.2,61.4,55.3,49.9,45.0,41.4,38.0,35.0,32.5,30.0,28.3,25.5,25.2,22.0,74.6,77.8,77.8,76.8,74.5,73.2,71.1,68.8,66.9,65.1,63.4,62.1,60.7,58.2,57.2,54.9,53.8,52.0])
  y = x/y # FWHM = lambda/R
  plt.plot(x,x/y,linewidth=3)
  plt.xlabel('wavelength (nm)')
  plt.ylabel('R')
  plt.show()
  
  x = x*10.
  y = y*10.

  vgsynth(synthfile,x,y,outsynthfile=outsynthfile,wrange=(3600.,9900),ppr=ppr)

  return()

def vgsynth(synthfile,wavelength,fwhm,outsynthfile=None,ppr=5,wrange=None,original=False):

  """Variable-width Gaussian convolution
  This is similar to gsynth but the FWHM of the Gaussian kernel can change
  with wavelength.
  
  Parameters
  ----------
  synthfile: str
      name of the input FERRE synth file
  wavelength: array of floats
      Wavelength (angstroms)
  fwhm: array of floats
      FWHM of the Gaussian kernel (in A) for convolution
  outsynthfile: str
      name of the output FERRE synth file
      (default is None, to set the same as synthfile, but starting with 'v')
  ppr: float, optional
      Points per resolution element to sample the convolved spectrum
      (default is 5, set to None to keep the original sampling)
  wrange: tuple
      Starting and ending wavelengths (if a smaller range that
      the input's is desired)
      (default None, to keep the original range of the wavelength/fwhm arrays)
  original: bool
      Switch to retain the wavelength sample of the original grid in the output one
      (default False, and the FWHM is resampled with ppr points)

  Returns
  -------
  writes outsynthfile with the smooth spectra
  
  """

  h, p, d = read_synth(synthfile)
  xx = lambda_synth(synthfile)
  
  ending = synthfile.find('.dat')
  if ending < -1: 
      ending = synthfile.find('.pickle')
  if ending < -1:
      ending = len(synthfile) + 1
  root = synthfile[2:ending]
  if outsynthfile is None: outsynthfile = 'v_'+root+'.dat'
    
  gg = np.zeros((len(xx),len(xx)))

  for i in range(len(xx)):
     sigma=np.interp(xx[i],wavelength,fwhm)/2.0/np.sqrt(-2.0*np.log(0.5))
     kernel = np.exp(-(xx-xx[i])**2/2./sigma**2)
     kernel = kernel/np.sum(kernel)
     #print(xx[i],sigma)
     gg[i,:] = kernel/np.sum(kernel)

  rr = np.matmul(gg,np.transpose(d))
  d2 = np.transpose(rr)

  if original:
     xff = xx.copy()
     d3 = d2.copy() 
  else:
    if wrange is None: wrange = (wavelength[0],wavelength[-1])
    xf = wrange[0]
    xff = [xf]
    while xf <= wrange[1]:
      xf = xf + np.interp(xf,wavelength,fwhm)/ppr
      xff.append(xf)

    xff = np.array(xff)    
    d3 = np.zeros( (len(d[:,0]),len(xff)) )
    for i in range(len(d[:,0])):
      d3[i,:] = np.interp(xff,xx,d2[i,:])
      

  h3 = h.copy()
  h3['NPIX'] = str(len(xff))
  h3['RESOLUTION'] = str(np.mean(wavelength/fwhm))
  h3['WAVELENGTHS'] = ' '.join(map(str,xff))
  #h3['LOGW'] = str(1)
  h3.pop('LOGW')
  #h3['WAVE'] = str(np.log10(np.min(xff)))+' '+str( ( np.log10(np.max(xff)) - np.log10(np.min(xff)) ) / len(xff) )
  h3.pop('WAVE')
  h3['COMMENTS5'] = "'Data smoothed by Variable-FWHM Gausssian convolution (synple.vgsynth)'"
  #h3['COMMENTS6'] = "'wavelenghts are NOT linear in loglambda'"
  #h3['COMMENTS7'] = "'wavelenghts are as following:"+' '.join(map(str,xff))+"'"


  write_synth(outsynthfile,p,d3,hdr=h3)
  xff.tofile(outsynthfile+'.lambda',sep=" ",format="%s")

  return()
  
  
  
def fit(xdata, ydata, modelfile, params, bounds,  
        vmicro=1.0, abu=None, vrot=0.0, fwhm=0.0, vmacro=0.0,
        dw=None, strength=1e-4, linelist=linelist0, atom='ap18',
        steprot=0.0, stepfwhm=0.0, lte=None, method='Powell', tol = 0.001):
  
  """
  Fitting a piece of the spectrum to find the optimal value of one or
  multiple params
  """      
  from scipy.optimize import minimize

  #must check that all params are included in possible_parnames

  atmostype, teff, logg, vmicro2, abu, nd, atmos = read_model(modelfile)   
  wrange = ( np.min(xdata) - 0.5 , np.max(xdata) + 0.5 )

  othernames = np.array(['modelfile','wrange','dw','strength', 
                         'linelist', 'atom', 'steprot', 'stepfwhm'])
  others = {}
  for entry in othernames:
	  others[entry] = locals()[entry]
 
  possible_parnames = np.array(['vmicro', 'vrot', 'fwhm', 'vmacro'])
  parnames = []
  parvalues = []
  for entry in possible_parnames:
      if entry in params:
         parnames.append(entry)
         parvalues.append(locals()[entry])
      else:
          others[entry] = locals()[entry]
       
  if len(params) > len(parnames):
      symbol, mass, sol = elements()
      for entry in symbol:
          if entry in params:
             elindex = int(np.where(np.array(symbol) == entry)[0])
             parnames.append(entry)
             parvalues.append(np.log10(abu[elindex])-np.log10(sol[elindex])+3.)
             print(entry,elindex)
     

  res = minimize(fun, np.array(parvalues), 
                 ( xdata, ydata, parnames, others), method=method ,
                 bounds = bounds, tol = tol )                 
  
  print('res=',res)
  return(res)
  		  

def fun(parvalues, *args ):
  """
  Auxiliary to the routine fit
  """

  (xdata, ydata, parnames, others) = args
  
  argnames = np.array(['modelfile','wrange','dw','strength', 
                         'linelist', 'atom', 'steprot', 'stepfwhm'])
  for entry in argnames: 
	  globals()[entry] = others[entry]
	  print('assigning ',entry,' to ',others[entry])
  
  possible_parnames = ('vmicro', 'vrot', 'fwhm', 'vmacro')
  k = 0
  for entry in possible_parnames:
      if entry in parnames: 
          globals()[entry] = parvalues[k]
          k = k + 1
      else:
          globals()[entry] = others[entry]
          
  atmostype, teff, logg, vmicro2, abu, nd, atmos = read_model(modelfile)

  if len(parvalues) > k:
      symbol, mass, sol = elements()
      for entry in symbol:
          if entry in parnames:
              elindex = int(np.where(np.array(symbol) == entry)[0])
              abu[elindex] = 10.**(parvalues[k] + np.log10(sol[elindex])-3.)
              k = k + 1
              #print('symbol[elindex]=',np.array(symbol)[elindex])
              #print('entry=',entry)
              #print('parnames=',parnames)


  print('parnames=',parnames)
  print('parvalues=',parvalues)		  
	    
  out = syn( modelfile, wrange , dw=dw, strength=strength , \
    vmicro=vmicro, abu=abu, \
    linelist=linelist, atom=atom, vrot=vrot, fwhm=fwhm, vmacro=vmacro, \
    steprot=steprot, stepfwhm=stepfwhm,  intensity=False, \
    lineid=False, tag=False,  \
    clean=True, save=False, synfile=None, lte=None, compute=True, tmpdir=None)

  chi = np.sum( (ydata - np.interp(xdata, out[0], out[1]/out[2]) )**2)
  
  #plt.clf()
  #plt.plot(xdata,ydata,xdata,np.interp(xdata, out[0], out[1]/out[2]) )
  #plt.show()
  print('chi=',chi)
  
  return(chi)


def cebas(p,d,flx,iva):

    """Contrast-Expansion for BAS
    Evaluates the chi-squared between an observed spectrum
    (flx, with inverse variance iva) and an entire model grid  (p,d)
    to determine the best-fitting parameters and model
    
    Parameters
    ----------
    p: 2D numpy array of floats 
      parameter table with as many columns as parameters and
      as many rows as spectra in the array d
      
    d: 2D numpy array of floats
      spectra table with as many columns as frequencies in the spectra
      and as many rows as spectra
      
    flx: 1D numpy array of floats
      observed spectrum 
    
    iva: 1D numpy array of float
      inverse variance for the observed spectrum
    
    Returns
    -------
    res: numpy array of floats
      best-fitting parameters (as many entries as columns in p)
      
    eres: numpy array of floats
      uncertainties for the best-fitting parameters
    
    cov: numpy array of floats
      top half of the covariance matrix for the best-fitting parameters
      
    bflx: numpy array of floats
      best-fitting model (same size as flx)

    weights: numpy array of floats
      likelihood for the model grid, with as many elements as rows in p and d 
      
    """

    chi = np.sum((d-flx)**2 * iva,-1)
    beta = np.median(chi) / 1490. / 5.
    #print('min/max/median chi=',np.min(chi),np.max(chi),np.median(chi))
    #print('beta=',beta)
    while np.exp(-np.min(chi)/2./beta) <= 0.0:
      beta = beta * 2.
      #print('-- new beta=',beta)

    #parameters
    ndim = len(p[0,:])
    res = np.zeros(ndim)
    eres = np.zeros(ndim)
    cov = np.zeros(ndim*(ndim+1)//2)
    likeli = np.exp(-chi/2./beta)
    den = np.sum(likeli)
    #print('den=',den)
    k = 0
    for i in range(ndim):
        #parameters
        res[i] = np.sum( likeli * p[:,i])/den
        
        #uncertainties
        for j in range(ndim-i):
          cov[k] = np.sum( likeli * (p[:,i] - res[i]) * (p[:,j+i] - res[j+i]) )/den
          if j == 0: eres[i] = np.sqrt(cov[k])
          k = k + 1
      
    #best-fitting model
    bflx = np.matmul(likeli,d)/den
    #bflx = [0.0,0.0]
        
    print('res=',res,'eres=',eres)
      
    return(res,eres,cov,bflx,likeli)


def bas(infile, synthfile=None, outfile=None, target=None, rv=None, ebv=None, 
        focus=False, star=True, conti=0, absolut=False, wrange=None):

    """Bayesian Algorithm in Synple
    
    Parameters
    ----------
    infile: str or list
      input FITS file, 
      root for the filenames of the input frd/err files,
      (and, if outfile is None, the output opf/mdl files,)
      or a string with wildcards (*?[]) that expand into multiple files,
      or a list of files
    synthfile: str
      name of the model grid
      (default is None and the code attempts to choose the appropriate
      grid according to the source of the input data)      
    outfile: str
      output FITS file or
      root for the filenames of the output opdf/mdf files
    target: iterable
      input list of numerals or targetids to select objects
      to process. If the list includes numbers < 10000, they are interpreted
      as the order of the targets in the input file(s). Otherwise they are 
      interpreted as a list of target ids. Either way only the target list
      will be analyzed and the others skipped.
      (default is None meaning that all the targets in the input file(s) are analyzed)
    rv: iterable
      this can be an iterable matching the length of target with the RVs to be corrected
      prior to the analysis. When equal to None, if RV is among the parameters in the synthfile, 
      RVs will be determined as such, but otherwise RVs are derived by the routine xxc.
      (default is None)
    ebv: iterable
      this can be an iterable matching the length of target with the E(B-V) values 
      to be corrected prior to the analysis. If equal to None, if E(B-V) is 
      among the parameters in the synthfile, it will be determined as such,  but 
      otherwise, only for DESI, it will be read from infile and corrected for
    focus: bool
      switch to activate a two-step algorithm in which a coarsely 
      subsampled version of the grid is used to identify first where  
      the optimal solution is, and then perform an focused analysis
      in that region (a +/- 3 sigma volume)
    star: bool
      switch to limit the analysis of DESI spectra to stars. It has no
      effect on other data sets. Activating target disables star.
      (default True)
    conti: int
      conti > 0 activates the continuum normalization (see 'continuum' function)
      by a running mean with a width of  'conti'
      NOTE that the default (0) is dividing the input/model fluxes in each 
      spectrum by their mean value
      (default 0)
    absolut: bool
      activates the output of the absolute fluxes for the best-fitting 
      model (.flx file) and the input (unnormalized) fluxes (.frd file)
    wrange: 2-element iterable
      spectral range to use in the fittings
      (default None, and sets wrange to the values of the adopted grid)
   
    Returns
    -------
    Creates output FERRE-formatted files with the normalized data (.nrd),
    best-fitting parameters (.opf) and best-fitting models (.mdl).
    
    """
    
    if type(infile) is list:
        infiles = infile
    else:
        if type(infile) is str:
            if '*' in infile or '?' in infile or '[' in infile:
              infiles = glob.glob(infile)
            else:
              infiles = [infile]

    instr, default_synthfile = identify_instrument(infiles[0])
    if synthfile is None: synthfile = default_synthfile
    print('data appear to be from ',instr)
    print('adopting synthfile ',synthfile)
    instr0 = instr

    #models    
    hd, p, d = read_synth(synthfile)      
    x = lambda_synth(synthfile)
    lenx = len(x)
    if type(hd) is list:
      hd0 = hd[1]
    else:
      hd0 = hd

    if wrange is not None:
      assert(len(wrange) == 2),'wrange must be a 2-element array' 
      w = (x >= wrange[0]) & (x <= wrange[1])
      d = d[:,w]
      x = x[w]
      lenx = len(x)
      
    #normalization
    print('normalizing grid...')
    if absolut: da = d.copy() # da keeps a copy of the original grid
    for entry in range(len(d[:,0])):
        if conti > 0:
          cc = continuum(d[entry,:],window_length=conti)
        else:
          cc = np.mean(d[entry,:])
        d[entry,:] = d[entry,:] / cc


    if focus:
      p2 = p
      d2 = d
      nmod = len(p[:,0])
      irnd = np.array(np.random.random_sample(int(nmod*0.1))*nmod,dtype=int)
      p = p[irnd,:]
      d = d[irnd,:]


    
    #sanity check  
    if len(infiles) > 1 and (rv is not None or ebv is not None):
      print('BAS warning: are the same rv/ebv arrays/lists intended for multiple infiles ...??')  

        
    for file in infiles:

      #data
      print('reading data from file '+file+'...')
      instr, default_synthfile = identify_instrument(file)
      if instr0 is None:
          print(rv,instr0)
          #FERRE files, expand lists of wavelengths arrays into a single array
          if type(x) is list: x = np.hstack(x)
      else:
          assert(instr == instr0),'all the input files must be from the same instrument'
      
      if ebv is not None:
        if 'EBV' in hd0.values() or 'E(B-V)' in hd0.values():
          ebv = 0


              
      ids, x2, obs, ivr, xtr = read_spec(file,wavelengths=x,target=target,rv=rv,
                                    ebv=ebv, star=star)
      lenx2 = len(x2)
      if ivr.ndim == 1: 
        obs = obs.reshape((1,lenx2))
        ivr = ivr.reshape((1,lenx2))
      nspec = len(obs[:,0])
      print('nspec in bas:',nspec)

      if outfile is None or len(infiles) > 1:
        outfile = file
      else:
        if type(outfile) is list:
          assert len(outfile) == 1,'outfile can only be specified when there is a single infile'

      opffile = outfile + '.opf'
      mdlfile = outfile + '.mdl'      
      nrdfile = outfile + '.nrd'
      errfile = outfile + '.err'
      wavfile = outfile + '.wav'
      fmpfile = outfile + '.fmp.fits'
      scrfile = outfile + '.scr.fits'
      if absolut:
        frdfile = outfile + '.frd'
        flxfile = outfile + '.flx'

      #open output parameter, observed and model file
      opf = open(opffile,'w')
      mdl = open(mdlfile,'w')
      nrd = open(nrdfile,'w')
      err = open(errfile,'w')
      wav = open(wavfile,'w')
      if absolut:
        frd = open(frdfile,'w')
        flx = open(flxfile,'w')
          

      for j in range(nspec):

        print('spectrum ',j,' of ',nspec,' in ',file)
        
        #clean the data
        spec = obs[j,:]
        www = np.where(np.isnan(spec))[0]
        #print('www:',www)
        if len(www) > 0:
          www2 = np.where(not np.isnan(spec))[0]
          xax = np.arange(lenx2)
          flx = np.interp(xax,xax[www2],spec[www2])

        #normalize
        if conti > 0:
          mspec = continuum(spec, window_length=conti)
          www = (mspec == 0.0)
          mspec[www] = 1. 
        else:
          mspec = np.mean(spec)
          if mspec == 0.0: mspec = np.median(spec)
          if mspec == 0.0: mspec = 1.

        if absolut: rawspec = spec.copy()
        spec = spec / mspec
        ivar = ivr[j,:] * mspec**2

        #analyze
        res, eres, cov, bmod, weights = cebas( p, d, spec, ivar )
        lchi = np.log10( np.sum((bmod-spec)**2 * ivar) / (len(bmod) - len(res)) )
        print('reduced lchi =',lchi)
                
        vrad = 0.0
        if rv is None and 'RV' not in hd0.values() and instr0 is not None:

          print('type(x2) is ',type(x2))
          vrad, evrad = xxc(x2,spec,ivar,x2,bmod)
          print('RV = ',vrad,' km/s')
        
          #correct RV and reanalyze
          spec = np.interp(x2, x2 * (1. - vrad/clight), spec)
          ivar = np.interp(x2, x2 * (1. - vrad/clight), ivar)
          res, eres, cov, bmod, weights = cebas( p, d, spec, ivar )
          
          lchi = np.log10( np.sum((bmod-spec)**2 * ivar) / (len(bmod) - len(res)) )
          print('reduced lchi =',lchi)

        if focus:
          eres[eres < 1e-17] = 1e-17 # avoid division by zero
          w = ( (abs(p2-res)/eres).max(1) < 3. )
          if len(np.where(w)[0]) > 0:
            res, eres, cov, bmod, weights = cebas( p2[w,:], d2[w,:], spec, ivar )
            lchi = np.log10( np.sum((bmod-spec)**2 * ivar) / (len(bmod) - len(res)) )
          print('focus selected ',len(np.where(w)[0]), 'points, giving a reduced lchi =',lchi)

        if absolut:
          den = np.sum(weights)
          abbmod = np.matmul(weights,da)/den

      
        opf.write(str(ids[j])+' '+' '.join(map(str,res))+' '+
            ' '.join(map(str,eres))+' '+
            str(vrad)+' '+str(np.median(spec*np.sqrt(ivar)))+' '+
            str(lchi)+' '+' '.join(map(str,cov))+'\n')
        nrd.write(' '.join(map(str,spec))+'\n')
        mdl.write(' '.join(map(str,bmod))+'\n')
        err.write(' '.join(map(str,1./np.sqrt(ivar)))+'\n')
        if absolut:
            frd.write(' '.join(map(str,rawspec))+'\n')
            flx.write(' '.join(map(str,abbmod))+'\n')
        if j == 0: wav.write(' '.join(map(str,x2))+'\n')
      

      print('closing opf file:',opffile)
      opf.close()
      mdl.close()
      nrd.close()
      err.close()
      wav.close()
      if absolut:
        frd.close()
        flx.close()
      
      if instr == 'DESI':
        head, fibermap, scores = xtr
        fmp = tbl.Table(fibermap)
        hdu0 = fits.BinTableHDU(fmp)
        hdu0.writeto(fmpfile)
        scr = tbl.Table(scores)
        hdu0 = fits.BinTableHDU(scr)
        hdu0.writeto(scrfile)

      
    return()
   
    
def identify_instrument(infile):
	
    """Identify the instrument that produced infile
       LAMOST, DESI, NGSL(STIS)/CALSPEC, Gaia XP, MILES, INT/IDS-R900V, 
       GTC/OSIRIS-R2500U
       
    Parameters
    ----------
    infile: str
      name of an input (FITS) data file
      
    Returns
    -------
    instr: str
      instrument/telescope from which the data comes
      
    grid: str
      default BAS grid to adopt for it (see config/bas-grids.yaml
      
    """
    
    conf = load_conf(config='bas-grids.yaml',confdir=confdir)

    instr = None
    path, filename = os.path.split(infile)
    if infile[-4:] == 'fits':
        fi = fits.open(infile)
        head = fi[0].header
        if 'TELESCOP' in head: 
            telescop = head['TELESCOP']
            if head['TELESCOP'] == 'LAMOST' and filename[:4] == 'spec' and len(fi) == 2:
               instr = 'LAMOST'
            if 'INSTRUME' in head:
                if head['TELESCOP'] == 'HST' and head['INSTRUME'][:4] == 'STIS':
                    instr = 'STIS'  
            if head['TELESCOP'][:3] == 'INT' and head['INSTRUME'][:3] == 'IDS' and \
               head['CAMERA'][:3] == '235' and head['GRATNAME'][:5] == 'R900V':
               instr = 'IDS-R900V'
            if head['TELESCOP'][:3] == 'GTC' and head['INSTRUME'][:6] == 'OSIRIS' and \
             head['GRISM'][:6] == 'R2500U':
               instr = 'OSIRIS-R2500U'
        else:
            if 'MAPKEY' in head:
                if head['MAPKEY'] == 'calspec':
                    instr = 'CALSPEC'
            if (filename[:5] == 'coadd' or filename[:7] == 'spectra') and len(fi) > 10:
                instr = 'DESI'
            if 'TIMEXPOS' in head:
                if head['TIMEXPOS'] == -999.0:
                  instr = 'MILES'
    
    if instr is None:
        grid = None
    else:            
        grid = conf[instr]
                
    return(instr,grid)

def read_knoao2005(wrange):
 
  """Solar Kurucz 2005 atlas"""

  file = os.path.join(atlasdir,'solarfluxintwl.fits')

  d = fits.open(file)[0].data
  w = d[:,0]*10.
  f = d[:,1]
  ind = np.where((w > wrange[0]) & (w < wrange[1]))

  return(w[ind],f[ind])


def read_iag(wrange):

  """Solar IAG atlas (Reiners et al. 2016 A&A 587, 65)"""
  
  if wrange[0] < 10300.:
    file = os.path.join(atlasdir,'iag-vis.fits')
  else:
    file = os.path.join(atlasdir,'iag-nir.fits')
  d = fits.open(file)[0].data    
  w = vac2air(d[:,0])
  f = d[:,1]
  ind = np.where((w > wrange[0]) & (w < wrange[1]))

  return(w[ind],f[ind])
  

def read_arcturus(wrange):

  """Arcturus atlas (Hinkle et al. 2000)"""

  file = os.path.join(atlasdir,'ardata.fits')

  d = fits.open(file)[1].data
  w = d['wavelength']
  f = d['arcturus']
  ind = np.where((w > wrange[0]) & (w < wrange[1]))

  return(w[ind],f[ind])


def read_spec(infile,wavelengths=None,target=None,rv=None,ebv=None,star=True):
    """Read and (if wavelengths is given) resample spectral observations
    
    Parameters
    ----------
    infile: str
      name of an input FITS file, or root for input FERRE-formatted files
      (frd,err)
    wavelenghts: numpy array of floats
      array with the wavelengths of a grid to resample the observations
    target: iterable of integers/longs
      list of targets to read -- can be either integers indicating
      the order of the targets of interest in the input infile or
      targetids (e.g. for DESI)
      (default is none)
    rv: iterable of floats
      this can be an iterable matching the length of target with the RVs to  
      be corrected. 
      When equal to None, velocities offsets are not considered.
      (default is None)
    ebv: iterable of floats
      this can be an iterable matching the length of target with the reddening 
      to be corrected. 
      When equal to None, no reddening correction is applied, except for 
      DESI data, for which the SFD values from the DESI files will be used.
      When equal to 0 (integer!), no reddening correction is applied for 
      DESI data, or data from any other source
    star: bool
      flag to pre-select only stars for DESI for having any of the 
      following targetting bits set:
      STD_FAINT, STD_WD, STD_BRIGHT, MWS_ANY or SCND_ANY. 
      It has no effect on other data sets apart from DESI.
      Cannot be combined with target; passing a list in target  
      sets star to False 
      (default is True)

      
    Returns
    -------
    ids: numpy array of str
      strings identifying the target(s). In some cases it may simply 
      be a number 
    wav: numpy array of floats
      common wavelength array for the spectra
      which will be identical to the input 'wavelengths' array if provided
      
    frd: numpy array of floats
      observed spectra (as many columns as frequencies and as many rows
      as spectra)
      
    ivr: numpy array of floats
      inverse variance for frd (same size as frd)
      
    xtr: tuple of objects
      the first element is usually a header dictionary
      for DESI it also contains fibermap and scores structures
      
    """

    from extinction import apply,remove,ccm89

    if target is not None:
      try:
        _ = (e for e in target)
      except TypeError:
        print('target must be None or an iterable')

    if rv is not None:
      try:
        _ = (e for e in rv)
        if target is not None:
          assert(len(rv) == len(target)),'rv and target must have the same length when both are input'
      except TypeError:
        pass
  
    if ebv is not None:
      try:
        _ = (e for e in ebv)
        if target is not None:
          assert(len(ebv) == len(target)),'ebv and target must have the same length when both are input'
      except TypeError:
        pass


    #data
    if infile[-4:] == 'fits':
      instr, synthfile = identify_instrument(infile)    
      if instr == 'LAMOST':
        #reading LAMOST spec file
        fi = fits.open(infile)
        head = fi[0].header
        s = fi[1].data        
        wav = np.transpose(s['WAVELENGTH'])[:,0]
        wav = vac2air(wav)
        lenwav = len(wav)
        flux = np.transpose(s['FLUX'])[:,0]
        ivar = np.transpose(s['IVAR'])[:,0]
        if 'OBJNAME' in head:
          ids = np.array([head['OBJNAME']])
        else:
          ids = np.array([infile])

        assert(target is None),'target must be None for single-target LAMOST files'
 
        xtr = (head) 
        wav, frd, ivr = single_target_prep(wav, flux, ivar, rv, ebv, wavelengths=wavelengths)

      elif instr == 'DESI':
         if wavelengths is not None and type(wavelengths) is not list:
           twavelengths = []
           twavelengths.append(wavelengths[(wavelengths < 5800.)])
           twavelengths.append(wavelengths[(wavelengths >= 5800.) & (wavelengths < 7600.)])
           twavelengths.append(wavelengths[wavelengths >= 7600.])
           wavelengths = twavelengths.copy()

         i = 0
         for band in ('B','R','Z'):
           wav1,flux1,ivar1,res1,head1,map1,scores1 = read_desispec(infile,band)
           wav1 = vac2air(wav1)

           if band == 'B': #check if there is a 'target' or 'star' preselection
             ind = []
             if target is None:
               ind = np.arange(len(map1))
               if len(ind) == 0:
                 print('no DESI targets in file ',infile)
                 return(None,None,None,None)
             else:
               if star:
                 print('passing a list in target disables star!')
                 star = False
               if np.max(target) < 10000:
                 #target is a list with the order of the desired spectra
                 ind = target
               else:
                 #target is a list with targetids
                 ind = np.where(np.isin(map1['targetid'],target))[0]

               if len(ind) == 0:
                 print('no DESI target in file ',infile,' matches the input target=',target)
                 return(None,None,None,None)

             if star:
               ind = []
               j = 0
               for entry in map1['desi_target']:
                 bits = desimask(entry)
                 print('bits=',bits)
                 if 'STD_FAINT' in bits or 'STD_WD' in bits or 'STD_BRIGHT' in bits or 'MWS_ANY' in bits or 'SCND_ANY' in bits:
                   ind.append(j)
                 j += 1

               if len(ind) == 0:
                 print('no DESI target in file ',infile,' has star targetting bits (STD_BRIGHT/FAINT, STD_WD, MWS_ANY or SCND_ANY)')
                 print('processing ALL targets')
                 ind = np.arange(len(map1)) 

             if len(ind) > 0:
               nspec = len(ind)
               print('read_spec: selecting targets with indices -- ',ind)
             else:
               nspec = len(map1)
                                   
             if rv is None:
               vrad = np.zeros(nspec)
             else:
               vrad = rv

             if type(ebv) is int and ebv == 0:
               red = np.zeros(nspec)
               print('E(B-V)=0 adopted for all the objects')
             else:
               if ebv is None:
                 red = map1['EBV']
                 if len(ind) > 0:
                   red = red[ind]
                 print('Correcting E(B-V) from SFD map')
               else:
                 red = ebv
                 print('Correcting E(B-V) from values provided by the user')
                            
           #limit the sample to target/star          
           if len(ind) > 0:
             #print('ind=',ind)
             flux1 = flux1[ind,:]
             ivar1 = ivar1[ind,:]
             res1 = res1[ind,:,:]
             map1 = map1[ind] 
             scores1 = scores1[ind]

           #correct reddening
           print('Reddening from SFD map has a mean E(B-V)=',
                 np.mean(map1['EBV']),' +/- ', np.std(map1['EBV']))

           for j in range(nspec):
             if np.abs(red[j]) > 1e-7:
               xtmp = np.array(wav1,dtype=float)
               ytmp = np.array(flux1[j,:],dtype=float)
               tmp = remove(ccm89(xtmp, red[j] * 3.1, 3.1), ytmp)
               ivar1[j,:] = ivar1[j,:] * (np.divide(ytmp,tmp,where=tmp>0))**2
               flux1[j,:] = tmp


           if wavelengths is None:
             nfreq = len(wav1)
             #print('nfreq=',nfreq)
             flux2 = np.zeros((nspec,nfreq))
             ivar2 = np.zeros((nspec,nfreq))
             for j in range(nspec):
               flux2[j,:] = np.interp(wav1,wav1*(1.+vrad[j]/clight),flux1[j,:])
               ivar2[j,:] = np.interp(wav1,wav1*(1.+vrad[j]/clight),ivar1[j,:]) 
             flux1 = flux2
             ivar1 = ivar2
           else:
             assert (type(wavelengths) is list),'A list is expected for the input wavelengths'
             nfreq = len(wavelengths[i])
             #print('nfreq=',nfreq)
             flux2 = np.zeros((nspec,nfreq))
             ivar2 = np.zeros((nspec,nfreq))
             for j in range(nspec):
               flux2[j,:] = np.interp(wavelengths[i],wav1*(1.+vrad[j]/clight),flux1[j,:])
               ivar2[j,:] = np.interp(wavelengths[i],wav1*(1.+vrad[j]/clight),ivar1[j,:]) 
             flux1 = flux2
             ivar1 = ivar2
             wav1 = wavelengths[i]
               
           if band == 'B':
             wav = wav1
             frd = flux1
             wbad = (wav1 >= 4300.) & (wav1 <= 4450.)
             ivar1[:,wbad] = 0.
             ivr = ivar1
             ids = map1['targetid']             
             xtr = (head1, map1, scores1)
           else:
             wav = np.concatenate((wav,wav1))
             frd = np.concatenate((frd,flux1),axis=1)
             ivr = np.concatenate((ivr,ivar1),axis=1)
           print(len(wav1),len(wav))

           i += 1

      elif instr == "CALSPEC" or instr == "STIS":
        fi = fits.open(infile)
        head = fi[0].header
        s = fi[1].data
        wav = s['WAVELENGTH']
        wav = vac2air(wav)
        lenwav = len(wav)
        flux = s['FLUX']
        if 'STATERROR' in s.names:
          err = s['STATERROR']
        elif 'STATERR' in s.names:
          err = s['STATERR']
        else:
          print('Warning: cannot find the statistical error in the file from HST')
          print('         assuming S/N = 20!')
          err = flux * 0.05
        if 'SYSERROR' in s.names: err = flux*0.000001 + err + s['SYSERROR']
        ivar4 = np.divide(1.,err**2, where = (err**2 > 0.) ,dtype = np.float128)
        ivar = np.float64(ivar4)
        if 'TARGETID' in head:
          ids = np.array([head['TARGETID']])
        elif 'TARGNAME' in head:
          ids = np.array([head['TARGNAME']])
        else:
          ids = np.array([infile])
        
        assert(target is None),'target must be None for STIS data (1 target per file)'
          
        xtr = (head)
        wav, frd, ivr = single_target_prep(wav, flux, ivar, rv, ebv, wavelengths=wavelengths)
        

      elif instr == "MILES":
        fi = fits.open(infile)
        head = fi[0].header
        s = fi[0].data
        flux = s[0,:]
        wav = np.arange(len(s[0,:]))*head['CDELTA1']+head['CRVAL1']
        lenwav = len(wav)        
        print('Warning: MILES files do not include uncertainties')
        print('         assuming S/N = 20!')
        err = flux * 0.05
        ivar4 = np.divide(1.,err**2, where = (err**2 > 0.) ,dtype = np.float128)
        ivar = np.float64(ivar4)


        if 'OBJECT' in head:
          ids = np.array([head['OBJECT']])
        else:
          ids = np.array([infile])

        assert(target is None),'target must be None for INT-IDS/MILES data (1 target per file)'

        xtr = (head)
        wav, frd, ivr = single_target_prep(wav, flux, ivar, rv, ebv, wavelengths=wavelengths)

      elif instr == "IDS-R900V":
        fi = fits.open(infile)
        head = fi[0].header
        s = fi[0].data
        if s.ndim == 1:
          flux = s
          print('Warning: the input IDS file does not include uncertainties')
          print('         assuming S/N = 20!')
          err = flux * 0.05
        elif s.ndim == 3:
          flux = s[0,0,:]
          err = s[3,0,:]
        else:
          print('Error: we cannot hundle this type of IDS file')
          sys.exit(1)
        wav = np.arange(len(flux))*head['CD1_1']+head['CRVAL1']
        lenwav = len(wav)        
        ivar4 = np.divide(1.,err**2, where = (err**2 > 0.) ,dtype = np.float128)
        ivar = np.float64(ivar4)

        if 'OBJECT' in head:
          ids = np.array([head['OBJECT']])
        else:
          ids = np.array([infile])

        assert(target is None),'target must be None for INT-IDS data (1 target per file)'

        xtr = (head)
        wav, frd, ivr = single_target_prep(wav, flux, ivar, rv, ebv, wavelengths=wavelengths)

      elif instr == "OSIRIS-R2500U":
        fi = fits.open(infile)
        head = fi[0].header
        s = fi[0].data
        flux = s
        wav = np.arange(len(s))*head['CD1_1']+head['CRVAL1']
        lenwav = len(wav)        
        print('Warning: OSIRIS files do not include uncertainties')
        print('         assuming S/N = 20!')
        err = flux * 0.05
        ivar4 = np.divide(1.,err**2, where = (err**2 > 0.) ,dtype = np.float128)
        ivar = np.float64(ivar4)

        if 'OBJECT' in head:
          ids = np.array([head['OBJECT']])
        else:
          ids = np.array([infile])

        assert(target is None),'target must be None for GTC-OSIRIS data (1 target per file)'
        
        xtr = (head)
        wav, frd, ivr =  single_target_prep(wav, flux, ivar, rv, ebv, wavelengths=wavelengths)

                  
    else: 

      assert(target is None),'target must be none for FERRE input files'
      assert(rv is None),'rv must be none for FERRE input files'
      instr = 'FERRE'
      if infile is None: infile = synthfile[2:synthfile.find('.dat')]
      frdfile = infile + '.frd'
      errfile = infile + '.err'
              
      #read frdfile and errfile
      frd = np.loadtxt(frdfile,dtype=float)
      err = (np.loadtxt(errfile,dtype=float)**2)
      ivr = np.divide(1.,err, where = (err > 0.) )
      wav = wavelengths
      if frd.ndim == 1:
        ids = np.array([0])
      else:
        ids = np.array(list(map(str,range(len(frd[:,0])))))
      xtr = (dict())

    return(ids,wav,frd,ivr,xtr)


def single_target_prep(wav,flux,ivar,rv,ebv,wavelengths=None):

    from extinction import apply,remove,ccm89

    if rv is None:
      vrad = 0.0
    else:
      try:
        _ = (e for e in rv)
        assert (len(rv) == 1),'rv must have a single value for LAMOST files with a single spectrum'
        vrad = float(rv[0])
      except TypeError:
        vrad = float(rv)

    if ebv is None:
      red = 0.0
    else:
      try:
        _ = (e for e in ebv)
        assert (len(ebv) == 1),'ebv must be have a single value for LAMOST files with a single spectrum'
        red = float(ebv[0])
      except TypeError:
        red = float(ebv)

    if np.abs(red) > 1e-7:
      xtmp = np.array(wav,dtype=float)
      ytmp = np.array(flux,dtype=float)
      tmp = remove(ccm89(xtmp, red * 3.1, 3.1), ytmp)
      ivar = ivar * (np.divide(ytmp,tmp,where=tmp>0))**2
      flux = tmp

    if wavelengths is None:
      frd = np.interp(wav,wav*(1. + vrad/clight),flux)
      ivr = np.interp(wav,wav*(1. + vrad/clight),ivar)
    else:
      lenx = len(wavelengths)
      frd = np.interp(wavelengths,wav,flux)
      ivr = np.interp(wavelengths,wav,ivar)
      wav = wavelengths

    return(wav,frd,ivr)


def read_desispec(filename,band=None):
  """Reads a DESI band spectrum, or a (full) SDSS/BOSS spectrum
  
  Parameters
  ----------
  filename: str
    name of an input DESI (spectra* or coadd*) or SDSS-BOSS (spPlate)
    FITS file
    
  band: str
    name of a band for multi-band DESI spectra
    (can be 'B', 'R' or 'Z')
    
  Returns
  -------
  wavelength: numpy array of floats
  
  flux: numpy array of floats
    observed fluxes
  
  ivar: numpy array of floats
    inverse variance for the observed fluxes
    
  res: structure
    resolution matrix for DESI
    
  header: dict
    very first header of the file

  fibermap: structure  (.data from the FIBERMAP extension)
    fibermap 
    
  scores: structure  (.data from the SCORES extension)
    scores
    
    
  """

  hdu=fits.open(filename)

  if filename.find('spectra-') > -1 or filename.find('exp_') > -1 or filename.find('coadd') > -1: #DESI
    header=hdu[0].header
    wavelength=hdu[band+'_WAVELENGTH'].data #wavelength array
    flux=hdu[band+'_FLUX'].data       #flux array (multiple spectra)
    ivar=hdu[band+'_IVAR'].data       #inverse variance (multiple spectra)
    #mask=hdu[band+'_MASK'].data       #mask (multiple spectra)
    res=hdu[band+'_RESOLUTION'].data  #resolution matrix (multiple spectra)
    #bintable=hdu['BINTABLE'].data  #bintable with info (incl. mag, ra_obs, dec_obs)
    fibermap=hdu['FIBERMAP'].data
    scores=hdu['SCORES'].data

  if filename.find('spPlate') > -1: #SDSS/BOSS
    header=hdu['PRIMARY'].header
    wavelength=header['CRVAL1']+arange(header['NAXIS1'])*header['CD1_1']
#wavelength array
    wavelength=10.**wavelength
    flux=hdu['PRIMARY'].data       #flux array (multiple spectra)
    #ivar=hdu['IVAR'].data       #inverse variance (multiple spectra)
    ivar=hdu[1].data       #inverse variance (multiple spectra)
    #andmask=hdu['ANDMASK'].data       #AND mask (multiple spectra)
    #ormask=hdu['ORMASK'].data       #OR mask (multiple spectra)
    #res=hdu['WAVEDISP'].data  #FWHM array (multiple spectra)
    res=hdu[4].data  #FWHM array (multiple spectra)
    #bintable=hdu['BINTABLE'].data  #bintable with info (incl. mag, ra, dec)
    fibermap=hdu['FIBERMAP'].data
    scores=hdu['FIBERMAP'].data

  return((wavelength,flux,ivar,res,header,fibermap,scores))


def plot_spec(root=None,x=None,n=None,m=None,o=None,xrange=None,yrange=None,nozero=None,res=False):

  """Plot one or multiple spectra
  """

  if root is not None:
    xx = np.loadtxt(root+'.wav') / 10.
    n = np.loadtxt(root+'.nrd')
    m = np.loadtxt(root+'.mdl')
    o = np.loadtxt(root+'.opf',dtype=str) 

  else:

    if type(x) is list: 
      xx = np.hstack(x)
      xx = xx/10.
      p = [0]
      for i in range(len(x)):
        p.append(p[i] + len(x[i])) 
    else:
      xx = x/10.

  if xrange is None: 
    xrange = ( np.min(xx), np.max(xx) )

  if xx.ndim == 1:
    nfreq = len(xx)
  else:
    nfreq = len(xx[0,:])
  if n.ndim == 1: 
    nspec = 1
    plt.clf()
    labels = []
    if type(x) is list:
      for i in range(len(x)):
        if nozero:
          w = (n[p[i]:p[i+1]] > 0.)
        else:
          w = range(p[i+1]-p[i])
        plt.plot(xx[p[i]:p[i+1]][w],n[p[i]:p[i+1]][w])
        labels.append('data')
        if m is not None: 
          plt.plot(xx[p[i]:p[i+1]][w],m[p[i]:p[i+1]][w])
          labels.append('model')
        if res:
          plt.plot(xx[p[i]:p[i+1]][w],m[p[i]:p[i+1]][w]-n[p[i]:p[i+1]][w])
          labels.append('residuals')
    else:
      if nozero:
        w = (n > 0.)
      else:
        w = range(len(xx))
      plt.plot(xx[w],n[w])
      labels.append('data')
      if m is not None:
        plt.plot(xx[w],m[w])
        labels.append('model')
      if res:
        plt.plot(xx[w],m[w]-n[w])
        labels.append('residuals')
    if yrange is None: yrange = [np.min(n)*0.95,np.max(n)*1.05]
    if res: yrange[0] = np.min(m[w]-n[w])*1.05
    plt.xlabel('wavelength (nm)')
    plt.ylabel('normalized flux')
    if o is not None:
      npar = len(o)-3
      if npar >= 16: 
        npar = int(np.sqrt(npar*1.0 - 1.))
      else:
        npar = npar// 2
      plt.title('params: '+' -- '.join(map("{:.2f}".format,np.array(o[1:npar+1],dtype=float))))
      xtext = 0.5*xrange[0]+0.5*xrange[1]
      ytext = 0.75*yrange[0]+0.25*yrange[1]
      ycurve = np.interp(xtext,xx[w],n[w])
      if abs(ycurve-ytext) < 0.5: ytext = ytext*1.2
      plt.text(xtext,ytext,o[0])
    plt.xlim(xrange)
    plt.ylim(yrange)
    if m is not None: plt.legend(labels)
    plt.savefig('fig1.png')
    plt.show()
  else:
    nspec = len(n[:,0])
    labels = []
    for j in range(nspec):
      plt.clf()
      if type(x) is list:
        for i in range(len(x)):
          if nozero:
            w = (n[j,p[i]:p[i+1]] > 0.)
          else:
            w = np.ones(p[i+1]-p[i],dtype=bool)
          plt.plot(xx[p[i]:p[i+1]][w],n[j,p[i]:p[i+1]][w])
          labels.append('data')
          if m is not None:
            plt.plot(xx[p[i]:p[i+1]][w],m[j,p[i]:p[i+1]][w])
            labels.append('model')
          if res:
            plt.plot(xx[p[i]:p[i+1]][w],m[j,p[i]:p[i+1]][w]-n[j,p[i]:p[i+1]][w])
            labels.append('residuals')
      else:
        if nozero:
          w = (n[j,:] > 0.)
        else:
          w = np.ones(nfreq,dtype=bool)
      
        if xx.ndim  == 1:
          xx2 = xx[w]
        else:
          xx2 = xx[j,w].transpose()
        plt.plot(xx2,n[j,w])
        labels.append('data')
        if m is not None:
          plt.plot(xx2,m[j,w])
          labels.append('model')
        if res:
          plt.plot(xx2,m[j,w]-n[j,w])
          labels.append('residuals')
      if yrange is None: 
        yrange2 = [np.min(n[j,:])*0.95,np.max(n[j,:])*1.05]
      else:
        yrange2 = yrange
      if res: yrange2[0] = np.min(m[j,w]-n[j,w])*1.05
      plt.xlabel('wavelength (nm)')
      plt.ylabel('normalized flux')
      plt.xlim(xrange)
      plt.ylim(yrange2)
      if m is not None: plt.legend(labels)
      if o is not None:
        npar = len(o[0,:])-3
        if npar >= 16: 
          npar = int(np.sqrt(npar*1.0 - 1.))
        else:
          npar = npar// 2
        plt.title('params: '+' -- '.join(map("{:.2f}".format,np.array(o[j,1:npar+1],dtype=float))))
        xtext = 0.5*xrange[0]+0.5*xrange[1]
        ytext = 0.75*yrange2[0]+0.25*yrange2[1]
        ycurve = np.interp(xtext,xx2,n[j,w])
        if abs(ycurve-ytext) < 0.5: ytext = 0.25*yrange2[0]+0.75*yrange2[1]
        plt.text(xtext,ytext,o[j,0])


      plt.savefig('fig'+str(j+1)+'.png')
        

  return()


def vac2air(wavelength):
    """Conversion from vacuum to air for wavelengths based on Ciddor (1996)
    
    Parameters
    ----------
    wavelength: float, or iterable
       vacuum wavelength(s) in AA
    
    Returns
    -------
    wavelength: numpy array
       corresponding air wavelengths in AA (keeps vacuum for lambda<=2000. A)
       
    Based on Ciddor (1996). Copied literally from the IDL Astro library
    """
    if type(wavelength) is list or type(wavelength) is float or type(wavelength) is tuple:
        wavelength = np.array(wavelength)
    g = (wavelength >= 2000.)
    sigma2 = (1e4/wavelength[g])**2
    fact = 1. +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/( 57.362 - sigma2)
    wavelength[g] = wavelength[g]/fact

    return(wavelength)
    
def xc(y1,y2, npoints=None, gaus=False, plot=False):
    """Determines a pixel offset between two arrays by 
    cross correlation.
    
    Parameters
    ----------
    y1 float array
       Array 1

    y2 float array
       Array 2    
       
    npoints int
       Number of points to use in fitting a model to the peak
       of the cross-correlation function
    
    gaus bool
       When True a Gaussian model is used instead of the default
       parabola
    
    plot bool
       When True the program plots the model fit to the peak of
       the cross-correlation function
    
    Returns
    -------
    delta float
       Shift to apply to y2 to match y1 (max of the cross-correlation
       function)
       (units are in pixels)
       
    edelta float
       Uncertainty in delta (pixels)
    
    """
    
    #set failure values for delta/edetla
    delta = None
    e_delta = None
    
    nel = len(y1)
    
    #basic checks
    assert nel == len(y2), 'error: the two vectors have different dimensions'
	
	#compute the cross-correlation
	
    #xlen = len(y1) + len(y2) - 1
    #x = np.arange(xlen)
    #ccf = np.correlate(y1,y2,mode='full')
    #w = np.where(ccf == np.max(ccf))[0]
        
    nrange = nel // 2 - 1
    ccf = np.zeros(2*nrange + 1)
    for i in np.arange(2*nrange + 1):
        yy = np.roll(y2,i-nrange)
        ccf[i] = np.sum(y1*yy)
                	
    x = np.arange(len(ccf),dtype=float)-nrange
    w = np.where(ccf == np.max(ccf))[0]
    w0 = w[0]
    
    if (gaus): 
        if npoints is None: npoints = 39
    else:
        if npoints is None: npoints = 7

    assert npoints <= nel, 'error: npoints > number of elements of the input array'
    
    #finding the central pixel from weighted average in a symmetric window
    nhalf = (npoints - 1) // 2
    w2 =  np.sum( ccf[w0-nhalf:w0+nhalf+1] * x[w0-nhalf:w0+nhalf+1] ) / np.sum( ccf[w0-nhalf:w0+nhalf+1] ) + nrange
    w = int(np.rint(w2))  


    #dealing with even values of npoints
    nhalf1 = nhalf
    nhalf2 = nhalf
    if (npoints // 2 == int(np.rint(npoints / 2))):
        if w2 > w:
            nhalf2 = nhalf2 + 1
        else:
            nhalf1 = nhalf1 + 1

    assert nhalf1 + nhalf2 + 1 == npoints, 'error: something is wrong!' 
    assert ((w-nhalf1 >= 0) and (w+nhalf2 <= len(x)-1)),'error: not enough points to fit'
    
    xx = x[w-nhalf1:w+nhalf2+1]
    yy = ccf[w-nhalf1:w+nhalf2+1]		

    if (gaus):
        p0 = [np.mean(yy)*2., x[w0], 2.0, np.min(yy)]
        coef, covar = curve_fit(gauss, xx, yy, p0=p0)
        delta = coef[1]
        edelta = covar[1,1]
        edelta = np.sqrt(edelta)
        model = gauss(xx,coef[0],coef[1],coef[2],coef[3])
    else:
        coef, covar = np.polyfit(xx, yy, 2, cov=True)
        delta = -coef[1]/2./coef[0]
        edelta = 1./4./coef[0]**2*(covar[1,1] + coef[1]**2/coef[0]**2*covar[0,0]) - coef[1]/2./coef[0]**3*covar[1,0]		
        edelta = np.sqrt(edelta)
        model = coef[0]*xx**2+coef[1]*xx+coef[2]

    if plot:
        plt.plot(xx,yy,'.')
        plt.plot(xx,model)
        plt.show()


    return (delta,edelta)


def xc_np(y1,y2, npoints=None, gaus=False, plot=False):
    """Determines a pixel offset between two arrays by 
    cross correlation. Variant of xc using numpy correlate 
    
    Parameters
    ----------
    y1 float array
       Array 1

    y2 float array
       Array 2    
       
    npoints int
       Number of points to use in fitting a model to the peak
       of the cross-correlation function
    
    gaus bool
       When True a Gaussian model is used instead of the default
       parabola
    
    plot bool
       When True the program plots the model fit to the peak of
       the cross-correlation function
    
    Returns
    -------
    delta float
       Shift to apply to y2 to match y1 (max of the cross-correlation
       function)
       (units are in pixels)
       
    edelta float
       Uncertainty in delta (pixels)
    
    """

    #set failure values for delta/edelta
    delta = None
    edelta = None
	
    xlen = len(y1) + len(y2) - 1
    x = np.arange(xlen, dtype=float) - len(y1) 
    ccf = np.correlate(y1,y2,mode='full')
	
	
    w = np.where(ccf == np.max(ccf))[0]
    w0 = w[0]
    
    if (gaus):
        if npoints is None: npoints = 39
        nhalf = (npoints - 1) // 2
        xx = x[w0-nhalf:w0+nhalf+1]
        yy = ccf[w0-nhalf:w0+nhalf+1]
        p0 = [np.mean(yy)*2., x[w0], 2.0, np.min(yy)]
        coef, covar = curve_fit(gauss, xx, yy, p0=p0)
        delta = coef[1]
        edelta = covar[1,1]
        edelta = np.sqrt(edelta)
        model = gauss(xx,coef[0],coef[1],coef[2],coef[3])
    else:
        if npoints is None: npoints = 7
        nhalf = npoints // 2
        xx = x[w0-nhalf:w0+nhalf+1]
        yy = ccf[w0-nhalf:w0+nhalf+1]		
        coef, covar = np.polyfit(xx, yy, 2, cov=True)
        delta = -coef[1]/2./coef[0]
        edelta = 1./4./coef[0]**2*(covar[1,1] + coef[1]**2/coef[0]**2*covar[0,0]) - coef[1]/2./coef[0]**3*covar[1,0]		
        edelta = np.sqrt(edelta)
        model = coef[0]*xx**2+coef[1]*xx+coef[2]

    if plot:
        plt.plot(xx,yy,'.')
        plt.plot(xx,model)
        plt.show()


    return (delta,edelta)
        
def gauss(x, *p):
    """Evaluate a Gaussian function from the input parameters
    """	
    A, mu, sigma, base = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + base
    
def continuum(x, window_length = 500, polyorder = 3):
	
    """Smoothing the data in the 1D array x using a 
    Saviztky-Golay filter

    parameters
    ----------
    
x : 1D array
    The data to be filtered. If `x` is not a single or double precision
    floating point array, it will be converted to type ``numpy.float64``
    before filtering.
window_length : int
    The length of the filter window (i.e., the number of coefficients).
    If `mode` is 'interp', `window_length` must be less than or equal
    to the size of `x`.
    (default 50)
polyorder : int
    The order of the polynomial used to fit the samples.
    `polyorder` must be less than `window_length`.
    (default 3)
    """    
    
    return(savgol_filter(x, window_length, polyorder))	


def xxc(x1,y1,iva1,x2,y2, maxv = 1000., plot=False):
	
    """Determines the velocity offset to apply to the template
    spectrum (x2,y2) so that it overlaps with the observed 
    spectrum (x1,y1,iva1)
    
    Parameters
    ----------
    x1: numpy array of floats
      wavelenghts for the observed spectrum  (A)
    y1: numpy array of floats
      fluxes for the observed spectrum
    iva1: numpy array of floats
      inverse variance for the fluxes of the observed spectrum
    x2: numpy array of floats
      wavelenghts for the model/template spectrum
    y2: numpy array of floats
      fluxes for the model/template spectrum
    maxv: float
      maximum velocity offset to consider between the observed
      and template spectra
      (default 1000 km/s)
    plot: bool
      switch to show a plot of the chi-square as a function of 
      velocity offset

    Returns
    -------
    delta: float
      velocity offset need to overlap the template on top of the
      observed (y1) spectrum in units of km/s
    edelta: float
      uncertainty in delta
    """

    lenx1 = len(x1)
    dv = np.diff(x1).mean()/x1.mean()*clight / 10.
    nv = int(2*maxv/dv)
    v = np.arange(nv)*dv - maxv
    lenv = len(v)
    chi = np.zeros(lenv)
    for i in range(lenv):
        yy2 = np.interp(x1,x2*(1.+v[i]/clight),y2)
        chi[i] = np.sum( ( y1 - yy2 )**2 * iva1 ) 
  
    beta = np.median(chi) / 1490. / 5. 
    while np.exp(-np.min(chi)/2./beta) <= 0.0:
        beta = beta * 2.

    likeli = np.exp(-chi/2./beta)
    den = np.sum(likeli)
    delta = np.sum(v*likeli) / den
    edelta = np.sum( ( v - delta)**2 * likeli) / den
    edelta = np.sqrt(edelta)
    
    if plot:
        plt.plot(v,chi)
        plt.show()
  
    return(delta,edelta)
 
 
def bas_build(synthfile):

    conf = load_conf(config='bas-build.yaml',confdir=confdir)

    ending = synthfile.find('.dat') 
    if ending < -1: 
        ending = synthfile.find('.pickle')
    if ending < -1: 
        ending = len(synthfile) + 1
    root = synthfile[2:ending]


    for entry in conf.keys():
        for task in conf[entry]:
            print(entry,task)
            for job in conf[entry][task]:
                if 'synthfile' not in job:
                  job['synthfile'] = synthfile
                if task == 'pickle_synth':
                  if 'synthfile' in job:
                    pos = job['synthfile'].find('_')
                    job['outsynthfile'] = job['synthfile'][:pos] + '_' + \
                      root + '-' + job['synthfile'][pos+1:-3] + 'pickle'
                print('kargs=',job)
                call(task,**job)
    
    return()

def call(func,**kargs):

    func_to_run = globals()[func]
    func_to_run(**kargs)

    return()


def rewrite_synth(synthfile,outsynthfile=None):

    """Read a FERRE grid and write it back to disk as an irregular grid
    """

    if outsynthfile is None:
        outsynthfile = 'i_'+synthfile[2:]

    h,p,d = read_synth(synthfile)
    write_synth(outsynthfile, p, d, hdr=h, irregular=True)


    return()



def bas_perfcheck(synthfile,n=1000,snr=1.e6,
    kernel='thin_plate_spline', neighbors=100, focus=False, edgemargin=0.05):

    """Carry out a full performance check using bas on a synthetic grid

    Parameters
    ---------
    synthfile: str
      Name of the input synth file
    n: int
      Number of mock spectra to produce for the test
      (default is 1000)
    snr: float
      Signal to noise ratio for the mock spectra
      (default is 1.e6)
    kernel: string
      Type of RBF function (linear, thin_plate_spline, cubic, gaussian ...)
    neighbors: int
       Number of nearest neighbors used to compute the interpolation
       coefficients for each grid point
    focus: bool
      switch to activate a two-step algorithm in which a coarsely 
      subsampled version of the grid is used to identify first where   
      the optimal solution is, and then perform an focused analysis
      in that region (a +/- 3 sigma volume) 
   edgemargin: float
      fraction from the min/max values of the input parameters to exclude
      from the range to sample
      (default 0.05 -- exclude 5% from the edges)


    Returns
    ------
    result: numpy array of floats
      16-50-86 percentiles for all the parameters
    
    """
    hd = head_synth(synthfile)
    if 'NTOT' in hd: 
      ntot = int(hd['NTOT'])
    else:
      if 'N_P' in hd:
         n_p = map(int,hd['N_P']) 
         ntot = np.product(n_p)
      else:
        ntot = 0

    checksynthfile=synthfile+'-check.dat'
    synth_rbf(synthfile,outsynthfile=checksynthfile,n=n,
              rv=False,ebv=False,kernel=kernel,neighbors=neighbors, 
              edgemargin=edgemargin)
    bas_test(checksynthfile,snr=snr)
    print('running ... ','bas(',checksynthfile[2:-4],'synthfile=',synthfile,')')
    now = time.time()
    bas(checksynthfile[2:-4],synthfile=synthfile,focus=focus)
    print('this run took ',time.time()-now,' seconds')
    result = fparams(checksynthfile[2:-4],synthfile=synthfile,
                     figure=checksynthfile[2:-4]+'-n'+str(n)+'-snr'+str(snr)+'.png')


    fh = open('-'.join((synthfile,str(n),kernel,str(neighbors),'bas_perfcheck.dat')),'w')
    fh.write(str(n)+' '+str(ntot)+' '+' '.join(map(str,np.concatenate(result)))+'\n')
    fh.close()
    

  
    return(result) 
     


def bas_test(synthfile,snr=1.e6):
    """Use the data in synthfile to create mock observations
	   for testing purposes
    """
    h,p,d = read_synth(synthfile)    
    ndim = len(p[0,:])
    npix = len(d[0,:])
    ending = synthfile.rfind('.dat')
    if ending < -1: 
        ending = synthfile.rfind('.pickle')
    if ending < -1:
        ending = len(synthfile) + 1
    root = synthfile[2:ending]
    vf = open(root+'.ipf','w')
    of = open(root+'.frd','w')
    ef = open(root+'.err','w')
    
    for i in range(len(p[:,0])):
        vf.write(str(i)+' '+' '.join(map(str,p[i,:]))+'\n')
        of.write(' '.join(map(str,d[i,:] * (1. + 
             1./snr*np.random.normal(size=npix))))+'\n')
        ef.write(' '.join(map(str,d[i,:]/snr))+'\n')
	
    vf.close()
    of.close()
    ef.close()
    
    return
        

def synth_rbf(synthfile,outsynthfile=None,n=None,rv=False,ebv=False,
              kernel='thin_plate_spline', neighbors=100, edgemargin=0.0):

    """Creates an irregular FERRE grid from a pre-existing regular or 
       irregular one
 
    Parameters
    ----------
    synthfile: str
      name of the input FERRE/BAS synthfile
    
    outsynthfile: str
      name of the output FERRE/BAS synthfile
    
    n: int
      number of spectra to create by RBF interpolation
    
    rv: bool
      if true, fold in a dimension with RV variations
      
    ebv: bool
      if true, fold in a dimension with E(B-V) variations

    kernel: string
      Type of RBF function (linear, thin_plate_spline, cubic, gaussian ...)

    neighbors: int
      Number of nearest neighbors used to compute the interpolation coefficients
      for each grid point

    edgemargin: float
      fraction from the min/max values of the input parameters to exclude
      from the range to sample
      (default 0.0 -- use the full range) 
       

      
    Returns
    -------
    Creates an output grid (outsynthfile)

    """
    
    from extinction import apply,ccm89

    if rv or ebv : 
      x = lambda_synth(synthfile)
      if type(x) is list: x = np.hstack(x)
    h,p,d = read_synth(synthfile)
    
    
    ndim = len(p[0,:])
    npix = len(d[0,:])
    ntot = len(p[:,0])

    if n is None: n = ntot
    
    for i in range(ndim):

        amin = np.min(p[:,i])*(1.0+edgemargin)
        amax = np.max(p[:,i])*(1.0-edgemargin)
        vals = np.random.random_sample(n)*(amax-amin) + amin
        if i == 0: 
            p2 = vals
        else:	
            p2 = np.vstack((p2,vals))
    
    c, pmin, ptp = rbf_get(synthfile, kernel=kernel, neighbors=neighbors)
    d2 = rbf_apply(c, pmin, ptp, np.transpose(p2))
    h2 = h

    
    ndim2 = 0
    if rv:		
        if type(h2) is list: 
          h0 = h2[1]
        else:
          h0 = h2  
        rvmax = 1000.
        if 'RESOLUTION' in h0: rvmax = clight/float(h0['RESOLUTION'])
        vals = np.random.random_sample(n)*2*rvmax - rvmax
        p2 = np.vstack((p2,vals))
        ndim2 += 1
        if type(h2) is list:
            for entry in range(len(h2)):
                h2[entry]['LABEL('+str(ndim+ndim2)+')'] = "'RV'"
        else:
            h2['LABEL('+str(ndim+ndim2)+')'] = "'RV'"  
    if ebv:	
        ebvmax = 0.25 # mag
        vals = np.random.random_sample(n)*ebvmax 
        p2 = np.vstack((p2,vals))
        ndim2 += 1
        if type(h2) is list:
           for entry in range(len(h2)):
               h2[entry]['LABEL('+str(ndim+ndim2)+')'] = "'E(B-V)'"
        else:
            h2['LABEL('+str(ndim+ndim2)+')'] = "'E(B-V)'"
        
    p2 = np.transpose(p2)
    ending = synthfile.rfind('.dat')
    if ending < -1: 
        ending = synthfile.rfind('.pickle')
    if ending < -1:
        ending = len(synthfile) + 1
    root = synthfile[2:ending]
    if outsynthfile is None: 
        outsynthfile = 'n_'+root+'rbf'
        if rv: outsynthfile += '-RV'
        if ebv: outsynthfile += '-EBV' 
        outsynthfile += '.dat'
    of = open(outsynthfile,'w')
    print('writing grid ',outsynthfile,'...')
    if type(h2) is not list:
      h2 = [h2]
    for block in h2:
      block['TYPE'] = "'irregular'"
      block['N_OF_DIM'] = str(ndim+ndim2)
      block['NTOT'] = str(n)
    for block in h2:
      of.write(' &SYNTH\n')
      for entry in block: of.write(' '+entry + ' = ' + block[entry] + '\n')
      of.write(' /\n')

    for i in range(len(p2[:,0])):
        flx = d2[i,:]
        if rv:
            #print(type(x),len(x),type(flx),len(flx))
            flx = np.interp(x,x*(1.+p2[i,ndim]/clight),flx)
        if ebv:
            flx = apply(ccm89(x, p2[i,ndim+ndim2-1]* 3.1, 3.1), flx)
        
        #normalization
        #flx = flx / np.mean(flx)
        of.write(' '.join(map(str,p2[i,:]))+' '+' '.join(map(str,flx))+'\n')	
    of.close()
    
    return

def rbf_test(synthfile,n=None, kernel='thin_plate_spline', neighbors=100):
    """Creates an irregular FERRE grid using RBF interpolation on a 
       pre-existing regular or irregular one, and from that one the
       interpolation is repeated to return to the original/input grid
       and estimate interpolation errors
 
    Parameters
    ----------
    synthfile: str
      name of the input FERRE/BAS synthfile
    
    n: int
      number of spectra to create by RBF interpolation
    (default is the number of spectra in the input grid)

    kernel: string
      Type of RBF function (linear, thin_plate_spline, cubic, gaussian ...)

    neighbors: int
       Number of nearest neighbors used to compute the interpolation 
       coefficients for each grid point
     
 
    Returns
    -------
    tuple of 3 floats
      16%, 50% and 84% percentiles between the original input and the output. 
    
    """
    

    h,p,d = read_synth(synthfile)    
    
    ndim = len(p[0,:])
    npix = len(d[0,:])
    ntot = len(p[:,0])

    if n is None: n = ntot
    
    #1st interpolation
    synth_rbf(synthfile,outsynthfile=synthfile+'-tmp',n=n,
              kernel=kernel, neighbors=neighbors)

    h2,p2,d2 = read_synth(synthfile+'-tmp')
    
    #2nd interpolation    
    c, pmin, ptp = rbf_get(synthfile+'-tmp',
                           kernel=kernel, neighbors=neighbors)
    d2 = rbf_apply(c, pmin, ptp, p)

    per = np.percentile( (d2-d)/d,[15.85,50.,84.15])

    fh = open('-'.join((synthfile,str(n),kernel,str(neighbors),'rbf_test.dat')),'w')
    fh.write(' '.join(map(str,(n, per[0], per[1], per[2])))+'\n')
    fh.close()
        
    return( per[0], per[1], per[2] )
    

def wferrefits(root, path=None):	
	
  """Packs FERRE output into a FITS file
  
  
  Parameters
  ----------
  root: str
      name of the root for input/output FERRE/BAS files 
      with extensions .opf, .wav, .nrd, .mdl

  path: string
      path to files
      (default is None, and the code looks for the FERRE files 
      in the current folder)
 
  Returns
  -------  
      A FITS files with the data, with the same root
      and a fits extension
      
  """

  if path is None: path=""
  proot=os.path.join(path,root)
  o=glob.glob(proot+".opf")

  
  xbandfiles = sorted(glob.glob(proot+'-*.wav'))
  band = []
  npix = []
  for entry in xbandfiles:
    print('entry=',entry)
    match = re.search('-[\w]*.wav',entry)
    tag = match.group()[1:-4]
    if match: band.append(tag.upper())
    x = loadtxt(proot+'-'+tag+'.wav')
    npix.append(len(x))
  
  print('proot+.wav=',proot+'.wav')  
  print('xbandfiles=',xbandfiles)
  x = np.loadtxt(proot+'.wav')
  if len(npix) == 0: npix.append(len(x))

  mdl=glob.glob(proot+".mdl")
  err=glob.glob(proot+".err")
  nrd=glob.glob(proot+".nrd")

 
  success=[]
  targetid=[]
  srcfile=[]
  teff=[]
  logg=[]
  feh=[]
  alphafe=[]
  micro=[]
  param = []
  covar=[]
  snr_med=[]
  chisq_tot=[]
  rv_adop=[]
  rv_err=[]
  bestgrid = []
  of=open(o[0],'r')


  for line in of:
    cells=line.split()
    bestgrid.append(o[0])
    id = cells[0]
    cells = cells[1:]

    m = len(cells)
    ndim = m - 3
    if ndim > 10:
        ndim = int(np.sqrt(4*ndim+1)-1)
    ndim = ndim // 2
    
    assert (m > 6), 'Error, the file '+o[0]+' has less than 7 columns, which would correspond to ndim=2'
    par = np.zeros(ndim)
    cov = np.zeros(ndim*ndim+ndim) 
    
    print('ndim=',ndim)
    print('m=',m)
    
    if (ndim == 2):
      #white dwarfs 2 dimensions: id, 2 par, 2err, 0., med_snr, lchi, 2x2 cov
      feh.append(-10.)
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(np.nan)
      micro.append(np.nan)


    elif (ndim == 3):
      #Kurucz grids with 3 dimensions: id, 3 par, 3 err, 0., med_snr, lchi, 3x3 cov
      #see Allende Prieto et al. (2018, A&A)
      feh.append(float(cells[2]))
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(np.nan)
      micro.append(np.nan)

    elif (ndim == 4):
      #Phoenix grid from Sergey or MARCS grid, with 4 dimensions: id, 4 par, 4err, 0., med_snr, lchi, 4x4 cov
      feh.append(float(cells[1]))
      teff.append(float(cells[2]))
      logg.append(float(cells[3]))
      alphafe.append(float(cells[0]))
      micro.append(np.nan)
   

    elif (ndim == 5):
      #Kurucz grids with 5 dimensions: id, 5 par, 5 err, 0., med_snr, lchi, 5x5 cov
      #see Allende Prieto et al. (2018, A&A)
      feh.append(float(cells[0]))
      teff.append(float(cells[3]))
      logg.append(float(cells[4]))
      alphafe.append(float(cells[1]))
      micro.append(float(cells[2]))


    chisq_tot.append(10.**float(cells[2+2*ndim]))
    snr_med.append(float(cells[1+2*ndim]))
    rv_adop.append(float(cells[0+2*ndim]))
    rv_err.append(np.nan)
    par = np.array(cells[0:ndim],dtype=float)
    cov = np.array(cells[2*ndim:],dtype=float)
    param.append(par)
    covar.append(cov)    
     


    if (chisq_tot[-1] < 1.5 and snr_med[-1] > 5.): # chi**2<1.5 and S/N>5
      success.append(1) 
    else: success.append(0)
    targetid.append(id)
    srcfile.append(root)
    #fiber.append(int32(tmp[1]))

  #primary extension
  hdu0=fits.PrimaryHDU()

  #find out processing date and add it to primary header
  now = datetime.datetime.fromtimestamp(time.time())
  nowstr = now.isoformat() 
  nowstr = nowstr[:nowstr.rfind('.')]
  hdu0.header['DATE'] = nowstr

  #find out host machine and add info to header
  try:
    host=os.environ['HOST']
  except:
    host='Unknown'
  hdu0.header['HOST'] = host
  #find out OS name/platform
  osname = os.name 
  platf = platform.system() + ' '+ platform.release()
  hdu0.header['OS'] = osname
  hdu0.header['PLATFORM'] = platf

  #keep track of the number of targets processed and the time it took
  nspec = len(targetid)
  hdu0.header['NSPEC'] = nspec

  
  hdulist = [hdu0]

  #sptab extension
  cols = {}
  cols['SUCCESS'] = success
  cols['TARGETID'] = targetid
  cols['SRCFILE'] = srcfile
  cols['BESTGRID'] = bestgrid
  cols['TEFF'] = np.array(teff)*units.K
  cols['LOGG'] = np.array(logg)
  cols['FEH'] = np.array(feh)
  cols['ALPHAFE'] = np.array(alphafe) 
  cols['LOG10MICRO'] = np.array(micro)
  cols['PARAM'] = np.array(param)
  cols['COVAR'] = np.array(covar)
  cols['CHISQ_TOT'] = np.array(chisq_tot)
  cols['SNR_MED'] = np.array(snr_med)
  cols['RV_ADOP'] = np.array(rv_adop)*units.km/units.s
  cols['RV_ERR'] = np.array(rv_err)*units.km/units.s

  colcomm = {
  'success': 'Bit indicating whether the code has likely produced useful results',
  'TARGETID': 'DESI targetid',
  'SRCFILE': 'DESI data file',
  'BESTGRID': 'Model grid that produced the best fit',
  'TEFF': 'Effective temperature (K)',
  'LOGG': 'Surface gravity (g in cm/s**2)',
  'FEH': 'Metallicity [Fe/H] = log10(N(Fe)/N(H)) - log10(N(Fe)/N(H))sun' ,
  'ALPHAFE': 'Alpha-to-iron ratio [alpha/Fe]',
  'LOG10MICRO': 'Log10 of Microturbulence (km/s)',
  'PARAM': 'Array of atmospheric parameters ([Fe/H], [a/Fe], log10micro, Teff,logg)',
  'COVAR': 'Covariance matrix for ([Fe/H], [a/Fe], log10micro, Teff,logg)',
  'ELEM': 'Elemental abundance ratios to hydrogen [elem/H]',
  'ELEM_ERR': 'Uncertainties in the elemental abundance ratios',
  'CHISQ_TOT': 'Total chi**2',
  'SNR_MED': 'Median signal-to-ratio',
  'RV_ADOP': 'Adopted Radial Velocity (km/s)',
  'RV_ERR': 'Uncertainty in the adopted Radial Velocity (km/s)'
  }      

  
  table = tbl.Table(cols)
  hdu=fits.BinTableHDU(table,name = 'SPTAB')
  #hdu.header['EXTNAME']= ('SPTAB', 'Stellar Parameter Table')
  k = 0
  for entry in colcomm.keys():
    print(entry) 
    hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
    k+=1
  hdulist.append(hdu)

  #aux extension
  #p = ['[Fe/H]','[a/Fe]','log10micro','Teff','logg']
  #cols = {}
  #colcomm = {}
  #cols['p'] = [p]
  #colcomm['p'] = 'PARAM tags'
  
  #table = tbl.Table(cols)
  #hdu=fits.BinTableHDU(table,name = 'AUX')

  #k = 0
  #for entry in colcomm.keys():
  #  print(entry) 
  #  hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
  #  k+=1
  #hdulist.append(hdu)


  #hdul=fits.HDUList(hdulist)
  #hdul.writeto('sptab_'+root+'.fits')

  #now we handle the extensions with the fluxes  

  edata=np.loadtxt(err[0])
  if (len(mdl) > 0): 
    mdata=np.loadtxt(mdl[0])
  if (len(nrd) > 0): 
    odata=np.loadtxt(nrd[0])


  i = 0
  j1 = 0


  if len(xbandfiles) > 0:
    for entry in band:
      j2 = j1 + npix[i] 
      print(entry,i,npix[i],j1,j2)
      #colx = fits.Column(name='wavelength',format='e8', array=array(x[j1:j2]))
      #coldefs = fits.ColDefs([colx])
      #hdu = fits.BinTableHDU.from_columns(coldefs)
      hdu = fits.ImageHDU(name=entry+'_WAVELENGTH', data=x[j1:j2])
      hdulist.append(hdu)
      
      cols = {}
      colcomm = {}
      if odata.ndim == 2: tdata = odata[:,j1:j2]
      else: tdata = odata[j1:j2][None,:]
      cols['obs'] = tdata
      colcomm['obs'] = 'Observed spectra as fit'
      if edata.ndim == 2: tdata = edata[:,j1:j2]
      else: tdata = edata[j1:j2][None,:]
      cols['err'] = tdata
      colcomm['err'] = 'Error in spectra as fit'
      if (len(mdl) > 0): 
        if mdata.ndim == 2: tdata = mdata[:,j1:j2]
        else: tdata = mdata[j1:j2][None,:]
        cols['fit'] = tdata
        colcomm['fit'] = 'Best-fitting model'      
      

      table = tbl.Table(cols)
      hdu=fits.BinTableHDU(table,name = entry+'_MODEL')
      k = 0
      for entry in colcomm.keys():
        print(entry) 
        hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
        k+=1
      hdulist.append(hdu)
      i += 1
      j1 = j2      
      
  else:
    print('single wavelength entry!')
    #colx = fits.Column(name='wavelength',format='e8', array=array(x[j1:j2]))
    #coldefs = fits.ColDefs([colx])
    #hdu = fits.BinTableHDU.from_columns(coldefs)
    hdu = fits.ImageHDU(name='WAVELENGTH', data=x)
    hdulist.append(hdu)
      
    cols = {}
    colcomm = {}
    if odata.ndim == 2: tdata = odata[:,:]
    else: tdata = odata[:][None,:]
    cols['obs'] = tdata
    colcomm['obs'] = 'Observed spectra as fit'
    if edata.ndim == 2: tdata = edata[:,:]
    else: tdata = edata[:][None,:]
    cols['err'] = tdata
    colcomm['err'] = 'Error in spectra as fit'
    if (len(mdl) > 0): 
      if mdata.ndim == 2: tdata = mdata[:,:]
      else: tdata = mdata[:][None,:]
      cols['fit'] = tdata
      colcomm['fit'] = 'Best-fitting model'      
      

    table = tbl.Table(cols)
    hdu=fits.BinTableHDU(table,name = 'MODEL')
    k = 0
    for entry in colcomm.keys():
      print(entry) 
      hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
      k+=1
    hdulist.append(hdu)


  hdul=fits.HDUList(hdulist)
  hdul.writeto(root+'.fits')
  
  return None
  
def wtabmodfits(root, path=None):
	
  """Write out DESI MWS SP pipeline output

  Parameters
  ----------
  root: str
      name of the root for input/output FERRE/BAS files 
      with extensions .opf, .wav, .nrd, .mdl

  path: string
      path to files
      (default is None, and the code looks for the FERRE/BAS files 
      in the current folder)
 
  Returns
  -------  
      An SPTAB file with the parameters, and an SPMOD with the spectra      
  
  """
  
  if path is None: path=""
  proot=os.path.join(path,root)

  o=glob.glob(proot+".opf")

  xbandfiles = sorted(glob.glob(proot+'-*.wav'))
  band = []
  npix = []
  for entry in xbandfiles:
    print('entry=',entry)
    match = re.search('-[\w]*.wav',entry)
    tag = match.group()[1:-4]
    if match: band.append(tag.upper())
    x = loadtxt(proot+'-'+tag+'.wav')
    npix.append(len(x))
  
  print('proot+.wav=',proot+'.wav')  
  print('xbandfiles=',xbandfiles)
  x = np.loadtxt(proot+'.wav')
  if len(npix) == 0: npix.append(len(x))
    
  m=glob.glob(proot+".mdl")
  e=glob.glob(proot+".err")
  n=glob.glob(proot+".nrd")

  fmp=glob.glob(proot+".fmp.fits")
  scr=glob.glob(proot+".scr.fits")
 
  if len(fmp) > 0:
    ff=fits.open(fmp[0])
    fibermap=ff[1]

  if len(scr) > 0:
    fs=fits.open(scr[0])
    scores=fs[1]
 
  success=[]
  targetid=[]
  target_ra=[]
  target_dec=[]
  ref_id=[]
  ref_cat=[]
  srcfile=[]
  bestgrid=[]
  teff=[]
  logg=[]
  feh=[]
  alphafe=[]
  micro=[]
  param=[]
  covar=[]
  elem=[]
  elem_err=[]
  snr_med=[]
  chisq_tot=[]
  rv_adop=[]
  rv_err=[]

  of=open(o[0],'r')
    
  for line in of:
	  
    cells=line.split()
    #k = int(cells[0])  # the very first line gives the index (1,2...) for the successful grid
    bestgrid.append(o[0])
    id = cells[0]
    cells = cells[1:]

    ncells = len(cells)
    ndim = ncells - 3
    if ndim > 10:
        ndim = int(np.sqrt(4*ndim+1)-1)
    ndim = ndim // 2
    
    assert (ncells > 6), 'Error, the file '+o[0]+' has less than 7 columns, which would correspond to ndim=2'
    par = np.zeros(ndim)
    #cov = np.zeros(ndim*ndim+ndim) 
    cov = np.zeros((5,5))
        
    print('ndim=',ndim)
    print('ncells=',ncells)


    if (ndim == 2):
      #white dwarfs 2 dimensions: id, 2 par, 2err, 0., med_snr, lchi, 2x2 cov
      feh.append(-10.)
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(np.nan)
      micro.append(np.nan)


    elif (ndim == 3):
      #synple grids with Teff, logg and [Fe/H]
      feh.append(float(cells[2]))
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(np.nan)
      micro.append(np.nan)

    elif (ndim == 4):
      #Teff, logg, [Fe/H] and [a/Fe]
      feh.append(float(cells[2]))
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(float(cells[3]))
      micro.append(np.nan)
   
    elif (ndim == 5):
      #Teff, logg, [Fe/H], [a/Fe] and ?
      #see Allende Prieto et al. (2018, A&A)
      feh.append(float(cells[2]))
      teff.append(float(cells[0]))
      logg.append(float(cells[1]))
      alphafe.append(float(cells[3]))
      micro.append(np.nan)


    chisq_tot.append(10.**float(cells[2+2*ndim]))
    snr_med.append(float(cells[1+2*ndim]))
    rv_adop.append(float(cells[0+2*ndim]))
    rv_err.append(np.nan)
    par = np.array(cells[0:ndim],dtype=float)
    cov = np.array(cells[2*ndim:],dtype=float)
    param.append(par)
    covar.append(cov)
    
    targetid.append(np.int64(id))
    srcfile.append(root)

    if (chisq_tot[-1] < 1.5 and snr_med[-1] > 5.): # chi**2<1.5 and S/N>5
      success.append(1) 
    else: success.append(0)



  #add info copied from fibermap
  nspec = len(targetid)
  try:
    #targetid=fibermap.data['targetid']
    target_ra=fibermap.data['target_ra']
    target_dec=fibermap.data['target_dec']
    ref_id=fibermap.data['ref_id']
    ref_cat=fibermap.data['ref_cat']
  except NameError:
    target_ra=np.zeros(nspec)
    target_dec=np.zeros(nspec)
    ref_id=np.zeros(nspec,dtype=int64)
    ref_cat=np.array(["" for x in range(nspec)])

  #primary extension
  hdu0=fits.PrimaryHDU()

  #find out processing date and add it to primary header
  now = datetime.datetime.fromtimestamp(time.time())
  nowstr = now.isoformat() 
  nowstr = nowstr[:nowstr.rfind('.')]
  hdu0.header['DATE'] = nowstr
  #hdu0.header['FCONFIG'] = config

  #find out host machine and add info to header
  try:
    host=os.environ['HOST']
  except:
    host='Unknown'
  hdu0.header['HOST'] = host
  #find out OS name/platform
  osname = os.name 
  platf = platform.system() + ' '+ platform.release()
  hdu0.header['OS'] = osname
  hdu0.header['PLATFORM'] = platf

  #keep track of the number of targets processed and the time it took
  hdu0.header['NSPEC'] = nspec
  #ftiming = get_ferre_timings(proot)
  #hdu0.header['FTIME'] = ftiming
  #stiming = get_slurm_timings(proot)
  #hdu0.header['STIME'] = stiming
  #ncores = get_slurm_cores(proot)
  #hdu0.header['NCORES'] = ncores

  #get versions and enter then in primary header
  ver = get_versions()
  for entry in ver.keys(): hdu0.header[entry] = ver[entry]
  
  hdulist = [hdu0]

  #sptab extension
  cols = {}
  cols['SUCCESS'] = success
  cols['TARGETID'] = targetid
  cols['TARGET_RA'] = target_ra
  cols['TARGET_DEC'] = target_dec
  cols['REF_ID'] = ref_id
  cols['REF_CAT'] = ref_cat
  cols['SRCFILE'] = srcfile
  cols['BESTGRID'] = bestgrid
  cols['TEFF'] = np.array(teff)*units.K
  cols['LOGG'] = np.array(logg)
  cols['FEH'] = np.array(feh)
  cols['ALPHAFE'] = np.array(alphafe) 
  cols['LOG10MICRO'] = np.array(micro)
  cols['PARAM'] = np.vstack ( (teff, logg, feh, alphafe, micro) ).T
  cols['COVAR'] = np.array(covar)  #.reshape(len(success),5,5)
  #cols['ELEM'] = np.array(elem)
  #cols['ELEM_ERR'] = np.array(elem_err)
  cols['CHISQ_TOT'] = np.array(chisq_tot)
  cols['SNR_MED'] = np.array(snr_med)
  cols['RV_ADOP'] = np.array(rv_adop)*units.km/units.s
  cols['RV_ERR'] = np.array(rv_err)*units.km/units.s

  colcomm = {
  'success': 'Bit indicating whether the code has likely produced useful results',
  'TARGETID': 'DESI targetid',
  'TARGET_RA': 'Target Right Ascension (deg) -- details in FIBERMAP',
  'TARGET_DEC': 'Target Declination (deg) -- details in FIBERMAP',
  'REF_ID': 'Astrometric cat refID (Gaia SOURCE_ID)',
  'REF_CAT': 'Astrometry reference catalog',
  'SRCFILE': 'DESI data file',
  'BESTGRID': 'Model grid that produced the best fit',
  'TEFF': 'Effective temperature (K)',
  'LOGG': 'Surface gravity (g in cm/s**2)',
  'FEH': 'Metallicity [Fe/H] = log10(N(Fe)/N(H)) - log10(N(Fe)/N(H))sun' ,
  'ALPHAFE': 'Alpha-to-iron ratio [alpha/Fe]',
  'LOG10MICRO': 'Log10 of Microturbulence (km/s)',
  'PARAM': 'Array of atmospheric parameters ([Fe/H], [a/Fe], log10micro, Teff,logg)',
  'COVAR': 'Covariance matrix for ([Fe/H], [a/Fe], log10micro, Teff,logg)',
  #'ELEM': 'Elemental abundance ratios to hydrogen [elem/H]',
  #'ELEM_ERR': 'Uncertainties in the elemental abundance ratios',
  'CHISQ_TOT': 'Total chi**2',
  'SNR_MED': 'Median signal-to-ratio',
  'RV_ADOP': 'Adopted Radial Velocity (km/s)',
  'RV_ERR': 'Uncertainty in the adopted Radial Velocity (km/s)'
  }      

  
  table = tbl.Table(cols)
  hdu=fits.BinTableHDU(table,name = 'SPTAB')
  #hdu.header['EXTNAME']= ('SPTAB', 'Stellar Parameter Table')
  k = 0
  for entry in colcomm.keys():
    print(entry) 
    hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
    k+=1
  hdulist.append(hdu)

  #fibermap extension
  if len(fmp) > 0:
    hdu=fits.BinTableHDU.from_columns(fibermap, name='FIBERMAP')
    hdulist.append(hdu)
    ff.close()

  #scores extension
  if len(scr) > 0:
    hdu=fits.BinTableHDU.from_columns(scores, name='SCORES')
    hdulist.append(hdu)
    fs.close()

  #aux extension
  #p = ['[Fe/H]','[a/Fe]','log10micro','Teff','logg']
  #if 'elem' in conf: e = conf['elem']
  #cols = {}
  #colcomm = {}
  #cols['p'] = [p]
  #colcomm['p'] = 'PARAM tags'
  #if 'elem' in conf:
  #  cols['e'] = [e]
  #  colcomm['e']= 'ELEM tags'
  
  #table = tbl.Table(cols)
  #hdu=fits.BinTableHDU(table,name = 'AUX')

  #k = 0
  #for entry in colcomm.keys():
  #  print(entry) 
  #  hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
  #  k+=1
  #hdulist.append(hdu)


  hdul=fits.HDUList(hdulist)
  hdul.writeto(os.path.join(path,'sptab_'+root+'.fits'))
  
  #now spmod
  hdulist = [hdu0]
 

  edata=np.loadtxt(e[0])
  mdata=np.loadtxt(m[0])
  odata=np.loadtxt(n[0])
  
  i = 0
  j1 = 0

  if len(xbandfiles) > 0:
    for entry in band:
      j2 = j1 + npix[i] 
      print(entry,i,npix[i],j1,j2)
      #colx = fits.Column(name='wavelength',format='e8', array=array(x[j1:j2]))
      #coldefs = fits.ColDefs([colx])
      #hdu = fits.BinTableHDU.from_columns(coldefs)
      hdu = fits.ImageHDU(name=entry+'_WAVELENGTH', data=x[j1:j2])
      hdulist.append(hdu)
    
      cols = {}
      colcomm = {}
      if odata.ndim == 2: tdata = odata[:,j1:j2]
      else: tdata = odata[j1:j2][None,:]
      cols['obs'] = tdata
      colcomm['obs'] = 'Observed spectra as fit'
      if edata.ndim == 2: tdata = edata[:,j1:j2]
      else: tdata = edata[j1:j2][None,:]
      cols['err'] = tdata
      colcomm['err'] = 'Error in spectra as fit'
      if mdata.ndim == 2: tdata = mdata[:,j1:j2]
      else: tdata = mdata[j1:j2][None,:]
      cols['fit'] = tdata
      colcomm['fit'] = 'Best-fitting model'

      table = tbl.Table(cols)
      hdu=fits.BinTableHDU(table,name = entry+'_MODEL')
      k = 0
      for entry in colcomm.keys():
        print(entry) 
        hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
        k+=1
      hdulist.append(hdu)
      i += 1
      j1 = j2
      
  else:
	  
    print('single wavelength entry!')
    #colx = fits.Column(name='wavelength',format='e8', array=array(x[j1:j2]))
    #coldefs = fits.ColDefs([colx])
    #hdu = fits.BinTableHDU.from_columns(coldefs)
    hdu = fits.ImageHDU(name='WAVELENGTH', data=x)
    hdulist.append(hdu)
      
    cols = {}
    colcomm = {}
    if odata.ndim == 2: tdata = odata[:,:]
    else: tdata = odata[:][None,:]
    cols['obs'] = tdata
    colcomm['obs'] = 'Observed spectra as fit'
    if edata.ndim == 2: tdata = edata[:,:]
    else: tdata = edata[:][None,:]
    cols['err'] = tdata
    colcomm['err'] = 'Error in spectra as fit'
    if mdata.ndim == 2: tdata = mdata[:,:]
    else: tdata = mdata[:][None,:]
    cols['fit'] = tdata
    colcomm['fit'] = 'Best-fitting model'      
      

    table = tbl.Table(cols)
    hdu=fits.BinTableHDU(table,name = 'MODEL')
    k = 0
    for entry in colcomm.keys():
      print(entry) 
      hdu.header['TCOMM'+"{:03d}".format(k+1)] = colcomm[entry]
      k+=1
    hdulist.append(hdu)
	  

  if len(fmp) > 0:
    ff=fits.open(fmp[0])
    fibermap=ff[1]
    hdu=fits.BinTableHDU.from_columns(fibermap, name='FIBERMAP')
    #hdu.header['EXTNAME']='FIBERMAP'
    hdulist.append(hdu)

  if len(scr) > 0:
    ff=fits.open(scr[0])
    scores=ff[1]
    hdu=fits.BinTableHDU.from_columns(scores, name='SCORES')
    hdulist.append(hdu)

  hdul=fits.HDUList(hdulist)
  hdul.writeto(os.path.join(path,'spmod_'+root+'.fits')) 
  
  return None

#get dependencies versions, shamelessly copied from rvspec (Koposov's code)
def get_dep_versions():
    """
    Get Packages versions
    """
    
    import importlib
  
    packages = [
        'numpy', 'astropy', 'matplotlib', 'scipy',
        'yaml'
    ]
    # Ideally you need to check that the list here matches the requirements.txt
    ret = {}
    for curp in packages:
        ret[curp[:8]] = importlib.import_module(curp).__version__
    ret['python'] = str.split(sys.version, ' ')[0]
    return ret



#find out versions 
def get_versions():

  ver = get_dep_versions()
  ver['synple'] = 0
  log1file = glob.glob("*.log_01")
  fversion = 'unknown'
  if len(log1file) < 1:
    print("Warning: cannot find any *.log_01 file in the working directory")
  else:
    l1 = open(log1file[0],'r')
    #while 1:
    #  line = l1.readline()
    for line in l1:
      if 'f e r r e' in line:
        entries = line.split()
        fversion = entries[-1][1:]
        break
    l1.close()
    ver['ferre'] = fversion

  return(ver)

    
def fparams(root,synthfile=None,figure=None):
    """"Evaluates the agreement between input/output parameters in 
    FERRE (root.ipf and root.opf) files
    
    Parameters
    ----------
    root: str
      root for the input/output parameter FERRE files
    synthfile: str
      associated FERRE/BAS synthfile to look for labels and other info
    figure: str
      name for the output file with a summary figure (usually a png file)
     (default is None, and the plot is only shown on the screen)
      
    Returns
    -------
    result: numpy array of floats
      16-50-84 percentiles for all the parameters
    
    """
    
    if synthfile is None:
      synthfile = 'f_'+root+'.dat'
      if not os.path.exists(synthfile): 
        synthfile = 'n_'+root+'.dat'
        if not os.path.exists(synthfile): 
          synthfile = 'l_'+root+'.dat'
          if not os.path.exists(synthfile): 
            synthfile = 'f_'+root+'.pickle'
            if not os.path.exists(synthfile): 
              synthfile = 'n_'+root+'.pickle'
              if not os.path.exists(synthfile):
                synthfile = 'l_'+root+'.pickle'
                assert os.path.exists(synthfile),'cannot find synthfile:'+synthfile
                
    h = head_synth(synthfile)
    if type(h) is list: h = h[1] 
    ndim = int(h['N_OF_DIM'])
    v = np.loadtxt(root + '.ipf',usecols=np.arange(ndim)+1)
    o = np.loadtxt(root + '.opf',usecols=np.arange(ndim)+1)
    result = np.percentile(o-v,[15.85,50.,84.15],axis=0)


    plt.clf()
    for i in range(ndim):
        plt.subplot(2,ndim,i+1)
        plt.plot(v[:,i],o[:,i],'.')
        if 'LABEL('+str(i+1)+')' in h: plt.title(h['LABEL('+str(i+1)+')'])
        if i == 0: plt.ylabel('output')
    for i in range(ndim):
        plt.subplot(2,ndim,i+1+ndim)
        hi = plt.hist(o[:,i]-v[:,i],bins=20)
        if i == 0: plt.xlabel('difference')
        plt.text(hi[1][0],max(hi[0])*0.6, "p16" % result[0,i] )
        plt.text(hi[1][0],max(hi[0])*0.4, "p50" % result[1,i] )
        plt.text(hi[1][0],max(hi[0])*0.2, "p84" % result[2,i] )
        plt.text(result[0,i],max(hi[0])*0.6, "%6.2f" % result[0,i] )
        plt.text(result[1,i],max(hi[0])*0.4, "%6.2f" % result[1,i] )
        plt.text(result[2,i],max(hi[0])*0.2, "%6.2f" % result[2,i] )
 
    manager = plt.get_current_fig_manager() 
    #manager.resize(*manager.window.maxsize())

    plt.ion()
    plt.show()

    if figure is not None:
        plt.savefig(figure)
    
    return(result)


def desida(path_to_data='healpix',path_to_output='sp_output',
           synthfile=None, seconds_per_target=2.,star=True):

  """ Prepare a DESI data for parallel processing
  """

  
  python_path1=os.environ['HOME']+"/synple"
  pwd=os.path.abspath(os.curdir)

  try:
    host=os.environ['HOST']
  except:
    host='Unknown'
  now=time.strftime("%c")

  if synthfile is None:
    synthfile1 = 'None'
  else:
    synthfile1 = "'"+str(synthfile)+"'"

 
  infiles = list(glob.iglob(os.path.join(path_to_data,'**',
            'coadd*fits'), recursive=True)) 

  folders = []
  for entry in infiles:
    parts = entry.split('/')
    infile = parts[-1]
    folder = ('/'.join(parts[-5:-1]))
    tpath = os.path.join(path_to_output,folder)

    os.makedirs(tpath,exist_ok=True) 

    minutes = 4000*seconds_per_target/60.
    root = infile[:-5]

    sfile = os.path.join(tpath,root+'.job') 
    outfile = os.path.join(tpath,root)

    print('infile=',entry)
    print('outfile=',outfile)
    s = open(sfile,'w')
    s.write("#!/bin/bash \n")
    s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
    s.write("#This script was written by synple.py on "+now+" \n")
    s.write("#SBATCH --time="+str(int(minutes)+1)+"\n") #minutes
    s.write("#SBATCH --ntasks=1" + "\n")
    s.write("#SBATCH --nodes=1" + "\n")
    if (host == 'login1'): #lapalma
      nthreads = 4
      s.write("#SBATCH  -J "+str(root)+" \n")
      s.write("#SBATCH  -o "+str(root)+"_%j.out"+" \n")
      s.write("#SBATCH  -e "+str(root)+"_%j.err"+" \n")
      s.write("#SBATCH --cpus-per-task="+str(16)+"\n")
      s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
      s.write("module load python"+"\n")
    else: # perlmutter
      nthreads = 4
      s.write("#SBATCH --qos=regular" + "\n")
      s.write("#SBATCH --constraint=cpu" + "\n")
      s.write("#SBATCH --account=desi \n")
      s.write("#SBATCH --cpus-per-task="+str(128*2)+"\n")
      s.write("#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# \n")
      s.write("module load python"+"\n")

    s.write("cd "+pwd+"\n\n")

    command="python3 -c \"import sys; " + \
     " sys.path.insert(0, '"+python_path1 + "'); " + \
     " nthreads = " + str(nthreads) + ";" + \
     " import os; os.environ['OMP_NUM_THREADS'] = str(nthreads) ; "+ \
     " os.environ['OPENBLAS_NUM_THREADS'] = str(nthreads) ;"  + \
     " os.environ['MKL_NUM_THREADS'] = str(nthreads) ; " + \
     " from synple import bas, wtabmodfits; " + \
     " bas(\'" + entry + "\'," + \
     " outfile=\'" + outfile + "\'," + \
     " synthfile=" + str(synthfile1) + ", star= " + str(star) + "); " + \
     " wtabmodfits(\'" + root + "'" + ", path= '" + tpath + "\'" + \
     ")\"" + "\n"

    s.write(command)
    s.close()
    os.chmod(sfile,0o755)
 


  return()


def desimask(desi_target):

  """Returns targetting classes from DESI_TARGET bits in the FIBERMAP
  extension of DESI data files
  """

  mask = [('LRG', 0, 1),
  ('ELG', 1, 2),
  ('QSO', 2, 4),
  ('LRG_1PASS', 3, 8),
  ('LRG_2PASS', 4, 16),
  ('LRG_NORTH', 8, 256),
  ('ELG_NORTH', 9, 512),
  ('QSO_NORTH', 10, 1024),
  ('LRG_SOUTH', 16, 65536),
  ('ELG_SOUTH', 17, 131072),
  ('QSO_SOUTH', 18, 262144),
  ('LRG_1PASS_NORTH', 24, 16777216),
  ('LRG_2PASS_NORTH', 25, 33554432),
  ('LRG_1PASS_SOUTH', 28, 268435456),
  ('LRG_2PASS_SOUTH', 29, 536870912),
  ('SKY', 32, 4294967296),
  ('STD_FAINT', 33, 8589934592),
  ('STD_WD', 34, 17179869184),
  ('STD_BRIGHT', 35, 34359738368),
  ('BAD_SKY', 36, 68719476736),
  ('SUPP_SKY', 37, 137438953472),
  ('NO_TARGET', 49, 562949953421312),
  ('BRIGHT_OBJECT', 50, 1125899906842624),
  ('IN_BRIGHT_OBJECT', 51, 2251799813685248),
  ('NEAR_BRIGHT_OBJECT', 52, 4503599627370496),
  ('BGS_ANY', 60, 1152921504606846976),
  ('MWS_ANY', 61, 2305843009213693952),
  ('SCND_ANY', 62, 4611686018427387904)]

  bits = []
  target = desi_target
  for entry in reversed(mask):
    if target - entry[2] >= 0: 
      target = target - entry[2]
      bits.append(entry[0])


  return(bits)


def mwsmask(mws_target):

  """Returns targetting classes from MWS_TARGET bits in the FIBERMAP
  extension of DESI data files
  """
  mask = [('MWS_BROAD', 0, 1),
 ('MWS_WD', 1, 2),
 ('MWS_NEARBY', 2, 4),
 ('MWS_BROAD_NORTH', 4, 16),
 ('MWS_BROAD_SOUTH', 5, 32),
 ('MWS_BHB', 6, 64),
 ('MWS_MAIN_BLUE', 8, 256),
 ('MWS_MAIN_BLUE_NORTH', 9, 512),
 ('MWS_MAIN_BLUE_SOUTH', 10, 1024),
 ('MWS_MAIN_RED', 11, 2048),
 ('MWS_MAIN_RED_NORTH', 12, 4096),
 ('MWS_MAIN_RED_SOUTH', 13, 8192),
 ('MWS_FAINT_BLUE', 14, 16384),
 ('MWS_FAINT_BLUE_NORTH', 15, 32768),
 ('MWS_FAINT_BLUE_SOUTH', 16, 65536),
 ('MWS_FAINT_RED', 17, 131072),
 ('MWS_FAINT_RED_NORTH', 18, 262144),
 ('MWS_FAINT_RED_SOUTH', 19, 524288),
 ('GAIA_STD_FAINT', 33, 8589934592),
 ('GAIA_STD_WD', 34, 17179869184),
 ('GAIA_STD_BRIGHT', 35, 34359738368),
 ('BACKUP_DIB', 57, 144115188075855872),
 ('BACKUP_GIANT_LOP', 58, 288230376151711744),
 ('BACKUP_GIANT', 59, 576460752303423488),
 ('BACKUP_BRIGHT', 60, 1152921504606846976),
 ('BACKUP_FAINT', 61, 2305843009213693952),
 ('BACKUP_VERY_FAINT', 62, 4611686018427387904)]

  bits = []
  target = mws_target
  for entry in reversed(mask):
    if target - entry[2] >= 0: 
      target = target - entry[2]
      bits.append(entry[0])


  return(bits)


if __name__ == "__main__":

  npar = len(sys.argv)
  assert (npar >= 4), 'Synple requires at least 3 input parameters (modelfile wstart wend)'
  assert (npar <= 7), 'Synple requires at maximum 6 input parameters (modelfile wstart wend vmicro vrot fwhm)'
  vmicro = None
  vrot = 0.0
  fwhm = 0.0
  modelfile = sys.argv[1]
  wstart = float(sys.argv[2])
  wend = float(sys.argv[3])
  if (npar > 4): 
    vmicro = float(sys.argv[4])
    if (npar > 5):
      fwhm = float(sys.argv[5])
      if (npar > 6):
        vrot = float(sys.argv[6])

  #symbol, mass, sol = elements()
  s = syn(modelfile, (wstart,wend), save=True, vmicro=vmicro, vrot=vrot, fwhm=fwhm)


