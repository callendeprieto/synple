#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Python wrapper for synspec 

Calculation of synthetic spectra of stars and convolution with a rotational/Gaussian kernel.
Makes the use of synspec simpler, and retains the main functionalities (when used from
python). The command line interface for the sheread_ll is even simpler but fairly limited. 

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
import sys
import stat
import string
import random
import subprocess
import numpy as np
import glob
import time
import copy
import gzip
from math import ceil
from scipy import interpolate
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from itertools import product

#configuration
#synpledir = /home/callende/synple
synpledir = os.path.dirname(os.path.realpath(__file__))

#relative paths
modeldir = synpledir + "/models"
modelatomdir = synpledir + "/data"
linelistdir = synpledir + "/linelists"
linelist0 = ['gfATO.19','gfMOLsun.20','gfTiO.20','H2O-8.20']
bindir = synpledir + "/bin"
synspec = bindir + "/synspec54"
rotin = bindir + "/rotin"

#internal synspec data files
isdf = ['CIA_H2H2.dat',  'CIA_H2H.dat', 'CIA_H2He.dat', 'CIA_HHe.dat', \
        'irwin_bc.dat',  'tremblay.dat', \
        'tsuji.molec_bc2']

#other stuff
clight = 299792.458
epsilon = 0.6 #clv coeff.
bolk = 1.38054e-16  # erg/ K
zero = " 0 "
one =  " 1 "
two =  " 2 "



def syn(modelfile, wrange, dw=None, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0,  lineid=False, tag=False,  \
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
  vmicro: float, optional
      microturbulence (km/s) 
      (default is None, which is overriden by the value from the model atmosphere)
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

  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2
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
  print ('teff,logg,vmicro=',teff,logg,vmicro)
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
        

  #link the data folder (synspec data + default tlusty model atoms)
  if dd != 'data': os.symlink(modelatomdir,'./data')          #data directory  

  write5(teff,logg,abu,atom,inlte=inlte,atommode=atommode,atominfo=atominfo)   #abundance/opacity file  
  write8(teff,logg,nd,atmos,atmostype)                    #model atmosphere

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
          linelist=linelist,atmostype=atmostype)
  #synspec control file
  writetas('tas',nd,linelist,nonstd=nonstd)               #non-std param. file
  create_links(linelist)                                  #auxiliary data

  if compute == False:

    wave = None
    flux = None  
    cont = None

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
    s = wave, flux, cont, [la,li,lo]
  else:
    s = wave, flux, cont

  if tag: tags(s)

  return(s)


def mpsyn(modelfile, wrange, dw=None, strength=1e-4, vmicro=None, abu=None, \
    linelist=linelist0, atom='ap18', vrot=0.0, fwhm=0.0, vmacro=0.0, \
    steprot=0.0, stepfwhm=0.0, lineid=False, tag=False,  \
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
  vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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

  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2


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
      steprot, stepfwhm,  lineid, tag, clean, False, None, lte, \
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
    steprot=0.0, stepfwhm=0.0,  lineid=False, tag=False, \
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
  vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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
    atom,vrot,fwhm,vmacro,steprot,stepfwhm,lineid,tag,clean,save,synfile,compute = cons

    s = syn(modelfile, wrange, dw, strength, vmicro, abu, \
              linelist, atom, vrot, fwhm, vmacro, \
              steprot, stepfwhm,  lineid, tag, clean, save, synfile, \
              lte, compute, tmpdir)

    return(s)


  #basic checks on the line list and model atmosphere
  linelist, modelfile = checksynspec(linelist,modelfile)

  #read model atmosphere
  atmostype, teff, logg, vmicro2, abu2, nd, atmos = read_model(modelfile)

  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2

  if nthreads == 0: 
    nthreads = int(psutil.cpu_count(logical=False) - 1)

  if tag: lineid = True

  print('nthreads=',nthreads)

  tmpdir = ''.join(random.choices(string.ascii_lowercase + string.digits, k = 16))

  ray.init(num_cpus=nthreads)

  rest = [ modelfile,dw,strength,vmicro,abu,linelist, \
    atom,vrot,fwhm,vmacro,steprot,stepfwhm,lineid,tag,clean,False,None,compute ]

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
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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
  vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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

  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2


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
        s.write("#SBATCH  -t 04:00:00"+" \n") #hh:mm:ss
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
  vmicro: float, optional
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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

  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2

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
    vmicro=None, vrot=0.0, fwhm=0.0, vmacro=0.0, nfe=0.0, \
    linelist=linelist0, atom='ap18', \
    steprot=0.0, stepfwhm=0.0,  clean=True, save=None, lte=True):

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
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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

  Returns
  -------
  wave: numpy array of floats (1D)
      wavelengths (angstroms)
  flux: numpy array of floats (2D -- as many rows as models input)
      flux (H_lambda in ergs/s/cm2/A)
  cont: numpy array of floats (2D -- as many rows as models input)
      continuum flux (same units as flux)

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
  try: 
    nnfe = len(nfe)
    nnfes = nfe
  except TypeError:
    nnfe = 1
    nnfes = [ nfe ] 


  idir = 0
  dirfile = open('dirtree.txt','w')
  for entry in modelfiles:
    for vmicro1 in vmicros:
      for nfe1 in nnfes:

        idir = idir + 1
        dir = ( "hyd%07d" % (idir) )
        dirfile.write(str(idir)+' folder='+dir+' model='+entry+' vmicro='+str(vmicro1)+' [N/Fe]='+str(nfe1)+' \n')
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
          if (abs(nfe1) > 1e-7):
            if (abu1 == None):
              linelist, entry = checksynspec(linelist,entry)
              atmostype, teff, logg, vmicro2, abu1, nd, atmos = read_model(entry)
            abu1[6] = abu1[6] * 10.**nfe1

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


  return(None,None,None)



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




def collect_marcs(modeldir=modeldir, tteff=None, tlogg=None, \
                  tfeh=(1,0.0,0.0), tafe=(1,0.0,0.0), \
                  tcfe=(1,0.0,0.0), tnfe=(1,0.0,0.0), \
                  tofe=(1,0.0,0.0), trfe=(1,0.0,0.0), tsfe=(1,0.0,0.0), \
                  ignore_missing_models=False, ext='mod'):

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

                    file = glob.glob(os.path.join(modeldir,filename))
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
  tcfe=(1,0.0,0.0), ignore_missing_models=False, ext='mod'):

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
                    print(file)
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
    ignore_missing_models=False):

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
           vmicro=1.0, nfe=0.0, vrot=0.0, fwhm=0.0, vmacro=0.0, 
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
           vmicro=None, vrot=None, fwhm=None, nfe=None,
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
                              if fwhm is not None:
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
              vmicro=1.0, nfe=0.0, vrot=0.0, fwhm=0.0, vmacro=0.0):	  
  
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
  if np.abs(np.max(nfe)) > 1e-7 and len(nfe) > 1:
    ndim = ndim + 1
    n_p.append(len(nfe))
    labels.append('[N/Fe]')
    llimits.append(nfe[0])
    steps.append(nfe[1]-nfe[0])
    dnfe=np.diff(nfe)
    assert np.max(dnfe) - np.min(dnfe) < 1.e-7, '[N/Fe] values are not linearly spaced!'
  if np.abs(np.max(vrot)) > 1e-7 and len(vrot) > 1:
    ndim = ndim + 1
    n_p.append(len(vrot))
    labels.append('vrot')
    llimits.append(vrot[0])
    steps.append(vrot[1]-vrot[0])
    dvrot=np.diff(vrot)
    if np.max(dvrot) - np.min(dvrot) > 1.e-7:
      vrot = np.log10(vrot)
      dvrot=np.diff(vrot)
      assert np.max(dvrot) - np.min(dvrot) < 1.e-7, 'Vrot values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10vrot'
      llimits[-1]=vrot[0]
      steps[-1]=vrot[1]-vrot[0]
  if np.abs(np.max(fwhm)) > 1e-7 and len(fwhm) > 1:
    ndim = ndim + 1
    n_p.append(len(fwhm))
    labels.append('FWHM')
    llimits.append(fwhm[0])
    steps.append(fwhm[1]-fwhm[0])
    dfwhm=np.diff(fwhm)
    if np.max(dfwhm) - np.min(dfwhm) > 1.e-7:
      fwhm = np.log10(fwhm)
      dfwhm=np.diff(fwhm)
      assert np.max(dfwhm) - np.min(dfwhm) < 1.e-7, 'FWHM values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10FWHM'
      llimits[-1]=fwhm[0]
      steps[-1]=fwhm[1]-fwhm[0]
  if np.abs(np.min(vmacro)) > 1e-7 and len(vmacro) > 1:
    ndim = ndim + 1
    n_p.append(len(vmacro))
    labels.append('VMACRO')
    llimits.append(vmacro[0])
    steps.append(vmacro[1]-vmacro[0])
    dvmacro=np.diff(vmacro)
    if np.max(dvmacro) - np.min(dvmacro) > 1.e-7:
      vmacro = np.log10(vmacro)
      dvmacro=np.diff(vmacro)
      assert np.max(dvmacro) - np.min(dvmacro) < 1.e-7, 'vmacro values are neither linearly spaced or linearly spaced in log!'
      labels[-1]='log10vmacro'
      llimits[-1]=vmacro[0]
      steps[-1]=vmacro[1]-vmacro[0]


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
           vmicro=1.0, vrot=0.0, fwhm=0.0, vmacro=0.0, 
           wrange=None, dw=None, logw=0, ignore_missing_models=False,**elements):



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
  feh: tuple
    Activate to track this parameter
  vmicro: float, optional, can be an iterable
      microturbulence (km/s) 
      (default is taken from the model atmosphere)
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
    
  **elements: booleans
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
  for entry in elements.keys():
    pars.append(entry)
  if nvrot > 1: pars.append('vrot')
  if nfwhm > 1: pars.append('fwhm')
  if nvmacro > 1: pars.append('vmacro')


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

                    teff,logg,vmicro2,abu = read_madaf(madaffile,startdir=entry)
                    imode, iprin, inmod, inlte, hydprf, wrange, cutoff,  \
                         strength, dw, molls, vmicro1 = read55(os.path.join(entry,'fort.55'))
                    feh = np.log10(abu[25])+12-7.50
	                                
                    print(teff,logg,feh,vmicro1)
                    
                    ntot = ntot + 1

                    iconv = 1
                    outconv = ("%07dfort.7" % (iconv) )
                    file = os.path.join(dir,outconv)
                     
                    if break_out == False and os.path.isfile(file):
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


                    teff,logg,vmicro2,abu = read_madaf(madaffile,startdir=entry)
                    imode, iprin, inmod, inlte, hydprf, wrange, cutoff, \
                         strength, dw, molls, vmicro1 = read55(os.path.join(entry,'fort.55'))
                    feh = np.log10(abu[25])+12-7.50
	                  
                    print(teff,logg,feh,vmicro1)

                    pars = []
                    if teff: pars.append(teff)
                    if logg: pars.append(logg)
                    if feh: pars.append(feh)
                    if nvmicro > 1: pars.append(vmicro1)
                    for el in elements.keys():
                      pars.append(el)

                    iconv = 0
                    for vrot1 in vrots:
                      for fwhm1 in fwhms:
                        for vmacro1 in vmacros:
								
                          iconv = iconv + 1
                          outconv = ("%07dfort.7" % (iconv) )
                          file = os.path.join(dir,outconv)
 
                          if nvrot > 1: pars.append(vrot1)
                          if nfwhm > 1: pars.append(fwhm1)
                          if nvmacro > 1: pars.append(vmacro1)

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
                          **kargs):
							  
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
                          **kargs):
							  
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

       print(comm)
       s.write(comm+'\n')
       st = os.stat(sfile)
       os.chmod(sfile, st.st_mode | stat.S_IEXEC)

    s.close()
    
    return()

  
#extract the header of a synthfile
def head_synth(synthfile):
    meta=0
    multi=0
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
            if (meta > 1): multi_header.append(header)
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
    xx=[]
    j=0
    for header in multi_header:
      tmp=header['WAVE'].split()
      npix=int(header['NPIX'])
      step=float(tmp[1])
      x0=float(tmp[0])
      x=np.arange(npix)*step+x0
      if header['LOGW']:
        if int(header['LOGW']) == 1: x=10.**x
        if int(header['LOGW']) == 2: x=np.exp(x)   
      j=j+1
      xx.append(x)

    if len(xx)>1: x=xx[:]

    return x

#read a synthfile
def read_synth(synthfile):
    """
    Reads a FERRE spectral grid from disk
    """
	
    meta=0
    multi=0
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
            if (meta > 1): multi_header.append(header)
            header={}
            line=file.readline()
            nlines+=1
        else:
          k=part[0].strip()
          v=part[1].strip()
          header[k]=v
          if k == 'MULTI': 
            multi=int(v)
            multi_header=[]
    if (multi > 1): header=multi_header
    file.close()

    if np.ndim(header) > 0: 
      snp=header[0]['N_P'] 
    else: 
      snp=header['N_P']

    n_p = tuple(np.array(snp.split(),dtype=int)) + (-1,)
    data=np.loadtxt(synthfile, skiprows=nlines, dtype=float)
    data = np.reshape( data, n_p)
    

    return header,data
    
def write_synth(synthfile,d,hdr=None):
    """
    Writes a FERRE spectral grid to disk
    """
	
    ndim = d.ndim-1
    n_p = d.shape[:-1]
    npix = d.shape[-1]
    
    pwd=os.path.abspath(os.curdir)
    nowtime=time.ctime(time.time())
    osinfo=os.uname()
   

    if hdr is None:
        #minimal header        
        hdr = {}
        hdr['DATE'] = "'"+nowtime+"'"
        hdr['N_OF_DIM'] = str(ndim)
        hdr['N_P'] = '  '.join(map(str,n_p))
        for i in range(ndim): hdr['LABEL('+str(i+1)+")"] = "'"+"unknown"+"'"
        hdr['LLIMITS'] = '  '.join(map(str,np.zeros(ndim)))
        hdr['STEPS'] = '  '.join(map(str,np.ones(ndim)))
        hdr['COMMENTS1'] = "'created by write_synth, without axis information'"
        hdr['COMMENTS2'] = "'"+osinfo[0]+' '+osinfo[2]+'.'+osinfo[4]+' running on '+osinfo[1]+"'"
        hdr['COMMENTS3'] = "'pwd is "+pwd+"'"

    fout = open(synthfile,'w')
    fout.write(' &SYNTH\n')
    for entry in hdr: fout.write(' '+entry + ' = ' + hdr[entry] + '\n')
    fout.write(' /\n')

    #now the data	
    if ndim > 1:
        dd = np.reshape( d, (np.product(n_p), npix) )
    for entry in range(np.product(n_p)):
        dd[entry,:].tofile(fout,sep=" ",format="%0.4e")
        fout.write("\n")

    return(None)			
			
def fill_synth(d,function='cubic'):
    """
    Completes data rows with zeros interpolating using Rbf from 
    non-zero rows 
    """
    
    from scipy.interpolate import Rbf
    
    ndim = d.ndim-1
    n_p = d.shape[:-1]
    nfreq = d.shape[-1]	
    
    dd = np.reshape(d.copy(), (np.product(n_p),nfreq) )
    dd2 = np.sum(dd, dd.ndim-1)
    wi = np.where(dd2 + 1e-31 > 1e-30)[0]
    wo = np.where(dd2 + 1e-31 < 1e-30)[0]

    print(wi)
    print(wo)

    print('ndim=',ndim)
    print('n_p=',n_p)

    #get loop indices for entries in grid
    iarr = getaa(n_p)
    
    for i in np.arange(nfreq):
      print('freq i=',i)
      print('coeff. calculation ...')
      
      rbfi  = Rbf(*np.transpose(iarr[wi]), dd [ wi, i ], 
                  function=function )
      print('interpolation...')
      dd [ wo , i ] = rbfi (*np.transpose(iarr[wo]))
    
    d2 = np.reshape(dd, tuple(n_p)+(nfreq,) ) 
    
    return(d2)
    

def getaa(n_p):
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
      ll.append(np.arange(n_p[i]))
    aa = np.array(list(product(*ll)))
  
    return(aa)
  
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
  plt.xlim([min(x),max(x)])
  plt.ylim([min(y)*0.1,max(y)*1.1])
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
  linelist=linelist0, atmostype='kurucz'):


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
    inlist = 10
    if file[-3:] == '.11' : inlist = 11
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
  if atom == 'hhm' : 
    zex = [1]  #atomic numbers of elements included explicitly (contributing cont. opacity)
  elif atom == 'yo19':
    zex = [1,11,12,19,20]
  elif atom == 'test': 
    zex = [1,26]
  elif atom == 'ap18': 
    zex = [1,2,6,7,8,11,12,13,14,20,26]

  for i in zex: ex[i-1] = 2

  #tlusty models provide atommode and atominfo that override the defaults for atom and ex
  if atommode is not None: ex = np.array(atommode)
  if atominfo is not None: atom = 'own'

  for i in range(natom):
    f.write(' %2d %e %i %s\n' %  (ex[i], abu[i], 0, '  ! ' +symbol[i]) )

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
    f.write("   8    0     54     0      0    0    ' O 1' 'data/o1.t'  \n")
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


def write8(teff, logg, nd, atmos, atmostype, ofile='fort.8'):

  """Writes the model atmosphere for synspec

     MARCS models can be passed in 'Tlusty' (default, after read with 
           read_marcs_models2) or 'Kurucz' format
     Phoenix and Kurucz models are passed to synspec formatted as 'Kurucz'

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
    else:
      pass

  else:

    if atmostype == 'marcs':
      f.write(" "+str(nd)+" "+str(-4)+"\n")
      for i in range(nd):
        f.write(' %e ' % atmos['dm'][i])
      f.write("\n")
      for i in range(nd):
        f.write( '%f %e %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i], atmos['rho'][i]/atmos['mmw'][i]/1.67333e-24 + atmos['ne'][i] ) )

    else:
      f.write( 'TEFF %7.0f  GRAVITY %7.5f  LTE \n' % (teff, logg) )
      for i in range(21): f.write('\n')
      f.write( 'READ DECK6%3i RHOX,T,P,XNE \n' % nd )
      for i in range(nd): 
        f.write( '%e %f %e %e \n' % (atmos['dm'][i], atmos['t'][i], atmos['p'][i], atmos['ne'][i]) )
      
  f.close()

  return()
  

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

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()
    dm.append( float(entries[0]))
    t.append(  float(entries[1]))
    p.append(  float(entries[2]))
    ne.append( float(entries[3]))

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 
  atmos['dm'] = dm
  atmos['t'] = t
  atmos['p'] = p
  atmos['ne'] = ne

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
      atmos['pop'] = atm [:,4:]
    else:
      atmos['dep'] = atm [:,4:]

  return (teff,logg,vmicro,abu,nd,atmos)

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
  abu : list
      abundances, number densities of nuclei relative to hydrogen N(X)/N(H)
      for elements Z=1,99 (H to Es)

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
      entries = entry.replace('\n','').split(',')
      for piece in entries:
        sides = piece.split('=')
        nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')

    print('Tlusty nonstd params=',nonstd)

  #the micro might be encoded as VTB in the nonstdfile!!
  #this is a temporary patch, but need to parse that file
  vmicro = 1.0
  if 'VTB' in nonstd: vmicro = float(nonstd['VTB'])

  line = f.readline()
  line = f.readline()
  entries = line.split()
  natoms = int(entries[0])
  
  abu = []
  for i in range(natoms):
    line = f.readline()
    entries = line.split()
    abu.append( float(entries[1]) )

  if i < 98: 
    for j in range(98-i):
      abu.append(1e-111)
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
      entries = entry.replace('\n','').split(',')
      for piece in entries:
        sides = piece.split('=')
        nonstd[sides[0].replace(' ','')]= sides[1].replace(' ','')


  line = f.readline()
  line = f.readline()
  entries = line.split()
  natoms = int(entries[0])
  
  atommode = []
  for i in range(natoms):
    line = f.readline()
    entries = line.split()
    atommode.append(int(entries[0]))

  atominfo = []
  #keep reading until you find 'dat' to identify data directory 
  line = f.readline()
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
      mf = os.path.join(startdir,datadir)
      if os.path.exists(mf): 
        datadir = mf
      else:
        mf = os.path.join(synpledir,datadir)
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


def gsynth(synthfile,fwhm=0.0,units='km/s',outsynthfile=None,ppr=5,wrange=None,freeze=None):

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
            line = " RESOLUTION = "+str(mean(wrange)/np.sqrt(mean(wrange)**2/resolution**2 + fwhm**2))+"\n"
    if line[1] != "/": fout.write(line)

  try: resolution
  except NameError: 
        if units == 'km/s': 
             line = " RESOLUTION = "+str(clight/fwhm)+"\n"
        else:
             line = " RESOLUTION = "+str(mean(wrange)/fwhm)+"\n"
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
  
  
def fit(xdata, ydata, modelfile, params, bounds,  
        vmicro=1.0, abu=None, vrot=0.0, fwhm=0.0, vmacro=0.0,
        dw=None, strength=1e-4, linelist=linelist0, atom='ap18',
        steprot=0.0, stepfwhm=0.0, lte=None, method='Powell'):
  
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
                 bounds = bounds )                 
  
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
    steprot=steprot, stepfwhm=stepfwhm,  lineid=False, tag=False,  \
    clean=True, save=False, synfile=None, lte=None, compute=True, tmpdir=None)

  chi = np.sum( (ydata - np.interp(xdata, out[0], out[1]/out[2]) )**2)
  
  #plt.clf()
  #plt.plot(xdata,ydata,xdata,np.interp(xdata, out[0], out[1]/out[2]) )
  #plt.show()
  print('chi=',chi)
  
  return(chi)
  
  

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


