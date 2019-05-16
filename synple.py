#!/usr/bin/python3
'''
Python wrapper for synspec 
'''
import os
import sys
import subprocess
import numpy as np
import time
import matplotlib.pyplot as plt


#configuration
synspecdir = "/home/callende/synspec"

#relative paths
modelatomdir = synspecdir + "/data"
linelistdir = synspecdir + "/linelists"
bindir = synspecdir + "/bin"
synspec = bindir + "/synspec53p"
rotin = bindir + "/rotin3"

#other stuff
bolk = 1.38054e-16 # erg/K
zero = " 0 "
one =  " 1 "
two =  " 2 "

#def mksyn(teff,logg,mh,am,cm,nm,wrange=[15100.,17000],dw=0.05,vmicro=2.0,solarisotopes=False,elemgrid='',welem=None,
#    els=None,atmod=None,kurucz=True,atmosroot=None,atmosdir=None,nskip=0,fill=True,
#    linelist='20150714',h2o=None,linelistdir=None,atoms=True,molec=True,
#    save=False,run=True)

def syn(modelfile, wrange, dw=1e-2, strength=1e-4, vmicro=None, abu=None, \
    linelist=['gfallx3_bpo.19','kmol3_0.01_30.20'],hhm=False, vrot=0.0, fwhm=0.0, \
    steprot=0.0, stepfwhm=0.0, save=False):
#
# high level routine to compute a synthetic spectrum
# input:  modelfile (str) filename for the model atmosphere
#         wrange (2-element tuple) wavelength range (angstroms)
#         dw (float) wavelength step for the calculations (angstroms)
#         strength (float) minimum estimated ratio line/continuum opacity for a line to be considered
#         vmicro (float) microturbulence velocity (km/s)
#         abu (99-element float array) abundances (number density relative to hydrogen)
#         linelist (list of strings) list of line lists to use (1st must be atomic, the rest molecular)
#         hhm (boolean) when True metal continuum opacities are ignored (H, H-, and molecules considered)
#
# output: wave (float array) wavelengths (angstroms)
#         flux (float array) flux (H_lambda, erg/cm2/s/A)
#         cont (float array) flux when no line opacities are considered (H_lambda)
#

  atmostype = identify_atmostype(modelfile)
  if atmostype == 'kurucz':
    teff, logg, vmicro2, abu2, nd, atmos = read_kurucz_model(modelfile) 
  if atmostype == 'marcs':
    teff, logg, vmicro2, abu2, nd, atmos = read_marcs_model(modelfile)
  if atmostype == 'phoenix':
    teff, logg, vmicro2, abu2, nd, atmos = read_phoenix_model(modelfile)
    
  if vmicro == None: vmicro = vmicro2
  if abu == None: abu = abu2

  checksynspec(linelist)  
  checkinput(wrange, modelfile, vmicro, linelist)

  print ('teff,logg,vmicro=',teff,logg,vmicro)
  #print ('abu=',abu)
  #print (len(abu))
  #print ('nd=',nd)
  cleanup()
  
  writetas('tas',nd,linelist)                           #non-std param. file
  write5(teff,logg,abu,hhm)                               #abundance/opacity file
  write8(teff,logg,nd,atmos,atmostype)                  #model atmosphere
  write55(wrange,dw,strength,vmicro,linelist,atmostype) #synspec control file
  create_links(modelfile,linelist)                      #auxiliary data


  synin = open('fort.5')
  synout = open('syn.log','w')
  start = time.time()
  p = subprocess.Popen([synspec], stdin=synin, stdout = synout)
  p.wait()
  synout.flush()
  wave, flux = np.loadtxt('fort.7', unpack=True)
  wave2, flux2 = np.loadtxt('fort.17', unpack=True)
  cont = np.interp(wave, wave2, flux2)
  end = time.time()
  print('syn ellapsed time ',end - start, 'seconds')

  if fwhm > 0. or vrot > 0.:
    start = time.time()
    flux = convol (wave, flux, cont, vrot, fwhm, dw, steprot, stepfwhm, save=True, reuseinputfiles=True)
    end = time.time()
    print('convol ellapsed time ',end - start, 'seconds')

  if save == False: cleanup()

  return(wave, flux, cont)

def convol(wave, flux, cont, vrot, fwhm, dw=1e-2, steprot=0.0, stepfwhm=0.0, save=False, reuseinputfiles=False):

  if reuseinputfiles == False:
    f = open('fort.7','w')
    for i in len(wave):
      f.write( '( %f %f \n)' % (wave[i], flux[i]) )
    f.close()

  f = open('fort.5','w')
  f.write( ' %s %s %s \n' % ("'fort.7'", "'fort.17'", "'fort.11'") )
  f.write( ' %f %f %f \n' % (vrot, dw, steprot) )
  f.write( ' %f %f \n' % (fwhm, stepfwhm) )
  print('stepfwhm=',stepfwhm)
  f.write( ' %f %f %i \n' % (np.min(wave), np.max(wave), 0) )
  f.close()

  synin = open('fort.5')
  synout = open('syn.log','a')
  p = subprocess.Popen([rotin], stdin=synin, stdout = synout)
  p.wait()
  synout.flush()
  wave2, flux2 = np.loadtxt('fort.11', unpack=True)
  print(len(wave),len(wave2))
  flux2 = np.interp(wave, wave2, flux2)
  
  if save == None: cleanup()

  return(flux2)


def identify_atmostype(modelfile):

  if ('PHOENIX' in modelfile and 'fits' in modelfile): atmostype = 'phoenix'
  else: 
    f = open(modelfile,'r')
    line = f.readline()
    if ('TEFF' in line): atmostype = 'kurucz'
    else: atmostype = 'marcs'
    f.close()
   
  return(atmostype)

def checksynspec(linelist):
#checking executable and data are where it should be

  dirs = [synspecdir,modelatomdir,linelistdir,bindir]
  for entry in dirs: assert (os.path.isdir(entry)), 'dir '+entry+' missing'

  files = [synspec,rotin]
  for entry in linelist: files.append(os.path.join(linelistdir,entry))
  for entry in files: assert (os.path.isfile(entry)), 'file '+entry+' missing'

  return()


def checkinput(wrange, modelfile, vmicro, linelist):
#checking input parameters from user

  #find range of atomic line list
  minlambda, maxlambda = getlinelistrange(os.path.join(linelistdir,linelist[0]))

  #check
  assert (wrange[0] > minlambda and wrange[1] < maxlambda),'wrange exceeds the allow range ('+str(minlambda)+' to '+str(maxlambda)+')'
  assert (os.path.isfile(modelfile)),'model atmosphere file '+modelfile+' missing'
  assert (vmicro >= 0.0),'vmicro = '+str(vmicro)+' but cannot < 0.'
  
  return()

def getlinelistrange(atomiclinelist):
#finds out min and max wavelengths for a line list

  f = open(atomiclinelist,'r')
  line = f.readline()
  entries = line.split()
  minlambda = float(entries[0])*10.
  f.seek(os.path.getsize(atomiclinelist)-103)
  line = f.readline()
  f.close()
  entries = line.split()
  maxlambda = float(entries[0])*10.

  return(minlambda,maxlambda)



def writetas(filename,nd,linelist):
#write non-std input parameters
# input: filename -- str -- name of the non-std. param. file to print
#        nd -- int -- number of layers in the model
#        nd -- list -- names of the linelist files (atomic first, then one 
#				or more molecular ones
  
  f = open(filename,'w')
  f.write("ND= "+str(nd)+" \n")
  if len(linelist) > 1:  f.write("IFMOL= "+one+" \n")
  f.close()

  return()


def write55(wrange,dw=1e-2,strength=1e-4,vmicro=0.0,linelist=['gfallx3_bpo.19','kmol3_0.01_30.20'], atmostype='kurucz'):

  imode = 0                        # default, atoms and molecules, at least 2 line lists
  if len(linelist) == 0: imode = 2  # no atomic or molecular line list -> pure continuum and no molecules


  #imode,idst,iprin
  #inmod,zero,ichang,ichemc
  #lyman,zero,zero,zero,zero
  #one,nlte,icontl,zero,ifhe2
  #ihydpr,ihe1pr,ihe2pr
  #wstart,wend,cutoff,zero,strength,wdist 

  if atmostype == 'tlusty': inmod = 1 
  else: inmod = 0

  f = open('fort.55','w')
  f.write(" "+str(imode)+" "+2*zero+"\n")
  f.write(" "+str(inmod)+3*zero+"\n")
  f.write(5*zero+"\n")
  f.write(one+4*zero+"\n")
  f.write(two+2*zero+"\n")
  f.write( ' %f %f %f %i %e %f \n ' % (wrange[0], wrange[1], 200., 0, strength,dw) )
  ll = len(linelist)
  if ll < 2: f.write(2*zero)
  else: f.write(str(ll-1) + ' ' + ' '.join(map(str,np.arange(ll-1)+20)))
  f.write("\n")
  f.write( ' %f  \n' % (vmicro) )
  f.close()

def write5(teff,logg,abu,hhm=False):

  symbol, mass, sol = elements()

  f = open('fort.5','w')
  f.write(' '+str(teff)+" "+str(logg)+"       ! TEFF, GRAV \n")
  f.write(" T  F               ! LTE, GRAY \n")
  f.write(" 'tas'              ! name of non-standard flags \n")
  f.write(" 50                 ! frequencies \n")
  f.write(" "+str(len(abu))+"        ! NATOMS \n")

  ex = np.ones(len(abu))
  if hhm : 
    zex = [1]  #atomic numbers of elements to be included explicitly (contributing cont. opacity)
  else: 
    zex = [1,2,6,7,8,11,12,13,14,20,26]

  for i in zex: ex[i-1] = 2

  for i in range(len(abu)):
    f.write(' %i %e %i %s\n' %  (ex[i], abu[i], 0, '  ! ' +symbol[i]) )

  for i in range(3): f.write("* \n")
  
  if hhm == True:
    f.write("   1   -1     1      0     0     1    ' H 1' 'data/hm.dat' \n" )
    f.write("   0    0     3      0 \n")
    f.write("   1    0     9      0     0     0    ' H 1' 'data/h1s.dat'  \n")
    f.write("   1    1     1      1     0     0    ' H 2' ' '  \n")
    f.write("   0    0     0     -1     0     0    '    ' ' '  \n")
  else:
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
  f.close()

def write8(teff, logg, nd, atmos, atmostype):

  if atmostype == 'tlusty':
    f = open('fort.8','w')
    f.write(" "+str(nd)+" "+str(3)+"\n")
    for i in range(nd):
      f.write(' %e ' % atmos['dm'][i])
    f.write("\n")
    for i in range(nd):
      f.write( '%f %e %e \n' % (atmos['t'][i], atmos['ne'][i], atmos['rho'][i] ) )
    f.close()

  else:
    f = open('fort.8','w')
    f.write( 'TEFF %7.0f  GRAVITY %7.5f  LTE \n' % (teff, logg) )
    for i in range(21): f.write('\n')
    f.write( 'READ DECK6%3i RHOX,T,P,XNE \n' % nd )
    for i in range(nd): 
      f.write( '%e %f %e %e \n' % (atmos['dm'][i], atmos['t'][i], atmos['p'][i], atmos['ne'][i]) )
    f.close()

  return()
  

def create_links(modelfile,linelist):
#create soft links for line lists, model atom dir and model atmosphere

  for i in range(len(linelist)): 
    if i == 0: os.symlink(os.path.join(linelistdir,linelist[0]),'fort.19')
    else: os.symlink(os.path.join(linelistdir,linelist[i]),'fort.'+str(20-1+i))

  os.symlink(modelatomdir,'./data')
  #os.symlink(modelfile,'fort.8')

  return()

def cleanup():
#cleanup all temporary files

  files = os.listdir('.')
  for entry in files: 
    if os.path.islink(entry) and entry.startswith('fort'): os.unlink(entry)
    if os.path.isfile(entry) and entry.startswith('fort'): os.remove(entry)

  if os.path.islink('data'): os.unlink('data')
  if os.path.isfile('tas'): os.remove('tas')

  return()


def read_kurucz_model(modelfile):

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

  dm = [ float(entries[7]) ]

  for i in range(nd-1):
    line = f.readline()
    entries = line.split()

    dm.append(  float(entries[7]))

  atmos = np.zeros(nd, dtype={'names':('dm', 't', 'p','ne'),
                          'formats':('f', 'f', 'f','f')}) 
  atmos['dm'] = dm
  atmos['t'] = t
  atmos['p'] = p
  atmos['ne'] = ne

  return (teff,logg,vmicro,abu,nd,atmos)


def read_phoenix_model(modelfile):

  from astropy.io import fits

  h = fits.open(modelfile)[0].header
  f = fits.open(modelfile)[1].data

  nd = len(f['temp'])

  teff = float(h['PHXTEFF'])
  logg = float(h['PHXLOGG'])
  vmicro = float(h['PHXXI_L'])

  m_h = float(h['PHXM_H'])
  alpha = float(h['PHXALPHA'])
  
  symbol, mass,sol = elements(husser=True) 
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

  vmicro = 0.0
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



def elements(husser=False):

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

  if not husser:
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
	      
    sol[0] = 1.

  else:
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

  sol[0] = 1.
  for i in range(len(sol)-1): sol[i+1] = 10.**(sol[i+1]-12.0)

  return (symbol,mass,sol)

if __name__ == "__main__":

  modelfile = sys.argv[1]
  wstart = float(sys.argv[2])
  wend = float(sys.argv[3])

  #symbol, mass, sol = elements()
  x,y,yc=syn(modelfile,[wstart,wend],  h=True, linelist=['gfallx3_bpo.19'], fwhm=0.2, vrot=10.)
  #x2,y2,yc2=syn('sun.mod',[4300,4400],linelist=['gfallx3_bpo.19'])
  plt.plot(x,y/yc,'b')
  plt.show()

