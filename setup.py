#!/usr/bin/env python

import os
import shutil
from setuptools import setup, find_packages, extension
from setuptools.command.install import install
#from numpy.distutils.core import setup, Extension


class PreInstallCommand(install):
        """Pre-installation for installation mode."""
        def run(self):
                # Step 1. Compile fortran code
                os.system('cd synspec; make clean; make; cd ..')
                # Step 2. Download the linelist
                
                
class PostInstallCommand(install):
        """Post-installation for installation mode."""
        def run(self):
                # Step 3. Copy binaries to scripts directory
                
                synspec = 'bin/synspec54'
                if os.path.exists(synspec):
                        self.copy_file(synspec,self.install_scripts)
                rotin = 'bin/rotin'
                if os.path.exists(rotin):
                        self.copy_file(rotin,self.install_scripts)
                

#mod = extension.Extension(name = 'synspec', sources = ['synspec/synspec54.f'])
setup(name='synple',
      version='1.2',
      description='An Easy-to-Use Python Wrapper for the Spectral Synthesis Code Synspec',
      author='Carlos Allende Prieto',
      author_email='callende@iac.es',
      url='https://github.com/callendeprieto/synple',
      packages=find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3.6',
      #ext_modules = [mod],
      scripts=['bin/synple'],
      requires=['gfortran','wget','numpy','astropy(>=4.0)','scipy'],
      include_package_data=True,
      #cmdclass={'install': PostInstallCommand},
      #package_data={'synple': ['data/']},
)
