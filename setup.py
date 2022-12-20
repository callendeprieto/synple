import os
import sys
import shutil
from setuptools import setup, find_packages, find_namespace_packages
from setuptools.command.install import install
import subprocess

def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
        return sys.prefix

    return None

def get_bin_path():
    # Get environment bin/ directory
    bindir = None
    # check if the --user option was set
    userdir = None
    for u in install.user_options:
        if u[0]=='user':
            # "install in user site-package '/Users/nidever/.local/lib/python3.7/site-packages'"
            uline = u[2]
            userdir = uline[uline.find('site-package')+13:]
            userdir = userdir.replace("'","")
    if userdir is not None:
        # /Users/nidever/.local/bin
        bindir = os.path.dirname(os.path.dirname(os.path.dirname(userdir)))+'/bin'
        if os.path.exists(bindir)==False:
            bindir = False
    # Try virtual environment using sys.prefix
    if bindir is None:
        venv = get_virtualenv_path()
        if venv is not None:
            bindir = venv+'/bin'
            if os.path.exists(bindir)==False:
                bindir = None
    # Get bin/ directory from python executable
    if bindir is None:
        out = subprocess.check_output('which python',shell=True)
        if type(out) is bytes:
            out = out.decode()
        bindir = os.path.dirname(out)
        if os.path.exists(bindir)==False:
            bindir = None

    if bindir is None:
        raise Exception('No bin/ directory found')

    return bindir
    
def compile_and_install_software():
    """Used the subprocess module to compile/install the C software."""
    src_path = './synspec/'
    
    # Install the software
    subprocess.check_call('make install', cwd=src_path, shell=True)

    # Get the path for the binaries
    bindir = get_bin_path()
        
    # Copy fortran binaries to bin/ directory
    for f in ['rotin','synspec54','list2bin']:
        if os.path.exists(bindir+f): os.remove(bindir+f)
        print('Copying bin/'+f+' -> '+bindir+'/'+f)
        shutil.copyfile('bin/'+f,bindir+'/'+f)
        
    # Download and convert linelists
    curdir = os.path.abspath(os.curdir)
    subprocess.check_call('make', cwd=curdir+'/linelists', shell=True)    
    
class CustomInstall(install):
    """Custom handler for the 'install' command."""
    def run(self):
        compile_and_install_software()
        super().run()

setup(name='synple',
      version='1.0',
      description='An Easy-to-Use Python Wrapper for the Spectral Synthesis Code Synspec',
      author='Carlos Allende Prieto',
      author_email='callende@iac.es',
      url='https://github.com/callendeprieto/synple',
      scripts=['bin/synple'],
      requires=['numpy','astropy(>=4.0)','scipy'],
      cmdclass={'install': CustomInstall},
      zip_safe = False,
      include_package_data=True,
      packages=find_namespace_packages(where="python"),
      package_dir={"": "python"}      
)

