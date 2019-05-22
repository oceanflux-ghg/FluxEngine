#!/bin/bash

#Installation script will install the following dependencies:
# Command Line Developer Tools for xcode. This is required for installation of other dependencies.
# App::cpanminus (perl module)
# AppConfig (perl module)
# pip (for easy installation of netCDF4 on OS X)
# netCDF4 (python module)

echo "This script will install the following python modules if they are not already installed:"
echo "  *  Homebrew"
echo "  *  A seperate version Python 2.7 (as not to interfere with the MacOS version)"
echo "  *  numpy"
echo "  *  scipy"
echo "  *  netCDF4"
echo "  *  pandas"
echo "  *  jupyter (for tutorials)"
echo "  *  matplotlib (for tutorials)"
echo ""
read -p "You may need to enter your password during the installation. Press <enter> to continue the installation or press ctrl+c at anytime to cancel the installation."
echo ""

echo "Checking for Homebrew installation"
if (command -v brew >/dev/null 2>&1)
then
	echo "Homebrew is already installed."
else
	read -p "Homebrew will be installed. Press <enter> to install brew or ctrl+c to cancel installation."
	echo "Installing Homebrew."
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
fi

#Don't want to mess with the macOS version of python so setup a seperate python environment
brew install python2

#This version comes with pip installed, but it needs updating.
pip install --upgrade pip setuptools wheel


#Install numpy
echo "Installing python dependencies."
if (python -c 'import numpy;'>/dev/null 2>&1)
then
	echo "numpy is already installed."
else
	read -p "numpy will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing numpy."
	pip install numpy
fi

#Install netCDF4
echo "Installing python dependencies."
if (python -c 'import netCDF4;'>/dev/null 2>&1)
then
	echo "netCDF4 is already installed."
else
	read -p "netCDF4 will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing netCDF4."
	pip install netCDF4==1.3.1
fi

#Install scipy
echo "Installing python dependencies."
if (python -c 'import scipy.integrate;'>/dev/null 2>&1)
then
	echo "scipy is already installed."
else
	read -p "scipy will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing scipy."
	pip install scipy
fi

echo "Dependencies are now installed. You should be good to go!"

#Install pandas
echo "Installing python dependencies."
if (python -c 'import pandas;'>/dev/null 2>&1)
then
	echo "pandas is already installed."
else
	read -p "pandas will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing pandas."
	pip install pandas
fi

#Install jupyter
echo "Installing python dependencies."
if (python -c 'import jupyter;'>/dev/null 2>&1)
then
	echo "jupyter is already installed."
else
	read -p "jupyter will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing jupyter."
	pip install jupyter
fi

#Install matplotlib
echo "Installing python dependencies."
if (python -c 'import matplotlib;'>/dev/null 2>&1)
then
	echo "matplotlib is already installed."
else
	read -p "matplotlib will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	echo "Installing matplotlib."
	pip install matplotlib
fi

echo "Dependencies are now installed. You should be good to go!"


