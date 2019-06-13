#!/bin/bash

#Installation script will install the following dependencies:
# Command Line Developer Tools for xcode. This is required for installation of other dependencies.
# App::cpanminus (perl module)
# AppConfig (perl module)
# pip (for easy installation of netCDF4 on OS X)
# netCDF4 (python module)

echo "This script will install the following python modules if they are not already installed:"
echo "  *  numpy"
echo "  *  scipy"
echo "  *  netCDF4"
echo "  *  pandas"
echo "  *  jupyter"
echo "  *  matplotlib"
echo ""
echo "Pip will be installed unless it is already present."
echo ""
#read -p "You may need to enter your password during the installation. Press <enter> to continue the installation or press ctrl+c at anytime to cancel the installation."
#echo ""


#Install numpy
echo "Installing python dependencies."
if (python -c 'import numpy;'>/dev/null)
then
	echo "numpy is already installed."
else
	read -p "numpy will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo easy_install pip
	fi
	echo "Installing numpy."
	sudo pip install numpy
fi

#Install netCDF4
echo "Installing python dependencies."
if (python -c 'import netCDF4;'>/dev/null)
then
	echo "netCDF4 is already installed."
else
	read -p "netCDF4 will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo easy_install pip
	fi
	echo "Installing netCDF4."
	sudo pip install netCDF4
fi

#Install scipy
echo "Installing python dependencies."
if (python -c 'import scipy.integrate;'>/dev/null)
then
	echo "scipy is already installed."
else
	read -p "scipy will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo apt-get update
		sudo apt-get install python-pip -y
	fi
	echo "Installing scipy."
	sudo pip install scipy
fi

#Install pandas
echo "Installing python dependencies."
if (python -c 'import pandas;'>/dev/null)
then
	echo "pandas is already installed."
else
	read -p "pandas will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo apt-get update
		sudo apt-get install python-pip -y
	fi
	echo "Installing pandas."
	sudo pip install pandas
fi

#Install jupyter
echo "Installing python dependencies."
if (python -c 'import jupyter;'>/dev/null)
then
	echo "jupyter is already installed."
else
	read -p "jupyter will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo apt-get update
		sudo apt-get install python-pip -y
	fi
	echo "Installing jupyter."
	sudo pip install jupyter
fi

#Install matplotlib
echo "Installing python dependencies."
if (python -c 'import matplotlib;'>/dev/null)
then
	echo "matplotlib is already installed."
else
	read -p "matplotlib will be installed using pip. Press <enter> to continue or ctrl+c to cancel installation."
	if ! command -v pip >/dev/null
	then
		read -p "pip must be installed first. Press <enter> to install pip or ctrl+c to cancel installation."
		echo "Installing pip."
		sudo apt-get update
		sudo apt-get install python-pip -y
	fi
	echo "Installing matplotlib."
	sudo pip install matplotlib
fi

echo "Dependencies are now installed. You should be good to go!"


