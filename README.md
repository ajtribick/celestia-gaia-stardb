Gaia DR2 for celestia.Sci
=========================

This repository contains Python scripts to generate a celestia.Sci star
database from the Gaia DR2 data.

In order to limit the download size required and to maintain compatibility
with the use of HIP/TYC2 identifiers as the primary key for stars in
celestia.Sci, only Gaia data for HIP and TYC2 stars is processed.

## Prerequisites

-  Internet connection for downloading the data
-  Gaia archive account (https://gea.esac.esa.int/archive/)
-  Python 3
-  celestia.Sci

## Folder contents

-  `download_data.py`: script to download the data files

## How to use

1.  Clone or download this repository
2.  Open a command window in the repository directory
3.  Set up a Python 3 virtual environment
    - Windows: py -3 -m venv myenv
    - Linux: python3 -m venv myenv
4.  Switch to the virtual environment and install the requirements
    - Windows (Powershell): .\myenv\bin\Activate.ps1
    - Windows (cmd): myenv\bin\activate.bat
    - Linux: source myenv/bin/activate
5.  Install the requirements:
    - Windows: py -m pip install -r requirements.txt
    - Linux: python -m pip install -r requirements.txt
6.  Run the download script. You will need your Gaia archive login. THIS STEP
    MAY TAKE SEVERAL HOURS
    - Windows: py download_data.py
    - Linux: python download_data.py
