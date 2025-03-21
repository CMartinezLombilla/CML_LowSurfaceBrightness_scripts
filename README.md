# cml_lsb_utils

`cml_lsb_utils` is a Python package to help astronomers to **analyse low surface brightness images** of astronomical objects such as galaxies of groups of galaxies. The tools should be used to analyse outer galaxy features, galaxy discs, intra-group light (IGL), or intra-cluster light (ICL), as they allow the user to get photometry and extract different types of surface brightness profiles. `cml_lsb_utils` provides a variety of functions and classes to perform tasks such as converting counts to surface brightness, extracting radial profiles using different aperture geometries, fitting isophotes, interactively masking images, and converting between units.

## Features

- **Counts to Surface Brightness Conversion:**  
  Convert image counts to surface brightness with customizable zero points, extinction, dimming, and K-corrections.

- **Radial Profile Extraction:**  
  Extract surface brightness profiles using rectangular, circular, and elliptical apertures. Includes support for shifted and projected profiles.

- **Isophote Fitting:**  
  Fit elliptical isophotes to images to analyze galaxy structures and determine parameters like centroids, ellipticities, and position angles.

- **Interactive Masking:**  
  Tools to interactively delete or add masked regions to images using GUI prompts.

- **Unit Conversions:**  
  Functions to convert between pixels, arcseconds, and kiloparsecs, as well as to calculate magnitudes and luminosities.

- **PSF Scaling and Subtraction:**  
  Utilities to scale and subtract PSF light from stellar images.

## Project Structure

This repository contains:

1. **Core package:** 
  - **Script:** `cml_lsb_utils` usage and explanation explained below. 
  
Two additional directories containing scripts for specific analyses (****DISCLAIMER**: these two directories haven't been updated neither revised ins a few years now and they might contain comments in both English and Spanish):

2. **ICL Analysis:**  
   - **Directory:** `/ICL_LSST_Mock`  
   - **Main Subfolder:** `/ICL_LSST_Mock/Codes_for_ICL_analysis`  
   These scripts are designed to extract and analyze the Intra-Cluster Light (ICL) in LSST mock images of clusters of galaxies. The scripts are numerically prefixed (e.g., `1_Image_ready`) to indicate the recommended execution order.

3. **IGL Analysis:**  
   - **Directory:** `/IGL-HSC`  
   This directory contains scripts for extracting and analyzing the Intra-Group Light (IGL) in groups of galaxies using HSC PDR2 images. There are folders for three different methods to extract the IGL; however, the fully completed and reviewed method from our research is the one based on the distance map (located in `/Method-DistancesMap`). The subfolders and scripts here are also numerically prefixed to suggest the proper execution sequence. 


## Contributors

This project was developed years ago by myself and Felipe. It represents work we did during our research, and we do not plan to actively improve or maintain it. The code is provided as-is to help other astronomers and researchers with the tools we needed at that time.

## Disclaimer

This software is provided "as-is", without any warranty or guarantee of support. Use it at your own risk. We make no guarantees that the code is bug-free or suitable for your research purposes. The project is not actively maintained, and no further improvements will be made.

## Requirements

- Python 3.x
- numpy
- matplotlib
- astropy
- photutils
- scipy
- tkinter (for interactive GUI features)

Install the required dependencies using pip:

```bash
pip install -r requirements.txt
```

## Installation

Clone the repository or download the project files, and ensure that the `cml_lsb_utils.py` file and the `__init__.py` file are in the same package directory. For example, if the package folder is named `cml_lsb_utils`, you can install or add it to your Python path:

```bash
git clone https://github.com/yourusername/cml_lsb_utils.git
cd cml_lsb_utils
```

## Usage

Once installed and on your Python path, you can import the functions and classes directly. For example:

```python
from cml_lsb_utils import Counts2Sb, CentralProfile

# Example: Convert counts to surface brightness
zp = 25.0         # zero point in mag/arcsec^2
pixscale = 0.4    # pixel scale in arcsec/pixel
counts = 1000.0

c2sb = Counts2Sb(zp, pixscale)
surface_brightness = c2sb.do(counts)
print("Surface brightness:", surface_brightness)

# Example: Create a central profile object
# (Assume you have an image array 'data' and center coordinates Xc, Yc)
profile = CentralProfile(data, Xc=100, Yc=100, nbins=10, npix=200, height=5, zp=zp, pixscale=pixscale)
```

Refer to the inline documentation within the code for detailed explanations of parameters and return values.

## Documentation

Each function and class in this package includes detailed docstrings. Use Python’s built-in `help()` function or consult the source code for more information.

## Contributing

Since this work was completed years ago by myself (Cristina Martinez-Lombilla) and [Felipe Jimenez-Ibarra](https://github.com/felipeji) and is provided "as-is" without plans for future improvements, contributions are not expected. However, if you have suggestions or improvements that may benefit other researchers, feel free to fork the repository.

## License

No license is provided for this project. All rights are reserved by the contributors.

## Contact

For questions or further information, please contact [crismarlom@gmail.com].
