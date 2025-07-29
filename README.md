# Cryo-FLM-aberrations

Thistoolbox is distributed as accompanying software for manuscript: Hongjia Li, Lauren Ann Metskas and Fang Huang, "Systematic Characterization of Optical Aberrations Reveals Cryo-FLM Localization Fidelity "


# Abstract

Cryo-correlative light and electron microscopy (cryo-CLEM) facilitates in situ imaging and structural analysis by combining the molecular specificity of fluorescence microscopy with the ultrastructural resolution of cryo-electron microscopy. By further combining single molecule localization with cryo-CLEM, molecular positions of individual emitters can be revealed in the context of the electron density map of a cell, providing unique insights to profound questions in cell biology and virology. However, cryogenic fluorescence light microscopy (cryo-FLM) suffers from severe and spatially heterogeneous optical aberrations that distort the point spread function, limiting the accuracy of molecular localizations as well as downstream cryo-transmission electron microscopy workflows. Here, we present a systematic and quantitative analysis of optical aberrations in a commercial cryo-FLM system, uncovering the sources of significant distortions such as system imperfections, refractive index mismatches, and sample-induced heterogeneities. These system and sample induced aberrations lead to localization errors up to 90 nm laterally and over 300 nm axially, challenging the feasibility of precise molecular positioning within the vitrified specimen. We demonstrate that these errors are partially mitigated by spatially matched or adaptive point spread function models pushing the error rate down to ten nanometers or less, offering practical guidance for aberration-aware cryo-FLM and cryo-CLEM strategies. Our findings highlight the necessity of accurate, in situ point spread function modeling to achieve nanometer-scale localization in cryo-FLM. The experimental pipeline developed in this work establishes a novel tool to assess optical performance in cryo-CLEM and cryogenic focused ion beam milling workflows as the field strives toward accurate and precise molecular localization.



# Codes

We provide three MATLAB toolboxes as supplementary software for this manuscript:

1. MLE-based Phase Retrieval Algorithm - For retrieving point spread function (PSF) models from bead image stacks.
2. Data Simulation for Biplane Setup - For generating simulated biplane datasets.
3. 3D Localization for Biplane Setup - For high-precision single-molecule localization in biplane microscopy.



# Installation Requirements

* Operating System: Windows 7 or later (64-bit).
* MATLAB: R2023b (64-bit) or newer. Download from MathWorks.
* GPU \& Drivers: CUDA 12-compatible graphics driver. Download from NVIDIA CUDA 12.0 Archive.
* DIPimage Toolbox: Version 3.5.2. Download from MATLAB File Exchange.



# Dataset and source codes

1\. Demonstration Datasets

   We provide example datasets for testing and demonstration:

* 1\. Phase Retrieval\\3. Test Data\\1. Bead Stack RT\\Beads\_stack.mat: Bead image stack collected under room-temperature conditions.
* 1\. Phase Retrieval\\3. Test Data\\2. Bead Stack Cryo\\Beads\_stack.mat: Bead image stack collected under cryogenic-temperature conditions.



2\. Phase Retrieval Source Codes

   Folder: 1. Crop Single Beads

* Crop\_single\_beads\_1p.m – Crops subregions (64 × 64 × n) containing a single bead stack for phase retrieval.
* Crop\_single\_noise\_1p.m – Crops subregions (64 × 64 × n) containing only background regions for noise estimation.

   Folder: 2. Phase Retrieval

* main.m – Main script for running MLE-based phase retrieval.
* visualize\_all\_results.m – Visualizes the retrieved results, including measured vs. fitted PSFs, pupil phase, and evaluation metrics.
* calc\_nLogLikelihood.m – Calculates the negative log-likelihood used in phase retrieval optimization.
* Support Folder – Contains all auxiliary functions, including PSF generation, pupil modeling, and other helper utilities.



3\. Data Simulation Source Codes

   Folder: 2. Data Simulation Biplane

* main.m – Generates biplane simulated datasets with specified photon numbers and aberrations.



3\. 3D Localization Source Codes

   Folder: 3. 3D Localization Biplane

* main.m – Performs GPU-based 3D localization for simulated PSFs.
* genIniguess.m – Estimates initial lateral positions.
* geniniBiplane\_z\_mat\_parfor.m – Estimates initial axial positions.
* genpsf\_biplane\_real.m – Generates PSF models directly from pupil functions.
* gensamplepsf\_biplane.m – Pre-generates channel-specific PSF models.
* cc2.m – Calculates 2D cross-correlation.
* catstruct.m – Concatenates results from the two biplane channels.
* genpsfstruct.m – Computes first-order partial derivatives of a 3D PSF for MLE-based optimization.



# Citation



# Acknowledgments

The authors thank Louise Bertrand from Leica Microsystems, Inc. for providing technical and applications support. This work was supported by NIAID award 1DP2AI164293-01 to L.A.M. and NIGMS MIRA award R35GM119785 to F.H.

# Copyrights
Users agree to use the script package as is.

This package is provided under Creative Commons copyright license Attribution-NonCommercial-ShareAlike 4.0 International (https://creativecommons.org/licenses/by-nc-sa/4.0/).

This copyright license allows use for noncommercial use only, with attribution to the original package and sharing of all changes under the same license as the original. 
