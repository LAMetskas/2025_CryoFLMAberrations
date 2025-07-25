# Cryo-FLM-aberrations

Thistoolbox is distributed as accompanying software for manuscript: Hongjia Li, Lauren Ann Metskas and Fang Huang, "Systematic Characterization of Optical Aberrations Reveals Cryo-FLM Localization Fidelity "

# Abstract
Cryo-correlative light and electron microscopy (cryo-CLEM) facilitates in situ imaging and structural analysis by combining the molecular specificity of fluorescence microscopy with the ultrastructural resolution of cryo-electron microscopy. By further combining single molecule localization with cryo-CLEM, molecular positions of individual emitters can be revealed in the context of the electron density map of a cell, providing unique insights to profound questions in cell biology and virology. However, cryogenic fluorescence light microscopy (cryo-FLM) suffers from severe and spatially heterogeneous optical aberrations that distort the point spread function, limiting the accuracy of molecular localizations as well as downstream cryo-transmission electron microscopy workflows. Here, we present a systematic and quantitative analysis of optical aberrations in a commercial cryo-FLM system, uncovering the sources of significant distortions such as system imperfections, refractive index mismatches, and sample-induced heterogeneities. These system and sample induced aberrations lead to localization errors up to 90 nm laterally and over 300 nm axially, challenging the feasibility of precise molecular positioning within the vitrified specimen. We demonstrate that these errors are partially mitigated by spatially matched or adaptive point spread function models pushing the error rate down to ten nanometers or less, offering practical guidance for aberration-aware cryo-FLM and cryo-CLEM strategies. Our findings highlight the necessity of accurate, in situ point spread function modeling to achieve nanometer-scale localization in cryo-FLM. The experimental pipeline developed in this work establishes a novel tool to assess optical performance in cryo-CLEM and cryogenic focused ion beam milling workflows as the field strives toward accurate and precise molecular localization.

# Codes
We provide three MATLAB toolboxes as supplementary software for this manuscript:
1.	MLE-based Phase Retrieval Algorithm – For retrieving point spread function (PSF) models from bead image stacks.
2.	Data Simulation for Biplane Setup – For generating simulated biplane datasets.
3.	3D Localization for Biplane Setup – For high-precision single-molecule localization in biplane microscopy.

# Installation Requirements
•	Operating System: Windows 7 or later (64-bit).
•	MATLAB: R2023b (64-bit) or newer. Download from MathWorks.
•	GPU & Drivers: CUDA 12-compatible graphics driver. Download from NVIDIA CUDA 12.0 Archive.
•	DIPimage Toolbox: Version 3.5.2. Download from MATLAB File Exchange.
