# Repeated-quantum-backflow-and-quantum-overflow spectral estimates

## Description

- This repository contains the approximations for the maximum and minimum of the spectra of specific backflow operators $C^{(M)}$ ($1\le M\le 4$) described in section 5 of (https://arxiv.org/abs/2505.13184) and the code that generated them.
- Format: for $1\le m\le 4$, specestsN500Mm.csv contains the estimates for $C^{(m)}$ as a csv file giving values keeping 150 digits after the decimal point. 
  - Each line after the first has format $N$ $b$ $o$, where $b = \lambda_{back}^{(m)}(N)$ and $o = \lambda_{over}^{(m)}(N)$, with $1\le N\le 500$.
- Alongside the numerical estimates, a Maple worksheet (spectral_estimate_reader.mw) is provided to read in the estimates for use.

## Eigenvalue generation code
- Prerequisites:
    Python packages: mpmath,flint,math,time,sys.
    Maple 2024 or higher.

- To generate your own eigenvalue estimates, clone or download the /code folder. The .py sheet should be run first to generate the matrix elements and then the .mpl file to generate the eigenvalue estimates. Instructions for each are contained in the source code.

## Contact

- Christopher J. Fewster: [chris.fewster@york.ac.uk]
- Harkan J. Kirk-Karakaya: [harkan.kirk-karakaya@york.ac.uk]


Date: July 2025

