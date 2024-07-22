![263416749-1abc1d2d-36bb-4cfb-9328-b57a96044f55](https://github.com/user-attachments/assets/b45a5caa-6b39-4291-b137-125132e5ade0)

## The Source Information Flow Toolbox

Developed and Maintained by: Tim Mullen and Arnaud Delorme (SCCN, INC, UCSD) 2009

SIFT is an EEGLAB-compatible toolbox for the analysis and visualization of
multivariate causality and information flow between sources of
electrophysiological (EEG/ECoG/MEG) activity. It consists of a suite of
command-line functions with an integrated Graphical User Interface for
easy access to multiple features. There are currently six modules: data
preprocessing, model fitting and connectivity estimation, statistical
analysis, visualization, group analysis, and neuronal data simulation.

Methods currently implemented include:

-   Preprocessing routines
-   Time-varying (adaptive) multivariate autoregessive modeling
    -   granger causality
    -   directed transfer function (DTF, dDTF)
    -   partial directed coherence (PDC, GPDC, PDCF, RPDC)
    -   multiple and partial coherence
    -   event-related spectral perturbation (ERSP)
    -   and many other measures...
-   Bootstrap/resampling and analytical statistics
    -   event-related (difference from baseline))
    -   between-condition (test for condition A = condition B)
-   A suite of programs for interactive visualization of information
    flow dynamics across time and frequency (with optional 3D
    visualization in MRI-coregistered source-space).

## Acknowledgements

- Arnaud Delorme was instrumental in the development of the SIFT framework and integration into EEGLAB as well as contributing initial BrainMovie3D code.
- Christian Kothe contributed the arg() framework for function I/O and auto-GUI generation
- Wes Thompson consulted on statistics and methods for bayesian smoothing and multi-subject analysis
- Alejandro Ojeda contributed routines for fast ridge regression

SIFT makes use of routines from (or is inspired by) the following open-source packages:

- [ARFIT](https://github.com/tapios/arfit) (Schneider et al)
- [TSA/Biosig](http://octave.sourceforge.net/tsa/) (Schlögl et al)
- [Chronux](https://chronux.org) (Mitra et al)
- [DAL/SCSA](https://ttic.uchicago.edu/~ryotat/softwares/dal/) (Tomioka / Haufe et al)
- [BCILAB](http://sccn.ucsd.edu/wiki/BCILAB) (Kothe et al)


## Official Website

[SIFT page in the SCCN wiki](http://sccn.ucsd.edu/wiki/SIFT)

## Citation

If you find this toolbox useful for your research, PLEASE include the following citations with any publications and/or presentations which make use of SIFT:

1. Mullen, T. R. (2014). The dynamic brain: Modeling neural dynamics and interactions from human electrophysiological recordings (Order No. 3639187). Available from Dissertations & Theses @ University of California; ProQuest Dissertations & Theses A&I. (1619637939)
2. Delorme, A., Mullen, T., Kothe C., Akalin Acar, Z., Bigdely Shamlo, N., Vankov, A., Makeig, S. (2011) "EEGLAB, SIFT, NFT, BCILAB, and ERICA: New tools for advanced EEG/MEG processing." Computational Intelligence and Neuroscience vol. 2011, Article ID 130714, 12 pages.

## License

SIFT is licensed under the GPL-2, see LICENSE.txt
ANY USE OF SIFT IMPLIES THAT YOU HAVE READ AND AGREE WITH THE TERMS AND CONDITIONS OF THE SIFT LICENSE AS STATED BELOW:

## ADDITIONAL NOTE

SIFT is designed and distributed for research purposes only. SIFT should not be used for medical purposes. The authors accept no responsibility for its use in this manner.

## Verions

v1.6 - fix conflict with BrainMovie plugin. Fix minor GUI issues.
