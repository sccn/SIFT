
Author : Zhilin Zhang
Date   : Sep 25, 2012
Version: 1.3.4


------------------------------------------------------------------------------
This package includes BSBL-EM, BSBL-BO, and EBSBL-BO algorithms. 
All the three algorithms are introduced in:

[1] Zhilin Zhang, Bhaskar D. Rao, Extension of SBL Algorithms for the Recovery
of Block Sparse Signals with Intra-Block Correlation, submitted to IEEE Trans.
on Signal Processing, 2012. [Online] arXiv:1201.0862v1 [stat.ML]

BSBL-EM is also presented in the conference paper:

[2] Zhilin Zhang, Bhaskar D. Rao, Recovery of Block Sparse Signals Using the 
Framework of Block Sparse Bayesian Learning, ICASSP 2012

with the name 'Cluster-SBL (Type I)'. 



------------------------------------------------------------------------------
In this package, there are four demo files. Their functions are described
as follows:

DEMO_knownPartition_noise.m
    Shows how to use BSBL-EM and BSBL-BO to carry out a noisy experiment 
    when the block partition is known.

DEMO_knownPartition_noiseless.m
    Shows how to use BSBL-EM and BSBL-BO to carry out a noiseless experiment
    when the block partition is known.

DEMO_unknownBlockPartition.m
    Shows how to use the four algorithms to carry out a noisy experiment when
    the block partition is UNKNOWN.

DEMO_nonSparse.m
    Shows hwo to use BSBL-BO to recover a non-sparse signal (Fetal ECG signal)
    Details can be found in the paper:

  [3]Zhilin Zhang, Tzyy-Ping Jung, Scott Makeig, Bhaskar D. Rao, 
     Compressed Sensing for Energy-Efficient Wireless Telemonitoring of 
     Non-Invasive Fetal ECG via Block Sparse Bayesian Learning, 
     submitted to IEEE Trans. on Biomedical Engineering, 2012. 
     [Online] http://arxiv.org/abs/1205.1287



------------------------------------------------------------------------------
There are two folds, which show how to use BSBL-BO to recover fetal ECG and EEG. 
The details are given below:

[Compressed Sensing of FECG]
   includes demo files using BSBL-BO to recover a 8-channel
   fetal ECG raw recordings and then perform ICA decomposition.

[Compressed Sensing of EEG]
   includes demo files using BSBL-BO to recover EEG



------------------------------------------------------------------------------
Welcome to send me any questions. 
I can be reached at: zhangzlacademy@gmail.com


------------------------------------------------------------------------------







