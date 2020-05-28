# robust-mpc-sls
Code needed to reproduce the examples in 

S. Chen*, H.Wang*, M. Morari, V. Preciado, and N. Matni, [Robust Closed-loop Model Predictive Control via System Level Synthesis](https://arxiv.org/pdf/1911.06842.pdf)

In the folder Robust_MPC, SLSDesign.m implements robust SLS MPC, tube MPC and dynamic programming approaches in the paper and CoarseSLSDesign.m implements coarse SLS MPC.

Use the codes in the folder Double_integrator_example to regenerate the results in the paper.

Installation of MPT3 (https://www.mpt3.org/) and Yalmip (https://yalmip.github.io/) required to run the Matlab codes.