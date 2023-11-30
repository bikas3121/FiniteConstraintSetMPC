Python implementation of the Finite Constraint Set MPC (Moving horizon Optimal Quantizers)
* Reference: 
    ** Moving horizon optimal quantizer for audio signal, Goodwin et al. 2003
    ** Finite constraint set receding horizon quadratic control,  Quevedo et al. 2003
The python implementation of the MPC is done in the context of the quantization. The quantization is the 
process of the mapping an analog (continouous) signal into finite discrete set to represent them with binary words. 
The quantization is a very important in digital signal processing. The quantization process causes of the loss of 
the information and thus introduces distortion into the signal upon reconstruction. This MPC based optimal control 
technique minimized the error between the reference signal and the reconstructed signal. 

The MPC algorithm is implemented in the python. The dynamics of the quantizers is modelled as linear mixed integer program 
and solved using GUROBI solver. 

 

