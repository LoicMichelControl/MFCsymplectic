    "A variational and symplectic framework for model-free control: preliminary results"
    
    Code associated with the paper accepted at the CCTA 2025 conference.
    
    (c) [2025]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes
    by Loïc MICHEL
    All rights reserved under MIT license.
    
    In the file 'ModelFreeControl_Symplectic_PwrConverter_buck.m', select which case to simulate 
    by setting 'PowerControl', 'VarIntegrator' and the Delta_R perturbation:
    
    -> Fig. 2 : PowerControl = 0, VarIntegrator = 0, Delta_R = 1 
    -> Fig. 3-4 : PowerControl = 0, VarIntegrator = 1, Delta_R = 1 
    -> Fig. 5 : PowerControl = 0, VarIntegrator = 0, Delta_R = 5 
    -> Fig. 6 : PowerControl = 0, VarIntegrator = 1, Delta_R = 5 
    -> Fig. 7 : PowerControl = 1, VarIntegrator = 0, Delta_R = 1 
    -> Fig. 8(top & bottom) : PowerControl = 0 or = 1, VarIntegrator = 1, Delta_R = 1 
       
    The file 'ModelFreeControl_Symplectic_PwrConverter_SweepEpsilon.m' computes the evolution of
    the tracking error with respect to several values of Epsilon_M (Fig. 9).
