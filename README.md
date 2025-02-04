    "A variational and symplectic framework for model-free control: preliminary results"
    
    Code associated to the work submitted to the CCTA'25 conference
    
    (c) [2025]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes
    by Loïc MICHEL
    All rights reserved under MIT license.
    
    In the file 'ModelFreeControl_Symplectic_PwrConverter_buck.m', select which case to simulate 
    by setting 'PowerControl', 'VarIntegrator' and the Delta_R perturbation:
    
     -> Fig. 2 : PowerControl = 0, VarIntegrator = 0, Delta_R = 1 (l.68)
     -> Fig. 3-4 : PowerControl = 0, VarIntegrator = 1, Delta_R = 1 (l.68)
     -> Fig. 5 : PowerControl = 0, VarIntegrator = 0, Delta_R = 5 (l.67)
     -> Fig. 6 : PowerControl = 0, VarIntegrator = 1, Delta_R = 5 (l.67)
     -> Fig. 7 : PowerControl = 1, VarIntegrator = 0, Delta_R = 1 (l.68)
     -> Fig. 8-9 : PowerControl = 1, VarIntegrator = 1, Delta_R = 1 (l.68)
     
     The file 'ModelFreeControl_Symplectic_PwrConverter_SweepGamma.m' computes both the evolution of gamma(t) 
     and the tracking error with respect to several values of Gamma (Figs. 10 and 11).
     
     The file 'ModelFreeControl_Symplectic_PwrConverter_SweepEpsilon.m' computes both the evolution of gamma(t) 
     and the tracking error with respect to several values of Epsilon_M (Figs. 12 and 13).
