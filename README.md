# RGDFoam_v0


This repository is used to store the first version of RGDFoam with look up table mechanism.

RGDFoam main solver is use
     p = L(rho,e)
     h = L(rho,e)
     A = L(T,p)
     T = L(h,p)
     tables to solve all the properties.
In the thermophysical models, use p and T based table to update all the properties.
This is a good version for RGDFoam, but using different table.
The error is from interpolations.


Updated 19-12-2017
