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

#这个Repository是用来存储第一个版本的 RGDFoam，使用的是
     p = L(rho,e)
     h = L(rho,e)
     A = L(T,p)
     T = L(h,p)
     表格去求解所有的气体参数。
 在热力学模块，用的是p,T去查表。
 这一个求解器是一个完整的可用的求解器，虽然使用不同的表格会导致误差，但是通过加密表格可以减少误差。



Updated 19-12-2017
