# Simple 2d unstructured CFD solver #
### Solves the Euler equations (inviscid + compressible) governing fluid flow around arbitrary 2d object using an unstructured grid comprised of triangles.  Uses a low order finite volume method with scalar diffusion for shock capturing. ###


### **Source includes:** ###
  1. Simple FV solver
  1. Mesh generation routines for Matlab (including setting boundary conditions)
  1. Sample meshes and input files

### **System Requirement:** ###
  1. Fortran 90 compiler
  1. gmake (makefile included)
  1. Matlab for pre-processing of your mesh (or use default meshes).
  1. Visualization software-write out Tecplot format data for either _Tecplot_ (commercial) or _Visit_ (freely available at https://wci.llnl.gov/codes/visit/home.html)

![http://cybo.googlecode.com/svn/trunk/doc/figs/montage.jpg](http://cybo.googlecode.com/svn/trunk/doc/figs/montage.jpg)