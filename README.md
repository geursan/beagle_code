# beagle_code

Hello world! This is beagle_code version 1, developed for pedagogical purposes to introduce students to astrophysical fluid dynamics. It uses the HLL method to solve the Euler equations in 1D and 2D, and is written in Fortran 90. For the 2D case, it is possible to generate a gnuplot file using a Fortran program included in the directory. The initial conditions are provided in the Fortran files.

To compile 1d or 2d version, just do the standard Fortran way: gfortran euler**.F90 -o beagle

To run: ./beagle

You can modify the code and use the compiler Fortran options as your convenience. 

For doubts or questions, feel free to contact Gerardo Urrutia: geursan[at]gmail.com / gurrutia[at]cft.edu.pl

Recommendations:

CFD methods:
For a detailed explanation of CFD methods, refer to Eleuterio F. Toro - "Riemann Solvers and Numerical Methods for Fluid Dynamics."

Astrophysical applications:
- A classical CFD application in the dynamics of the interstellar medium by Alex Raga and his group can be found in Chapter III of their book Astroplasmas here: https://bigbang.nucleares.unam.mx/astroplasmas/images/stories/pdf/libro-mi.pdf.

- For special relativistic fluid dynamics, refer to Fabio de Colle's paper here: https://iopscience.iop.org/article/10.1088/0004-637X/746/2/122/meta.

- For general relativistic cases (e.g., FM torus) including tabulated equations of state, see Agnieszka Janiuk's prescription here:https://github.com/agnieszkajaniuk/HARM_COOL.

