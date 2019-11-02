project: slsqp
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/slsqp
summary: Modern Fortran Implementation of the SLSQP Optimization Method
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
display: private
display: protected
source: true
graph: true
exclude: pyplot_module.f90
exclude_dir: ./src/tests
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

Brief description
---------------

This is a modern object-oriented Fortran implementation of the SLSQP Optimization Method.  It can be used to solve nonlinear programming problems that seek to minimize a scalar performance index subject to nonlinear equality and inequality constraints as well as bounds on the variables.
