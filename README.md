# lammps-utils

A collection of utility scripts for easing the setup of Molecular Mechanics 
simulations on chemical systems

## Motivation
When I started working with [LAMMPS](https://lammps.sandia.gov), I noticed I
needed some way to generate parameters for new molecules which are only
marginally described in the literature. I also needed to account for the
numerous combinations of internal coordinates describing a molecular system.
I could have used utilities like
[PackMol](http://m3g.iqm.unicamp.br/packmol/home.shtml) and/or
[MolTemplate](https://moltemplate.org/), but neither actually help when you need
to create and derive force field parameters from (little more than) scratch.



