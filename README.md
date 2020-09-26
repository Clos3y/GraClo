# GraClo
A Maple package for finding the independent components of a given tensor, and corresponding equations.

Please refer to the Documentation_v01.pdf for more details.

# Installation
1. Download the GraClo.mpl file
2. Save it to your working Maple directory
3. Load it into Maple by using
  read("GraClo vN.mpl"); where N is the version number
4. Load the module with
  with(GraClo);

# Limitations
* Cannot work with ODEs or PDEs (maybe?)
* ASYM and SYM cannot permute over the same indicies (e.g., SYM([a,a],Psi[a,a,b]) would cause an error).
* There is no way to work explicitly with covariant and contravariant indices.
* It cannot work with components explicitly stated as functions as Maple seems to always try to solve for the independent variable, even when it is specified not to.
