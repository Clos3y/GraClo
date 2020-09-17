# GraClo
A Maple package for finding the independent components of a given tensor, and corresponding equations.

Please refer to the Documentation.pdf for more details.

# Installation
1. Download the GraClo.mpl file
2. Save it to your working Maple directory
3. Load it into Maple by using
  read("GraClo.mpl");
4. Load the module with
  with(GraClo);

# Limitations
* Cannot work with ODEs or PDEs
* The code determining the oscillating sum for the ASYM command is 'unstable'. It relies on the way Maple orders terms. If that were to be changed, then I believe it would break
* There is no way to work explicitly with covariant and contravariant indices.
* It cannot work with components explicitly stated as functions as Maple seems to always try to solve for the independent variable, even when it is specified not to.
