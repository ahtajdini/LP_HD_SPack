# Linear Programming Bound on Sphere Packing 
This repository included Julia and Mathematica code for calculating the upper bound on the packing density of spheres in different dimensions. This is based on a paper that my collaborators and I published for the papers cited in references.

## Quick Introduction
Sphere packing is an old problem in discrete geometry with a rich history. The problem is seeking to maximize the fraction of $\mathrm{R}^{D}$ covered with non-overlapping spheres. Sphere packing has applications in communication theory, error correcting codes, toy models of granular material, modular forms, and even biology. Despite how simple is to state the problem, it is a notoriously hard problem to solve. The only dimensions that the sphere packing has been solved are $D=1,2,3,8,24$. The solution in two dimensions $D=2$ is the Hexagonal lattice and was rigorously proved in 1942 by Toth. $D=3$ is the famous Kepler conjecture. Being asked about the optimal arrangement of cannonballs, Kepler conjectured that a a family of lattices consists of fcc and hpc lattices are the optimal packing in three dimensions. The Kepler conjecture was proved almost three centuries later by Thomas Hales in 1998. The solutions in $D=8, and $24$ correspond to well-known $E_8$ and Leech lattice, respectively. $D=8$ was proved by Viazofska in 2016. A week later after Viazofska published her work,  Cohn, Kumar, Miller, Radchenko, and Viazovska proved the optimal packing in 24 dimensions.

In all other dimensions, the sphere packing has not been solved. There is not much geometric intuition about higher dimensions. However, there are various bounds known on packing density. Some bounds are based on explict constructions, some of more abstract. In general, every dimension appears to be different and the less results are known in high dimensions. An important class of bounds for all dimensions are based on linear programming (LP). In the LP method, the space of functions acting on an arbitrary hypothetical packings is considered.  Using Fourier analysis and working in dual space of these functions, one finds an upper bound on packing density. See references for more details.

## What is being computed in the code?

The mathematical formulation of this repository is the following: Given dimension $D$, solve $2N$ non-linear equations below for $2N$ variables,
$$f_k(0) + \sum_{n=1}^{N} d_n f_k(\Delta_n) = 0 \qquad \text{for}\qquad 1 \le k \le 2N$$,
where 
$$ f_k(\Delta) = L^{D/2-1}_{2*k-1}(4\pi \Delta) e^{-2 \pi \Delta}$$











