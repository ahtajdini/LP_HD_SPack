# Linear Programming Bound on Sphere Packing 
This repository included Julia and Mathematica code for calculating the upper bound on the packing density of spheres in different dimensions $D$. I used the similar code here to obtain the best known linear programming bound on packing density in high dimensions, $D \sim 1-2000$. The work is based on a paper that my collaborators and I published cited in references.

## Quick Introduction
Sphere packing is an old problem in discrete geometry with a rich history. The problem is seeking to maximize the fraction of $\mathrm{R}^{D}$ covered with non-overlapping spheres. Sphere packing has applications in communication theory, error correcting codes, toy models of granular material, modular forms, and even biology. Despite how simple is to state the problem, it is a notoriously hard problem to solve. The only dimensions that the sphere packing has been solved are $D=1,2,3,8,24$. The solution in two dimensions $D=2$ is the Hexagonal lattice and was rigorously proved in 1942 by Toth. $D=3$ is the famous Kepler conjecture. Being asked about the optimal arrangement of cannonballs, Kepler conjectured that a a family of lattices consists of fcc and hpc lattices are the optimal packing in three dimensions. The Kepler conjecture was proved almost three centuries later by Thomas Hales in 1998. The solutions in $D=8, and $24$ correspond to well-known $E_8$ and Leech lattice, respectively. $D=8$ was proved by Viazofska in 2016. A week later after Viazofska published her work,  Cohn, Kumar, Miller, Radchenko, and Viazovska proved the optimal packing in 24 dimensions.

In all other dimensions, the sphere packing has not been solved. There is not much geometric intuition about higher dimensions. However, there are various bounds known on packing density. Some bounds are based on explict constructions, some of more abstract. In general, every dimension appears to be different and the less results are known in high dimensions. An important class of bounds for all dimensions are based on linear programming (LP). In the LP method, the space of functions acting on an arbitrary hypothetical packings is considered.  Using Fourier analysis and working in dual space of these functions, one finds an upper bound on packing density. See references for more details.

## What is being computed in the code?

The mathematical formulation of this repository is the following: Given dimension $D$, solve $2N$ non-linear equations below for $2N$ variables,
$$f_k(0) + \sum_{n=1}^{N} d_n f_k(\Delta_n) = 0 \qquad \text{for}\qquad 1 \le k \le 2N$$,
where 
$$f_k(\Delta) = L\_{2k-1}^{D/2-1} (4\pi \Delta) e^{-2 \pi \Delta}$$,
and ${L^{D/2-1}}\_{2k-1}$ are Legendre Polynomials. The unknowns are $\Delta_n, d_{\Delta_n}$ for $n=1,2,\cdots, N$. Solving these equations for $\Delta_1$ gives an upper bound on the any packing density $\rho$,
$$\rho \le \frac{\pi^{D/2}}{(D/2)!} \left(\frac{r}{2}\right)^{D},\qquad r = \sqrt{2\Delta_1}$$.
We choose the variable $N$. The upper bound on density improves as $N$ increases, however we need to solve more non-linear equations.

## A few considerations for using the code

The code computes packing density by solving the above-mentioned non-linear equations. For any $D$, the goal is to solve $\Delta_1$ for higher $N$.  I solve them by Newton's method. Newton's method is fast and stable, however it requires matrix inversion and every element of matrix has to be stored to high precision. I have tried Gradient descent (without momentum of Adam) in the past, but the convergence was very slow.

Solving equations is extremely sensitive to initial guess for the values of $\Delta_i$. In particular, if we sort $\Delta_i$, I had to guess the first few $\Delta_i$ with $\mathcal(O)(10^{-6} - 10^{-2})$ accuracy of the result, otherwise the equations do not converge. In practice, for a given $D$, the bounds are converged once $N\sim D$. Due to the non-linear nature of equations, they have to be solved with high precision computing. For $D=2000$, I used $~1000$ digit precision. I used originally Mathematica file for this problem which can do arbitrary-precision arithmetic. I did the calculation using Python's mpmath library. However, I found matrix inversion in mpmath is much slower than Mathematica and Julia. Since Mathematica is not open- source, I share the Julia notebook here. 











