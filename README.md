# QC_bielliptic

This repository has an edited version of Francesca Bianchi's SageMath code, originally written for genus 2 bielliptic curves over Q whose Jacobians have Mordell--Weil rank 2. Bianchi's code is discussed in the paper 

[BP22]: "Rational points on rank 2 genus 2 bielliptic curves in the LMFDB" by F. Bianchi and O. Padurariu. 

The original code is available at: https://github.com/bianchifrancesca/quadratic_chabauty/blob/master/quadratic_chabauty_bielliptic.sage.

The first set of small edits here are to handle the case of X_0(37) over Q(i), as J_0(37)(Q(i)) has rank 2, for the paper

[BM23]: "Ogg’s Torsion conjecture: Fifty years later" by J. S. Balakrishnan and B. Mazur, with an appendix by N. Dogra. 

The version of code by J. S. Balakrishnan is available at: https://github.com/jbalakrishnan/QC_bielliptic/tree/main

The second set of small edits are to generalize further to the case of curves whose Jacobians have Mordell--Weil rank 1 and check for algebraic points, for the paper

[Fin25]: "Quadratic Chabauty Experiments on Genus 2 Bielliptic Modular Curves in the LMFDB" by K. Finnerty.

The files qc_g2_bielliptic.sage and qc_g2_r1_bielliptic.sage are modified from the previous repositories as described below. 

We add a check on the points that do not immediately correspond to rational points to see if they appear to satisfy algebraic relations. In the event we find such a point, we print equations that describe its coordinates.

If the curve has potential good reduction at a bad prime q, then the local height contribution at that prime can be understood entirely by information about the local height at a single Q_q-point on the curve. We check in MAGMA if the curve is locally soluble at q to confirm that such a point will exist. In the event that the code fails to obtain a Q_q point, we proceed to the general case to compute the height contribution at q. 

At each residue disc, we need to compute the function rho. For the given point P, which is a F_p-point, we let Q denote the lift to the curve over Q_p. We let f_1 and f_2 denote the images of Q to E_1 and E_2 over Q_p. If, for one f_i, we have f_i multiplied by the order of E_i over F_p is the point at infinity, then 

There are some circumstances in which we must compute the Coleman integral from the point at infinity to a p-adic points on an elliptic curve. This step can be prohibitively computationally expensive for high precision. We implement a simplification to circumvent these difficulties with the point at infinity.

In the rank 1 case, we consider the height contributions of primes that divide the discriminant of the number field in use.

Description of files:

- g2r1run.sage: loads a list of equations that are rank 1 over the rationals and runs the analysis for each of these. Suitable number fields are computed such that the Mordell--Weil rank of the Jacobian increases to 2. 

- g2r2run.sage: loads a list of equations that are rank 2 over the rationals and runs the analysis for each of these.

- qc_g2_bielliptic.sage: defines functionality to perform quadratic Chabauty on any genus 2 bielliptic curve given by an equation of the form `y^2 = a_6*x^6 + a_4*x^4 + a_2*x^2 + a_0`, with `a_i\in \ZZ`, such that the corresponding elliptic curves each have rank 1 and using a prime of good reduction.

- qc_g2_r1_bielliptic.sage: modifies the rank 2 version of the file to allow for one of the corresponding elliptic curves to have rank 0 provided that a suitable number field and a prime of good reduction is provided.

- rank1eqs.sage: a list of equations that are rank 1 over the rationals.
