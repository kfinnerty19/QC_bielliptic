# QC_bielliptic

This repository has a mildly edited version of Francesca Bianchi's SageMath code, originally written for genus 2 bielliptic curves over Q whose Jacobians have Mordell--Weil rank 2. Bianchi's code is discussed in the paper 

[BP22]: "Rational points on rank 2 genus 2 bielliptic curves in the LMFDB" by F. Bianchi and O. Padurariu. 

The first set of small edits here are to handle the case of X_0(37) over Q(i), as J_0(37)(Q(i)) has rank 2, for the paper

[BM23]: "Ogg’s Torsion conjecture: Fifty years later" by J. S. Balakrishnan and B. Mazur, with an appendix by N. Dogra. 

The second set of small edits are to generalize further to the case of curves whose Jacobians have Mordell--Weil rank 1 and check for algebraic points, for the paper

[Fin25]: "Quadratic Chabauty Experiments on Genus 2 Bielliptic Modular Curves in the LMFDB" by K. Finnerty.

The main files ares `qc_g2_bielliptic.sage` which does quadratic Chabauty on any genus 2 bielliptic curve given by an equation of the form `y^2 = a_6*x^6 + a_4*x^4 + a_2*x^2 + a_0`, with `a_i\in \ZZ`, such that the corresponding elliptic curves each have rank 1, and `qc_g2_r1_bielliptic.sage`, which allows for one of the corresponding elliptic curves to have rank 0 provided that a suitable number field is provided. The prime should be of good ordinary reduction.
This code is based on an earlier version by F. Bianchi, available at: https://github.com/bianchifrancesca/quadratic_chabauty/blob/master/quadratic_chabauty_bielliptic.sage.
The version of code by J. S. Balakrishnan is available at: https://github.com/jbalakrishnan/QC_bielliptic/tree/main

The other files provide a list of potential equations as well as files to run the analysis over these equations. In the rank 1 case, suitable number fields are computed such that the Mordell--Weil rank of the Jacobian of the curve increases to 2.

NB: For the prime 3, the code relies on the following repository: https://github.com/jbalakrishnan/AWS
