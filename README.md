# QC_bielliptic

This repository has a very mildly edited version of Francesca Bianchi's SageMath code, originally written for genus 2 bielliptic curves over Q whose Jacobians have Mordell--Weil rank 2. Bianchi's code is discussed in the paper 

[BP22]: "Rational points on rank 2 genus 2 bielliptic curves in the LMFDB" by F. Bianchi and O. Padurariu. 

The small edits here are to handle the case of X_0(37) over Q(i), as J_0(37)(Q(i)) has rank 2, for the paper

[BM23]: "Ogg’s Torsion conjecture: Fifty years later" by J. S. Balakrishnan and B. Mazur, with an appendix by N. Dogra. 

Note that no further edits have been made to handle other curves over other quadratic fields.

The main file is `qc_g2_bielliptic.sage` which does quadratic Chabauty on any genus 2 bielliptic curve given by an equation of the form `y^2 = a_6*x^6 + a_4*x^4 + a_2*x^2 + a_0`, with `a_i\in \ZZ`, such that the corresponding elliptic curves each have rank 1. The prime should be of good ordinary reduction.
This code is based on an earlier version by F. Bianchi, available at: https://github.com/bianchifrancesca/quadratic_chabauty/blob/master/quadratic_chabauty_bielliptic.sage

The file also contains code for turning the output of the quadratic Chabauty computation into an input for a Mordell--Weil sieve (as described in [BP22, S.4.2]).

NB: For the prime 3, the code relies on the following repository: https://github.com/jbalakrishnan/AWS
