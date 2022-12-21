r"""
For (a sublist of) the polynomials in `allcurves` (`allcurves.sage`), this is code to do the following: 
- do quadratic Chabauty for 3 good ordinary primes (which primes is specified in ``aux-data-allcurves.sage``);
- check that all the known rational points in ``aux-data-allcurves.sage`` are recovered in each quadratic Chabauty computation;
- for all the extra p-adic points, compute `p`-adic coefficients with respect to a Mordell--Weil basis 
  (which is specified by `projE1E2` in ``aux-data-allcurves.sage``), and base-point the first listed known rational point (see `pts`).
For details, see [BP22, 4.2].

REFERENCES:
- [BP22] \F. Bianchi and O. Padurariu, "Rational points on rank 2 genus 2 bielliptic curves in the LMFDB", 2022.
"""

load("./qc_g2_bielliptic.sage")
load("./allcurves.sage")
load("./aux-data-allcurves.sage")

R = [0..5] # to consider allcurves[i] for i in R
n = 25 # quadratic Chabauty precision
N = 4  # desired p-adic precision for Mordell--Weil sieve.

def make_affine(P,X):
    x,y,z = P
    if z == 0:
        if x*y > 0:
            return "inf+"
        else:
            return "inf-"
    else:
        return X(x/z, y/z^3,1)

t1 = cputime()
coeffs_fake_lists = []
for i in R:
    print("Curve", i)
    sys.stdout.flush()
    f = allcurves[i]
    X = HyperellipticCurve(f)
    primes = ordinary_pr[i]
    allpoints = [make_affine(P,X) for P in pts[i]]
    if len(allpoints) == 0: # curves for which we don't know any rational points.
        # we considered these curves separately (see [BP22, 4.4]).
        coeffs_fake_lists.append([[[],[],[]]])
        continue
    im_div = projE1E2[i]
    base_pt = allpoints[0]
    rat_pts_lists = []
    coeffs_fake_lists_c = [[] for k in range(len(primes))]
    for p in primes:
        if primes.index(p) == 0:
            b = True
        else:
            b = False
        rat_points, other_points = quadratic_chabauty_bielliptic(f, p, n, omega_info = True, up_to_auto = b)
        if b:
            for l in range(len(rat_points)):
                rat_points_new_complete = []
                for P in rat_points[l]:
                    rat_points_new_complete.extend([P,X(P[0],-P[1], P[2]), X(-P[0],P[1],P[2]),X(-P[0],-P[1],P[2])])
                rat_points[l] = list(Set(rat_points_new_complete))
                rat_points[l].sort()
        rat_pts_lists.append(rat_points)

        print("Other points =",len(list(set().union(*other_points))))
        sys.stdout.flush()

        rat_points_new = list(set().union(*rat_points))
        rat_points_new_new = rat_points_new[:]
        if X(0,1,0) in rat_points_new_new:
            rat_points_new.remove(X(0,1,0))
            rat_points_new.extend(['inf+','inf-'])
        assert Set(rat_points_new) == Set(allpoints), "Not all rational points detected"
        coeffs_rat = coefficients_mod_pN_v2(f, allpoints, im_div, base_pt, p, N, k=5)
        coeffs_rat_int = [[ZZ(c[0]), ZZ(c[1])] for c in coeffs_rat] # this is essentially a sanity check
        # namely, it checks that the MW coefficients of the rational points are p-adically integral.
        for r in range(len(other_points)):
            if other_points[r] == []:
                coeffs_fake_lists_c[primes.index(p)].append([])
            else:
                coeffs_fake = coefficients_mod_pN_v2(f, other_points[r], im_div, base_pt, p, N, k=5)
                if min([min([coeffs_fake[i][0].precision_absolute(), coeffs_fake[i][1].precision_absolute()]) for i in range(len(coeffs_fake))]) < N:
                    coeffs_fake = coefficients_mod_pN_v2(f, other_points[r], im_div, base_pt, p, N, k=10)
                assert min([min([coeffs_fake[i][0].precision_absolute(), coeffs_fake[i][1].precision_absolute()]) for i in range(len(coeffs_fake))]) >= N,"Problem with precision of coefficients."
                coeffs_fake_int = [[ZZ(c[0]), ZZ(c[1])] for c in coeffs_fake if c[0].valuation(p)>= 0 and c[1].valuation(p) >= 0]
                coeffs_fake_lists_c[primes.index(p)].append(coeffs_fake_int)
    coeffs_fake_lists.append(coeffs_fake_lists_c)
    print(".......")
    sys.stdout.flush()
    
    #Check that not only were all the rational points detected, but also that 
    #they correspond to compatible elements of Omega, as we vary the prime. 
    assert (rat_pts_lists[0] == rat_pts_lists[1]) and (rat_pts_lists[1] == rat_pts_lists[2]), "BAD"

t2 = cputime()
print("This took", t2-t1)
sys.stdout.flush()
print(coeffs_fake_lists)
sys.stdout.flush()
