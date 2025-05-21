load("qc_g2_bielliptic.sage")
load("rank2eqs.sage")
primes = [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
i = 1
for f in equations:
    print(50*"*")
    print("Curve # %s"%(i))
    print("testing for y^2 = %s"%f)
    for p in primes:
        try:
            print("used p = %s"%p)    
            sys.stdout.flush()
            rat_pts,other_pts = quadratic_chabauty_bielliptic(f,p,25)
        except(NotImplementedError,ArithmeticError,AssertionError):
            print("NotImplemented or bad")
            pass
    i = i+1
