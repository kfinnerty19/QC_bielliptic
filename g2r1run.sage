load("rank1eqs.sage")
load("qc_g2_r1_bielliptic.sage")

import sys

def Dtest(E1,E2,D,p,r1,r2):
	try:
		r1d = silent_wrapper(E1.quadratic_twist(-D))
		r2d = silent_wrapper(E2.quadratic_twist(-D))
	except Exception:
		return False
	if r1+r1d != 1 or r2+r2d != 1:
		return False
	S.<a> = NumberField(x^2+D)
	I = S.ideal(p)
	F = I.factor()
	if len(F)==2:
		return True
	else:
		return False
		
from contextlib import redirect_stdout
from io import StringIO
import multiprocessing

def targetfunc(E,queue):
    with redirect_stdout(StringIO()): 
        try:
            queue.put(E.rank()) 
        except Exception as e:
            queue.put(-100) #Dummy number to ensure failure
    

def silent_wrapper(E):
    queue = multiprocessing.Queue()
    process = multiprocessing.Process(target=targetfunc, args=(E,queue))
    process.start()
    process.join(10)
    
    if process.is_alive():       
        process.terminate()
        process.join()
        return -100 
    
    if not queue.empty():
        return queue.get()
    return -100

i = 1
for f in equations:
	print(50*"*")
	print("Curve #%s"%(i))
	print("testing for y^2 = %s"%f)
	a6 = f[6]
	a4 = f[4]
	a2 = f[2]
	a0 = f[0]
	E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6^2])
	r1 = silent_wrapper(E1)
	E2 = EllipticCurve([0, a2, 0, a0*a4, a0^2*a6])
	r2 = silent_wrapper(E2)
	goodp = []
	for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
		check = E1.has_good_reduction(p) and E2.has_good_reduction(p) and E1.is_ordinary(p) and E2.is_ordinary(p)
		if check:
			goodp.append(p)
	print("good ps are %s"%goodp)
	pairs = []
	for p in goodp:
		for D in [1..20,43,67,163]:
			check = Dtest(E1,E2,D,p,r1,r2)
			if check:
				S.<a>=NumberField(x^2+D)
				print("[%s,%s]"%(D,p))
				try:
				    rat_pts,other_pts = quadratic_chabauty_bielliptic(f,p,25,F=S)
				    pairs.append([[D,p]])
				except Exception as e:
				    #print(e)
				    pass
		if len(pairs)>=5:
		    break
	i = i+1
