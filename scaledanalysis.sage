from collections import defaultdict
load('qc_g2_bielliptic.sage')

def rankjumpcheck(E1,E2,D):
	r"""
    Check whether or not the sum of the ranks of two elliptic curves increases from 1 to 2 when considered over `K:=\QQ[\sqrt{-D}]' instead of over `\QQ' 
    INPUT:
    - ``E1`` and ``E2`` -- Two elliptic curves over `\QQ'. One should be rank 1 and one should be rank 0, but the order doesn't matter.
    - ``D`` -- An integer indicating the number field `\QQ[\sqrt{-D}]'.
    OUTPUT:
    A boolean True/False indicating if the sum of the ranks of E1/K and E2/K is 2 or not.
    """
	S.<a> = NumberField(x^2+D)
 	E1S = E1.change_ring(S)
	E2S = E2.change_ring(S)
 	sum = Rank(E1S)+Rank(E2S)
  	if sum == 2:
   		return True
	else:
 		return False

def pfinder(E1, E2, startp, D=None,pbound=15):
	r"""
    Find a prime p greater than or equal to a given starting prime such that p is of good and ordinary reduction for two given and elliptic curves and, if D is not none, splits in a given number field. 
    INPUT:
    - ``E1`` and ``E2`` -- Two elliptic curves over `\QQ'. 
    - ``startp`` -- A prime integer or None (None passed in if a previous call to the function returned None, in which case None is returned automatically)
	- ``D`` -- If not n
 	- ``pbound`` -- An integer indicating an upper bound up to which to check primes. It is disadvantageous to work with too large primes.
    OUTPUT:
    A prime p suitable for use in the given quadratic Chabauty experiment, or None indicating no such prime below the bound exists.
    """
	p = startp
	if p = None:
		return None
	checked = False
 	if D != None:
		while checked == False and p < pbound:
 			if kronecker(-D,p)==1 and E1.has_good_reduction(p1) and E2.has_good_reduction(p1) and E1.is_ordinary(p1) and E2.is_ordinary(p1):
  				checked = True
			else:
 				p = next_prime(p)
    else:
   		while checked == False and p < pbound:
 			if E1.has_good_reduction(p1) and E2.has_good_reduction(p1) and E1.is_ordinary(p1) and E2.is_ordinary(p1):
  				checked = True
			else:
 				p = next_prime(p)
	if checked == True:
 			return p
    else:
			return None
 
def qcanalysis(file,degree=1,field=None,up_to_auto=True):
    r"""
    Perform quadratic Chabauty analysis on a list of bielliptic curves all of the same rank. The analysis is performed for 2 different primes p for each curve to make sieving more feasible. 
    INPUT:
    - ``file` -- a text file in which each line contains a polynomial over `\QQ` of the form `a_6*x^6 + a_4*x^4 + a_2*x^2 + a_0`,
      `a_i \in \ZZ`. The bielliptic curve `y^2=f(x)' should have rank 1 or 2 over `\QQ' (1 if degree=2, 2 if degree=1), and all curves should have the same rank.
    - ``degree`` -- either 1 or 2, denoting if the analysis should be performed over `\QQ' or a degree 2 extension
    - ``field`` -- A fixed integer D indicating that if degree=2, the analysis should be attempted over `\QQ[\sqrt{-D}]'.
    - ``up_to_auto`` -- True/False (default True): if True, the points in the output are returned up to the automorphisms `(x,y) \mapsto (\pm x,\pm y)`
    OUTPUT:
    A dictionary where each key is a curve and has the following values:
    - Two primes `p1' and `p2' that are good and ordinary
    - Four lists: `rat_points_1', `rat_points_2', `other_points_1', and `other_points_2' that store the results of the quadratic Chabauty analysis for each p.
    - If degree=2 and field=None, there is also `D', indicating that the analysis for that curve was performed over `\QQ[\sqrt{-D}]'
    If degree=2 and field!=None, some values may be errors, indicating the reason why the analysis failed for the specified field.
    The results are presented in this order: [(D,)p1,rat_points_1,other_points_1,p2,rat_points_2,other_points_2]
    """
R.<x> = PolynomialRing(Rationals())
result = defaultdict(list)
for line in open('rank2forsage.txt').readlines():
	C = eval(line)
	a6 = C[6]
	a4 = C[4]
	a2 = C[2]
	a0 = C[0]
	E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6^2])
	E2 = EllipticCurve([0, a2, 0, a0*a4, a0^2*a6])
	if degree = 1:
 		assert Rank(E1)==1 and Rank(E2)==1, "Both elliptic curves need rank 1 to perform analysis over QQ"
 		p1 = pfinder(E1,E2,5)
	    	p2 = pfinder(E1,E2,next_prime(p1))
		if p2 = None:
			result[C].append(["No 2 suitable primes"])
			continue
		rat_points_1, other_points_1 = quadratic_chabauty_bielliptic(C,p1,20,up_to_auto=up_to_auto)
  		rat_points_2, other_points_2 = quadratic_chabauty_bielliptic(C,p2,20,up_to_auto=up_to_auto)
		result[C].append([p1,rat_points_1,other_points_1,p2,rat_points_2,other_points_2])
	else:
 		assert degree=2, "Only degree 1 and 2 supported"
   		if field = None:
			#find a D value that works (checking up to 20), ie over which the ranks jump
   			D = 1
			p2 = None
			while D<20 and p2==None:
	  			while rankjumpcheck(E1,E2,D)==False and D<20: 
   					D = D+1
	   			if rankjumpcheck(E1,E2,D)==False:
					result[C].append(["no good (D,p1,p2)"])
					continue
				S.<a> = NumberField(x^2+D)
	  			p1 = pfinder(E1,E2,5,D=D)
				p2 = pfinder(E1,E2,next_prime(p1),D=D)
				if p2==None:
					D=D+1
   			rat_points_1, other_points_1 = quadratic_chabauty_bielliptic(C,p1,20,up_to_auto=up_to_auto,F=S)
  			rat_points_2, other_points_2 = quadratic_chabauty_bielliptic(C,p2,20,up_to_auto=up_to_auto,F=S)
			result[C].append([D,p1,rat_points_1,other_points_1,p2,rat_points_2,other_points_2])
  		else:
			D = field
   			S.<a> = NumberField(x^2+D)
	  		p1 = pfinder(E1,E2,5,D=D)
			p2 = pfinder(E1,E2,next_prime(p1),D=D)
			if p2 == None:
   				result[C].append(["No 2 suitable primes"])
	   		else:
   				rat_points_1, other_points_1 = quadratic_chabauty_bielliptic(C,p1,20,up_to_auto=up_to_auto,F=S)
  				rat_points_2, other_points_2 = quadratic_chabauty_bielliptic(C,p2,20,up_to_auto=up_to_auto,F=S)
				result[C].append([p1,rat_points_1,other_points_1,p2,rat_points_2,other_points_2])

   return result
   
   			
