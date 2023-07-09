from mpmath import loggamma, sqrt, sin, pi, fabs,\
exp, log, mp, mpf
from math import exp as mexp
from sys import exit
import numpy as np
from scipy import integrate

partials = {}

def startsum(mysum):
    partials[mysum] = []

def sumthis(x,mysum):
    i = 0

    for y in partials[mysum]:
        if abs(x) < abs(y):
            x, y = y, x
        pass

        hi = x + y
        lo = y - (hi - x)

        if lo:
            partials[mysum][i] = lo
            i += 1
        pass

        x = hi
    
    partials[mysum][i:] = [x]

    return

def sumall(mysum):
    return sum(partials[mysum])

## Constants
Psi_a = -5.7573749548413
Psi_b = 5.74962006661501
Psi_c = 1.39742867043842

## Tweedie
bitprecision = 1000

control = True          # by default returns zero if z is too small

def setcontrol(bool): # sets/unsets control
    global control
    control = bool

dbgoutput = False # in the calling program, change this var to
# True with setdbg() if you want debugging
# output (Figure 1)

def setdbg(bool): # sets/ unsets debugging info in auxiliary
                    # files
    global dbgoutput
    dbgoutput = bool

def testdbg(): # tests the variable dbgoutput
    if dbgoutput:
        print('True')
    else:
        print('False')
    pass

kconv = 0
ra = mpf(0)

def kcon():
    return kconv

def rate():
    return ra

def kappa(theta,alpha):
    res = ((alpha-1)/alpha )*((theta/(alpha-1))**alpha)
    return res

def bkshort(alpha,BB,k):
    xx = loggamma(1 + alpha*k)
    yy = k*log(BB)
    zz = loggamma(1 + k)
    bk = exp(xx + yy - zz)
    return bk

def pdfx_tweedie(x,theta,alpha,lambdaa):
    return (1.0/lambdaa)*pdfz_tweedie(x/lambdaa,theta,alpha)

def pdfz_tweedie(z,theta,alpha):
    Nz = Nz_tweedie(z,theta,alpha)
    return Nz*mexp(theta*z - kappa(theta,alpha))

def Nz_tweedie(z,theta,alpha):
    # --------------------------------------------------------------------
    # sum the series
    # --------------------------------------------------------------------
    # set the precision for mpf
    # --------------------------------------------------------------------
    mp.prec = bitprecision
    # --------------------------------------------------------------------
    # convert from floats to mpfs with exactly the same digits: hence the
    # need for str()
    # --------------------------------------------------------------------
    z = mpf(str(z))
    theta = mpf(str(theta))
    alpha = mpf(str(alpha))
    # --------------------------------------------------------------------
    # I can't do miracles
    # --------------------------------------------------------------------
    #assert mpf('0.01') < alpha < mpf('0.99')
    # --------------------------------------------------------------------
    # calculate p initially as a float
    # --------------------------------------------------------------------
    p = (alpha-2)/(alpha-1)
    if dbgoutput:
        bckname = "mp-s-p%6.4fz%6.4f.dat" % (p,z)
        fbck = open(bckname,'wt')
        fbck.write("#234567890"+"1234567890"*11)
        fbck.write("\n")
        header = '# k '+\
                    'bk '+\
                    'sumd '+\
                    'relerr'
    pass
    eps = mpf('1.0e-10') # a relative error
    nmax = 10000 # maximum number of terms to sum
    larg = -1/z # auxiliary
    BB = fabs(kappa(larg,alpha)) # it is possible that BB < 0
    oB = 1/BB # the reciprocal

    def abc(alpha): # is it possible?
        return exp(Psi_a + Psi_b*alpha**Psi_c)

    if (control) and (oB < abc(alpha)): # returns zero
        return(float(0.0))
    pass

    startsum('sumd') # the sum of dk
    sumd = mpf(0) # initialize
    relerr = mpf(1) # a relative error
    k = 1 # k
    sign = -1 # the sign of (-1)**k
    maxter = mpf(0) # the maximum term so far
    fkmax = z**(mpf(2)-p)/(p-mpf(2))# the estimated locus of the max
    kmax = int(fkmax)
# --------------------------------------------------------------------
# for double precision, check if kmax is less than 10000
# --------------------------------------------------------------------
    if (kmax > nmax) :
        print('too many expected iterations (> %d)' % nmax)
        return(float('nan'))
    pass
    lbmax = (1-alpha)*fkmax + mpf('0.5')*log(alpha)# the log of the max
    bmax = exp(lbmax)
# --------------------------------------------------------------------
# check for debugging output
# --------------------------------------------------------------------
    if dbgoutput :
        fbck.write(
        ("# theta=%10.6f alpha=%10.6f BB = %15.5e period = %d"+
        " kmax = %d, bmax = %15.8e eps=%15.8e\n") %
        (theta,alpha,BB,9999,kmax,bmax,eps))
    pass
# --------------------------------------------------------------------
# loop to sum the series
# --------------------------------------------------------------------
    while ((sumd < 0) or (fabs(relerr) > eps)):
        bk = bkshort(alpha,BB,k)/bmax
        ck = bk*sign
        dk = ck*sin(-k*pi*alpha)
        sign = -sign
# --------------------------------------------------------------------
# sum without loss of precision
# --------------------------------------------------------------------
        sumthis(dk,'sumd')
        sumd = sumall('sumd')
        relerr = fabs(bk/sumd)
        if dbgoutput :
            fbck.write("%4d %25.15e %25.15e %25.15e\n" %
            (k,bk,sumd,relerr))
        pass
        k += 1
        if (k > nmax):
            kconv = k - 1
            print('Nz_tweedie reached %d terms without convergence' % nmax)
            return(float('nan'))
        pass
    pass # end of while
    sumd *= mpf(1)/(pi*z)
    sumd *= bmax
    if dbgoutput:
        fbck.close()
    pass
    kconv = k - 1
    return float(sumd)

def quadrature_tweedie(alpha):
    f = lambda x: pdfz_tweedie(x, -1/2, alpha)

    return integrate.quad(f, 1e-6, 50, points=1000)

def quadrature(f, a, b, p):
    samples = [None] * p

    for i in range(p):
        r = np.random.rand() * (b - a) + a

        samples[i] = f(r)
    
    return np.mean(samples)

if __name__ == '__main__':
    #print("This is a module.  Do not run it directly.")
    #exit(1)

    # test
    #print(pdfz_tweedie(1, -1/2, 1/2))
    #print(quadrature_tweedie(0.5))

    for i in [0.5]:
        print(quadrature_tweedie(i))