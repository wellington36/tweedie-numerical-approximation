from mpmath import loggamma, sqrt, sin, pi, fabs, exp, log, mp, mpf
from math import exp as mexp
from sys import exit
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

## nsum.py
partials = {}

def startsum(mysum):
    partials[mysum] = []

def sumthis(x, mysum):
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

def setcontrol(bool):
    global control
    control = bool

dbgoutput = False

def setdbg(bool):
    global dbgoutput
    dbgoutput = bool

def testdbg():
    if dbgoutput:
        print('True')
    else:
        print('False')
    pass

kconv = 0
ra = mpf(0)

def kcon():
    return ra

def rate():
    return ra

def kappa(theta, alpha):
    res = ((alpha - 1) / alpha) * ((theta/(alpha - 1))**alpha)
    return res

def bkshort(alpha, BB, k):
    xx = loggamma(1 + alpha*k)
    yy = k*log(BB)
    zz = loggamma(1 + k)

    bk = exp(xx - yy - zz)
    return bk

def pdfx_tweedie(x, theta, alpha, lambdaa):
    return (1.0/lambdaa) * pdfz_tweedie(x/lambdaa, theta, alpha)

def pdfz_tweedie(z, theta, alpha):
    Nz = Nz_tweedie(z, theta, alpha)
    return Nz*mexp(theta*z - kappa(theta, alpha))

def Nz_tweedie(z, theta, alpha):
    #print('Nz_tweedie:', z, theta, alpha)
    # sum of series
    mp.prec = bitprecision

    # convert to mpfs
    z = mpf(str(z))
    theta = mpf(str(theta))
    alpha = mpf(str(alpha))

    # can't do miracles
    assert mpf('0.01') < alpha < mpf('0.99')

    # calculate p
    p = (alpha - 2) / (alpha - 1)
    if dbgoutput:
        bckname = "mp-s-p%6.4fz%6.4f.dat" % (p, z)
        fbck = open(bckname, 'wt')
        fbck.write("#234567890"+"1234567890"*11)
        fbck.write("\n")
        header =  '#  k                 '+\
                  'bk                   '+\
                  'sumd                 '+\
                  'relerr'
    pass

    eps = mpf('1e-10')
    nmax = 1000
    larg = -1/z
    BB = fabs(kappa(larg, alpha))
    oB = 1/BB

    def abc(alpha):
        return exp(Psi_a + Psi_b*alpha**Psi_c)
    
    if (control) and (oB < abs(alpha)):
        return float(0)
    pass

    startsum('sumd')
    sumd = mpf(0)
    relerr = mpf(1)
    k = 1
    sign = -1
    maxter = mpf(0)
    fkmax = z**(mpf(2)-p)/(p-mpf(2))
    kmax = int(fkmax)

    # check if kmax is less than 10000
    if (kmax > nmax):
        print('too many expected iterations (> %d' % nmax)
        return float('nan')
    pass

    lbmax = (1 - alpha)*fkmax + mpf('0.5')*log(alpha)
    bmax = exp(lbmax)

    # check for debugging output
    if dbgoutput:
        fbck.write(
            ("# theta = %10.6f, alpha = %10.6f BB = %15.5e period = %d"+
             " kmax = %d bmax = %15.5e eps = %15.8e\n") % (theta, alpha, BB, p, kmax, bmax, eps)
        )
    pass

    # loop to sum the series
    while ((sumd < 0) or (fabs(relerr) > eps)):
        bk = bkshort(alpha, BB, k)/bmax
        ck = bk*sign
        dk = ck*sin(-k * pi * alpha)
        sign = -sign

        # sum without loss of precision
        sumthis(dk, 'sumd')
        sumd = sumall('sumd')
        relerr = fabs(bk/sumd)

        if dbgoutput:
            fbck.write("%4d %25.15e %25.15e %25.15e\n" % (k, bk, sumd, relerr))
        pass

        k += 1
        if (k > nmax):
            kconv = k - 1
            print(f'Nz_tweedie reached {nmax} terms without convergence (z: {z})')
            return float('nan')
        pass
    pass

    # calculate the result
    sumd *= mpf(1)/(pi * z)
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
    print(quadrature_tweedie(0.5))

    #plt.plot(np.linspace(0.1, 5, 300), [quadrature(lambda x: pdfz_tweedie(x, -1/2, 1/2), 0, 20, 1000) for x in np.linspace(0.1, 5, 300)])
    #plt.show()

    #print(quadrature(lambda x: pdfz_tweedie(x, -1/2, 1/2), 1e-6, 50, 100))