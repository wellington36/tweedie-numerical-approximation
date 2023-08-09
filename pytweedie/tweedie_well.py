from extrapolation import acelsum, esum
from mpmath import mp, mpf, fabs, loggamma, pi, exp, log

### Constants ###
prec = 300
Psi_a = -5.7573749548413
Psi_b = 5.74962006661501
Psi_c = 1.39742867043842

### Tweedie ###
mp.prec = prec

def x(theta, alpha):
    return (alpha - mpf(1)) / alpha * (theta / (alpha - mpf(1))) ** alpha

def b_k(alpha, B, bmax, k):
    return exp(loggamma(mpf(1) + alpha * k) + k * log(B) - loggamma(1 + k)) / bmax

def c_k(alpha, B, bmax, k):
    return (-1)**k * b_k(alpha, B, bmax, k)

def d_k(alpha, B, bmax, k):
    return c_k(alpha, B, bmax, k) * mp.sin(-k * mp.pi * alpha)

def S_d(alpha, B, bmax, n):
    
    acel = acelsum(lambda k: d_k(alpha, B, bmax, k), transform='Richardson', n=n, logarithm=False, precision=prec)
    
    return acel[-1]

def N(z, alpha, n):
    # dunn estimate to find bmax
    p = (alpha-2)/(alpha-1)
    fkmax = z**(mpf(2)-p)/(p-mpf(2))
    lbmax = (mpf(1)-alpha)*fkmax + mpf('0.5')*log(alpha)
    bmax = exp(lbmax)
    
    B = fabs(x(-1/z, alpha))
    oB = 1/B # the reciprocal

    def abc(alpha): # from Dias (check later)
        return exp(Psi_a + Psi_b*alpha**Psi_c)

    if oB < abc(alpha): # returns zero
        return(float(0.0))
    pass
    
    sumd = S_d(alpha, B, bmax, n)
    
    sumd *= mpf(1)/(pi*z)
    sumd *= bmax

    return sumd

def pdfz_tweedie(z, theta, alpha, n=1000):
    assert alpha >= 0.01 and alpha <= 0.99
    
    z = mpf(str(z))
    theta = mpf(str(theta))
    alpha = mpf(str(alpha))
    
    return N(z, alpha, n) * exp(z * theta - x(theta, alpha))

def pdfx_tweedie(x, theta, alpha, lambdaa, n):
    x = mpf(str(x))
    z = mpf(str(z))
    theta = mpf(str(theta))
    alpha = mpf(str(alpha))
    lambdaa = mpf(str(lambdaa))
    
    return (mpf(1)/lambdaa) * pdfz_tweedie(x/lambdaa, theta, alpha, n)


if __name__ == '__main__':
    print("This is a module.  Do not run it directly.")
    exit(1)