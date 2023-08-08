from extrapolation import acelsum, esum
from mpmath import mp, mpf, fabs, loggamma, pi, exp

### Constants ###
prec = 100

mp.prec = prec

def x(theta, alpha):
    return (alpha - mpf(1)) / alpha * (theta / (alpha - mpf(1))) ** alpha

def b_k(alpha, B, k):
    return mp.exp(loggamma(mpf(1) + alpha * k) + k * mp.log(B) - loggamma(1 + k))

def c_k(alpha, B, k):
    return (-1)**k * b_k(alpha, B, k)

def d_k(alpha, B, k):
    return c_k(alpha, B, k) * mp.sin(-k * mp.pi * alpha)

def S_d(alpha, B, n):
    
    acel = acelsum(lambda k: d_k(alpha, B, k), transform='None', n=n, logarithm=False, precision=prec)
    
    print(f'acel: {acel[-1]}')
    
    return acel[-1]

def N(z, alpha, n):
    B = fabs(x(-1/z, alpha))

    return S_d(alpha, B, n) / (pi * z)

def pdfz_tweedie(z, theta, alpha, n=1000):
    assert alpha >= 0.01 and alpha <= 0.99
    
    z = mpf(str(z))
    theta = mpf(str(theta))
    alpha = mpf(str(alpha))

    return N(z, alpha, n) * exp(z * theta - x(theta, alpha))

def pdfx_tweedie(x, theta, alpha, lambdaa, n):
    return (1.0/lambdaa) * pdfz_tweedie(x/lambdaa, theta, alpha, n)


if __name__ == '__main__':
    print("This is a module.  Do not run it directly.")
    exit(1)