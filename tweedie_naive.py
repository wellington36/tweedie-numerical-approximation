from math import sin, pi, exp
from scipy.special import gamma

def x(theta, alpha):
    return (alpha - 1) / alpha * (theta / (alpha - 1)) ** alpha

def N(z, alpha, n):
    terms = [None] * n

    for k in range(1, n+1):
        terms[k-1] = gamma(1 + alpha * k) / (gamma(1 + k)) * x(-1/z, alpha) ** k * sin(- k * pi * alpha)
        
    
    return sum(terms) / (pi * z)


# z, alpha, theta, n
def pdfz_tweedie_naive(z, alpha, theta, n):
    return N(z, alpha, n) * exp(z * theta - x(theta, alpha))

def pdfx_tweedie(x, theta, alpha, lambdaa, n):
    return (1.0/lambdaa) * pdfz_tweedie_naive(x/lambdaa, theta, alpha, n)


if __name__ == '__main__':
    print(pdfz_tweedie_naive(1, 1/2, -1/2, 100))
    #print(gamma(179))