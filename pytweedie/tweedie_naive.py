from math import sin, pi, exp
from scipy.special import gamma
import numpy as np
from scipy import integrate

def x(theta, alpha):
    return (alpha - 1) / alpha * (theta / (alpha - 1)) ** alpha

def N(z, alpha, n):
    terms = [None] * n

    for k in range(1, n+1):
        terms[k-1] = gamma(1 + alpha * k) / (gamma(1 + k)) * x(-1/z, alpha) ** k * sin(- k * pi * alpha)

    return sum(terms) / (pi * z)


# z, alpha, theta, n
def pdfz_tweedie(z, theta, alpha, n=1000):
    assert alpha >= 0.01 and alpha <= 0.99

    return N(z, alpha, n) * exp(z * theta - x(theta, alpha))

def pdfx_tweedie(x, theta, alpha, lambdaa, n):
    return (1.0/lambdaa) * pdfz_tweedie(x/lambdaa, theta, alpha, n)


if __name__ == '__main__':
    #for alpha in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    #    f = lambda x: pdfz_tweedie_naive(x, -1/2, alpha, n=100)

    #    print(f'alpha: {alpha}, I: {integrate.quad(f, 1e-6, 50, points=1000)}')

    print("This is a module.  Do not run it directly.")
    exit(1)
