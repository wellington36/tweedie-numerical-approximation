import tweedie
from scipy.special import roots_legendre

def gauss_legendre_integration(f, a, b, n_points):
    x, w = roots_legendre(n_points)
    integral = 0.0
    
    for i in range(n_points):
        integral += w[i] * f(0.5 * (b - a) * x[i] + 0.5 * (b + a))
    integral *= 0.5 * (b - a)
    
    return integral

def alpha_to_p(alpha):
    return (alpha - 2) / (alpha - 1)

def alpha_to_points(alpha):
    if alpha >= 0.99:
        return 10_000
    return 1_000

for alpha in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    f = lambda x : tweedie.tweedie(mu=1, phi=1, p=alpha_to_p(alpha)).pdf(x)

    if alpha == 0.99:
        value = gauss_legendre_integration(f, 1e-6, 20, n_points=alpha_to_points(alpha))

        print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')
        continue

    value = gauss_legendre_integration(f, 1e-6, 50, n_points=alpha_to_points(alpha))

    print(f'alpha: {alpha}, I: {(value, abs(value - 1))}')