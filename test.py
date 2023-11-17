from pytweedie.tweedie_dias import pdfz_tweedie as pdfz_tweedie_dias
from pytweedie.tweedie_well import pdfz_tweedie as pdfz_tweedie_well
import numpy as np
import matplotlib.pyplot as plt
from time import time
from scipy.special import roots_legendre

def gauss_legendre_integration(f, a, b, n_points):
    x, w = roots_legendre(n_points)
    integral = 0.0
    
    for i in range(n_points):
        integral += w[i] * f(0.5 * (b - a) * x[i] + 0.5 * (b + a))
    integral *= 0.5 * (b - a)
    
    return integral

def alpha_to_points(alpha):
    if alpha >= 0.99:
        return 10_000
    else:
        return 1_000

############ visualization ############
def visualization(tweedie_pdf=pdfz_tweedie_dias, alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: tweedie_pdf(x, -1/2, alpha)

        x = np.linspace(1e-6, 10, 100)
        y = [f(i) for i in x]
        plt.plot(x, y, label=f'alpha={alpha}')

    plt.legend()
    plt.show()


############ tweedie naive ############
def test_well(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: pdfz_tweedie_well(x, -1/2, alpha)

        if alpha == 0.99:
            value = gauss_legendre_integration(f, 1e-6, 20, n_points=alpha_to_points(alpha))

            print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')
            continue

        value = gauss_legendre_integration(f, 1e-6, 50, n_points=alpha_to_points(alpha))

        print(f'alpha: {alpha}, I: {(value, abs(value - 1))}')

############ tweedie dias ############
def test_dias(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: pdfz_tweedie_dias(x, -1/2, alpha)

        if alpha == 0.99:
            value = gauss_legendre_integration(f, 1e-6, 20, n_points=alpha_to_points(alpha))

            print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')
            continue

        value = gauss_legendre_integration(f, 1e-6, 50, n_points=alpha_to_points(alpha))

        print(f'alpha: {alpha}, I: {(value, abs(value - 1))}')


############ generate table ############
def generate_table(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    print("alpha  |  error dias  |  error dias-with-extrapolation  |  time dias  |  time dias-with-extrapolation  |  score  ")

    for alpha in alphas:
        f_well = lambda x: pdfz_tweedie_well(x, -1/2, alpha)
        f_dias = lambda x: pdfz_tweedie_dias(x, -1/2, alpha)
        steps = 3

        if alpha >= 0.99:
            time_well = 0
            for _ in range(steps):
                t0 = time()
                value_well = gauss_legendre_integration(f_well, 1e-6, 20, n_points=alpha_to_points(alpha))
                time_well += time() - t0

            time_well /= 3

            time_dias = 0
            for _ in range(steps):
                t0 = time()
                value_dias = gauss_legendre_integration(f_dias, 1e-6, 20, n_points=alpha_to_points(alpha))
                time_dias += time() - t0

            time_dias /= 3
        
        else:
            time_well = 0
            for _ in range(steps):
                t0 = time()
                value_well = gauss_legendre_integration(f_well, 1e-6, 50, n_points=alpha_to_points(alpha))
                time_well += time() - t0

            time_well /= 3

            time_dias = 0
            for _ in range(steps):
                t0 = time()
                value_dias = gauss_legendre_integration(f_dias, 1e-6, 50, n_points=alpha_to_points(alpha))
                time_dias += time() - t0

            time_dias /= 3
        
        error_dias = abs(value_dias - 1)
        error_well = abs(value_well - 1)

        print(f"{alpha} | {error_dias} | {error_well} | {100 * time_dias/alpha_to_points(alpha)} | {100 * time_well/alpha_to_points(alpha)} | {(error_dias - error_well)/min(error_dias, error_well)}")



############ main ############
if __name__ == '__main__':

    t = time()

    #test_well()
    #test_dias()
    #visualization()
    generate_table()

    print(f'time: {time() - t}')
