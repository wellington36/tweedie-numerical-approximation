from pytweedie.tweedie_dias import pdfz_tweedie as pdfz_tweedie_dias
from pytweedie.tweedie_naive import pdfz_tweedie as pdfz_tweedie_naive
from pytweedie.tweedie_well import pdfz_tweedie as pdfz_tweedie_well
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

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
def test_naive(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: pdfz_tweedie_naive(x, -1/2, alpha, n=80)    # n not to large

        value, _ = integrate.quad(f, 1e-6, 50, points=1000)

        print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')


############ tweedie naive ############
def test_well(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: pdfz_tweedie_well(x, -1/2, alpha, n=300)    # n not to large

        value, _ = integrate.quad(f, 1e-6, 50, points=1000)

        print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')

############ tweedie dias ############
def test_dias(alphas=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]):
    for alpha in alphas:
        f = lambda x: pdfz_tweedie_dias(x, -1/2, alpha)

        if alpha == 0.99:
            value, _ = integrate.quad(f, 1e-6, 20, points=10000)

            print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')
            continue

        value, _ = integrate.quad(f, 1e-6, 50, points=1000)

        print(f'alpha: {alpha}, I: {(value, abs(1 - value))}')


############ main ############
if __name__ == '__main__':

    t = time.time()

    #test_naive()
    test_well(alphas=[0.01])
    #test_dias()
    #visualization()

    print(f'time: {time.time() - t}')
