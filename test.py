import numpy as np
from scipy import integrate
from scipy.stats import norm
import tweedie
from glum import TweedieDistribution

#print(TweedieDistribution(1.9).log_likelihood(2, 1, 1, 1))

#print(np.random.normal())

#print(integrate.quad(lambda x: TweedieDistribution(1.9).log_likelihood(x, 1, 1, 1), 0, 10))

def alpha_to_p(alpha):
    return (alpha - 2) / (alpha - 1)

#print(tweedie.tweedie(mu=1, phi=1.9, p=alpha_to_p(0.5)).pdf(1))



#print(integrate.quad(lambda x: tweedie.tweedie(mu=1, phi=1, p=alpha_to_p(0.5)).pdf(x), 1e-6, 50, points=1000))


for i in [0.005, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.995]:
    print(i, tweedie.tweedie(mu=1, phi=1.9, p=alpha_to_p(i)).pdf(1))

    print(integrate.quad(lambda x: tweedie.tweedie(mu=1, phi=1, p=alpha_to_p(0.5)).pdf(x), 1e-6, 50, points=1000))