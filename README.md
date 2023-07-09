# Tweedie numerical approximation

## Tweedie's pdf implementations
- `tweedie_naive`: Naive implementation of the tweedie's pdf. Receiving $z$, $\alpha$, $\theta$ and n (number of terms to sum in N function), and applys in the following equations:

$$
  f_Z(z; \theta, \alpha) = N(z; \alpha) \times exp[z \theta - x(\theta, \alpha)],\ z \geq 0,
$$

with

$$
  N(z; \alpha) = \frac{1}{\pi z} \times \sum_{k = 1}^\infty \frac{\Gamma(1 + \alpha k)}{\Gamma(1 + k)} (x(-1/z, \alpha))^k \sin(-k \pi \alpha),
$$

and

$$
  x(\theta, \alpha) = \frac{\alpha - 1}{\alpha} \left( \frac{\theta}{\alpha - 1} \right)^\alpha
$$

In general form $f_X(x) = \frac{1}{\sigma} f_Z(x/\sigma)$. Testing $0 < \alpha < 1$, with Gaussian quadrature to integrate Tweedie densities:

|alpha | points | I                   | error           |
|---------|-----------|-----------------|-----------------|
|0.01   | 1000     |1.0000e+00| 3.3409e-08|
|0.1     | 1000     |1.0000e+00| 1.0938e-08|
|0.2     | 1000     |1.0000e+00| 1.5166e-09|
|0.3     | 1000     |1.0000e+00| 2.6935e-09|
|0.4     | 1000     |1.0000e+00| 1.6544e-10|
|0.5     | 1000     |1.0000e+00| 3.3830e-09|
|0.6     | 1000     |1.0000e+00| 2.5255e-12|
|0.7     | 1000     |1.0000e+00| 9.3295e-14|
|0.8     | 1000     |1.0000e+00| 6.0013e-09|
|0.9     | 1000     |1.0000e+00| 8.4247e-11|
|0.99   | 10000     |1.0000e+00| 1.2097e-08|


- `tweedie_dias`: Implementation of tweedie's pdf in [1]. Following the specifications in the paper:

|alpha | points | I                   | error           |
|---------|-----------|-----------------|-----------------|
|0.01   | 1000     |1.0000e+00| 3.3409e-08|
|0.1     | 1000     |1.0000e+00| 1.0938e-08|
|0.2     | 1000     |1.0000e+00| 1.5166e-09|
|0.3     | 1000     |1.0000e+00| 2.6935e-09|
|0.4     | 1000     |1.0000e+00| 1.6544e-10|
|0.5     | 1000     |1.0000e+00| 3.3830e-09|
|0.6     | 1000     |1.0000e+00| 2.5255e-12|
|0.7     | 1000     |1.0000e+00| 9.3295e-14|
|0.8     | 1000     |1.0000e+00| 6.0013e-09|
|0.9     | 1000     |1.0000e+00| 8.4247e-11|
|0.99   | 10000     |1.0000e+00| 1.2097e-08|

  # Referencies
  [1]: DIAS NL, RIBEIRO JR PJ. Practical rules for summing the series of the Tweedie probability density function with high-precision arithmetic. An Acad Bras Ciênc [Internet]. 2019;91(4):e20180268. Available from: https://doi.org/10.1590/0001-3765201920180268
