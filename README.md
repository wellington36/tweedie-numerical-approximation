# Tweedie numerical approximation
Following [[2](#dunn2005)] and [[3](#jorgensen1987)] the Tweedie distributions are a family of probability distributions. This family of distributions is characterized by the formula:

$$
  Var(X) = \sigma E[X]^p.
$$

When X is a random variable with Tweedie's distribution (or $X \sim T_p(\theta, \sigma))$. The tweedie distribution includes some famous distributions like: the normal ($p = 0$), Poisson ($p = 1$), gamma ($p = 2$) and the inverse Gaussian ($p = 3$). Here we are interested in a specific interval for which it has certain calculation difficulties, since the Tweedie distribution does not have a closed formula. Following [[1](#dias)], we can write the Tweedie's pdf distribution:

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

When $0 < \alpha < 1$, when $\alpha = \frac{p - 2}{p - 1}$ .In general form $f_X(x) = \frac{1}{\sigma} f_Z(x/\sigma)$. For some values of $\alpha$ we have this graphic:

![image](https://github.com/wellington36/tweedie-numerical-approximation/assets/61877847/7844721e-f0d7-4fa4-a10a-a8b85c34f175)
*Note:* For visualization, $\sigma = 1$ and $\theta = -1/2$


## Tweedie's pdf implementations
#### `tweedie_naive`: Naive implementation of the tweedie's pdf. Receiving $z$, $\alpha$, $\theta$ and n (number of terms to sum in N function), and apply it directly to the formula above, with Gaussian quadrature to integrate Tweedie densities (with n=80, to avoid mistakes):

| alpha   | points    | I               | error ($\|I - 1\|$)|
|---------|-----------|-----------------|-----------------|
|0.01   | 1000     |$\infty$| $\infty$|
|0.1     | 1000     |1.0000e+00| 1.5666e-05|
|0.2     | 1000     |1.0000e+00| 3.1971e-06|
|0.3     | 1000     |1.0000e+00| 9.6754e-05|
|0.4     | 1000     |0.9670e+00| 3.2927e-02|
|0.5     | 1000     |1.0085e+00| 8.5399e-03|
|0.6     | 1000     |$-\infty$| $\infty$|
|0.7     | 1000     |0.0000e+00| 1.0000e-00|
|0.8     | 1000     |0.0000e+00| 1.0000e-00|
|0.9     | 1000     |0.0000e+00| 1.0000e-00|
|0.99   | 1000     |0.0000e+00| 1.0000e-00|

*(time: ~1.5 seconds)*


#### `tweedie_well`: Our tweedie's pdf implementation, starting from the naive version with modifications:

- Using mpmath with bit precision equal to 300 and evaluating 1000 series terms;
- The `acelsum` function from the `extrapolation` library with the Richardson's method;
- Use a trick from [[2]](#dunn2005), divide the terms of the series $b_k$ by $max b_k$;
- Using a trick from [[1]](#dias), return zero for series terms that can generate errors.

| alpha   | points    | I               | error ($\|I - 1\|$)|
|---------|-----------|-----------------|-----------------|
|0.01   | 1000     |0.9982e+00| 1.7359e-03|
|0.1     | 1000     |1.0000e+00| 7.6686e-11|
|0.2     | 1000     |1.0000e+00| 1.99595e-12|
|0.3     | 1000     |1.0000e+00| 6.3193e-13|
|0.4     | 1000     |$\infty$| $\infty$|
|0.5     | 1000     |$\infty$| $\infty$|
|0.6     | 1000     |1.0003e-00| 3.7573e-04|
|0.7     | 1000     |1.0139e+00| 1.3973e-02|
|0.8     | 1000     |1.03375e+00| 3.3754-02|
|0.9     | 1000     |$-\infty$| $\infty$|
|0.99   | 1000     |1.3521e-02| 0.9864e-00|

*(time: ~568 seconds)*


#### `tweedie_dias`: Implementation of tweedie's pdf in [[1](#dias)]. Following the specifications in the paper:

| alpha   | points    | I               | error ($\|I - 1\|$)|
|---------|-----------|-----------------|-----------------|
|0.01   | 1000     |1.0000e+00| 3.0163e-07|
|0.1     | 1000     |1.0000e+00| 7.6941e-11|
|0.2     | 1000     |1.0000e+00| 2.0745e-12|
|0.3     | 1000     |1.0000e+00| 3.8102e-13|
|0.4     | 1000     |1.0000e+00| 2.5213e-13|
|0.5     | 1000     |1.0000e+00| 9.1526e-13|
|0.6     | 1000     |1.0000e+00| 1.6128e-12|
|0.7     | 1000     |1.0000e+00| 8.9972e-13|
|0.8     | 1000     |1.0000e+00| 1.7730e-12|
|0.9     | 1000     |1.0000e+00| 1.2505e-12|
|0.99   | 10000     |1.0000e+00| 3.3139e-09|

*(time: ~510 seconds)*

  # Referencies
  [1]<a id="dias"></a>: DIAS NL, RIBEIRO JR PJ. Practical rules for summing the series of the Tweedie probability density function with high-precision arithmetic. An Acad Bras Ciênc [Internet]. 2019;91(4):e20180268. Available from: https://doi.org/10.1590/0001-3765201920180268
  
  [2]<a id="dunn2005"></a>: Dunn, P. K. and Smyth, G. K. (2005). Series evaluation of tweedie exponential dispersion model densities. Statistics and Computing, 15(4):267–280.
  
  [3]<a id="jorgensen1987"></a>: Jorgensen, B. (1987). Exponential dispersion models. Journal of the Royal Statistical Society. Series B (Methodological), 49(2):127–162.
