"""
Functions obtanied in:

- https://www.sfu.ca/~ssurjano/optimization.html
- Molga, M., & Smutnicki, C. (2005). Test functions for optimization needs. Test functions for optimization needs, 101, 48.

"""

# N-Dimensional


# ACKLEY FUNCTION
"""
# Description
*Dimensions*: d

The **Ackley function** is widely used for testing optimization algorithms.
In its two-dimensional form, as shown in the plot above, it is characterized
by a nearly flat outer region, and a large hole at the centre.
The function poses a risk for optimization algorithms, particularly hillclimbing algorithms,
to be trapped in one of its many local minima.

Recommended variable values are: **a = 20, b = 0.2 and c = 2π**.

# Input Domain:

The function is usually evaluated on the hypercube **xi ∈ [-32.768, 32.768],
for all i = 1, …, d,** although it may also be restricted to a smaller domain.

# Global Minimum:

**f(X\\*) = 0**, at **X\\* = (0,...,0)**

# Link:

*https://www.sfu.ca/~ssurjano/ackley.html*
"""
function Ackley(X...)
    d = length(X)
    return -20*exp(-0.2*sqrt((1/d)*sum(X.^2))) - exp((1/d)*sum(cos.(2 .*pi.*X))) + 20 + exp(1)
end

# SCHWEFEL’S FUNCTION
"""
# Description
*Dimensions*: d

**Schwefel’s function** is deceptive in that the global minimum is geometrically distant,
over the parameter space, from the next best local minima. Therefore, the search algorithms
are potentially prone to convergence in the wrong direction.

# Input Domain:

Test area is usually restricted to hyphercube **xi ∈ [−500, 500], for all i = 1, . . . , d**.

# Global Minimum:

**f(X\\*) = -418.9829d**, at **X\\* = (420.9687,...,420.9687)**
"""
Schwefel(X...) = sum(@. -X*sin(sqrt(abs(X))))

# GRIEWANK FUNCTION
"""
# Description:
*Dimensions*: d

The **Griewank function** has many widespread local minima, which are regularly distributed. The complexity is shown in the zoomed-in plots.

# Input Domain:
The function is usually evaluated on the hypercube **xi ∈ [-600, 600], for all i = 1, …, d**.

# Global Minimum:
**f(X\\*) = 0**, at **X\\* = (0,...,0)**
"""
function Griewank(X...)
    i = 1:length(X)
    A = sum(X.^2 ./4000)
    B = prod(cos.(X./sqrt.(i)))
    return A - B + 1
end

# ROSENBROCK FUNCTION
"""
*Description*:
# Dimensions: d

The **Rosenbrock function**, also referred to as the Valley or Banana function, is a popular test problem for gradient-based optimization algorithms. It is shown in the plot above in its two-dimensional form.

The function is unimodal, and the global minimum lies in a narrow, parabolic valley.
However, even though this valley is easy to find, convergence to the minimum is difficult (Picheny et al., 2012).

# Input Domain:
The function is usually evaluated on the hypercube **xi ∈ [-5, 10], for all i = 1, …, d**, although it may be restricted
to the hypercube **xi ∈ [-2.048, 2.048], for all i = 1, …, d**.

# Global Minimum:

**f(X\\*) = 0**, at **X\\* = (1,...,1)**

# Link:

https://www.sfu.ca/~ssurjano/rosen.html
"""
function Rosenbrock(X...)
    x1 = X[1:length(X)-1]
    x2 = X[2:length(X)]
    return sum(@. 100*(x2 - x1^2)^2 + (x1 - 1)^2)
end

# MICHALEWICZ FUNCTION
"""
# Description:
*Dimensions*: d

The **Michalewicz function** has d! local minima, and it is multimodal.
The parameter m defines the steepness of they valleys and ridges; a larger m leads to
a more difficult search. The recommended value of m is m = 10.

# Input Domain:
The function is usually evaluated on the hypercube **xi ∈ [0, π], for all i = 1, …, d**.

Global Minima:

at **d = 2**: **f(X\\*) = -1.8013**, at **X\\* = (2.20, 1.57)**
at **d = 5**: **f(X\\*) = -4.687658**
at **d = 10**: **f(X\\*) = -9.66015**

# Link:

*https://www.sfu.ca/~ssurjano/michal.html*
"""
function Michalewicz(X...)
    m = 10
    i = 1:length(X)
    return -sum(@. sin(X)*(sin(i*X^2 / pi)^(2*m)))
end

# 2-Dimensional


# BUKIN FUNCTION N. 6
"""
# Description:
*Dimensions*: 2

The **sixth Bukin function** has many local minima, all of which lie in a ridge.

# Input Domain:

The function is usually evaluated on the rectangle **x1 ∈ [-15, -5], x2 ∈ [-3, 3]**.

# Global Minimum:

**f(X\\*) = 0**, at **X\\* = (-10,1)**

# Link:

*https://www.sfu.ca/~ssurjano/bukin6.html*
"""
Bukin6(x,y) = 100*sqrt(abs(y - 0.01*x^2)) + 0.01*abs(x + 10)

# CROSS IN TRAY FUNCTION
"""
# Description:
*Dimensions*: 2

The **Cross-in-Tray function** has multiple global minima. It is shown here with a smaller domain in the second plot, so that its characteristic "cross" will be visible.

# Input Domain:

The function is usually evaluated on the square **xi ∈ [-10, 10], for all i = 1, 2**.

# Global Minima:

**f(X\\*) = -2.06261**, at **X\\* = (1.3491, -1.3491), (-1.3491, 1.3491), (-1.3491, -1.3491)** and **(1.3491, 1.3491)**

# Link:

*https://www.sfu.ca/~ssurjano/crossit.html*
"""
function Cross_in_tray(x,y)
    s = sqrt(x^2 + y^2)/pi
    return -0.0001*(abs(sin(x)*sin(y)*exp(abs(100-s)))+1)^0.1
end

# DROP WAVE FUNCTION
"""
# Description:
*Dimensions*: 2

The **Drop-Wave function** is multimodal and highly complex. The second plot above shows the function on a smaller input domain, to illustrate its characteristic features.

# Input Domain:

The function is usually evaluated on the square **xi ∈ [-5.12, 5.12], for all i = 1, 2**.

# Global Minimum:

**f(X\\*) = -1**, at **X\\* = (0,0)**

# Link:

*https://www.sfu.ca/~ssurjano/drop.html*
"""
function Drop_wave(X...)
    s = sum(X.^2)
    return -(1 + cos(12*sqrt(s)))/(0.5*(s) + 2)
end

# EGGHOLDER FUNCTION
"""
# Description:
*Dimensions*: 2

The **Eggholder function** is a difficult function to optimize, because of the large number of local minima.

# Input Domain:
The function is usually evaluated on the square **xi ∈ [-512, 512], for all i = 1, 2**.

# Global Minimum:

**f(X\\*) = -959.6407**, at **X\\* = (512,404.2319)**

# Link:

*https://www.sfu.ca/~ssurjano/egg.html*
"""
function Eggholder(x,y)
    A = (y + 47)*sin(sqrt(abs(y + x/2 + 47)))
    B = x*sin(sqrt(abs(x - (y + 47))))
    return -(A + B)
end