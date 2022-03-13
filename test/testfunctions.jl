"""
Functions obtanied in https://www.sfu.ca/~ssurjano/optimization.html

"""


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

# BUKIN FUNCTION N. 6
"""
# Description:
*Dimensions*: 2

The **sixth Bukin function** has many local minima, all of which lie in a ridge.

# Input Domain:

The function is usually evaluated on the rectangle **x1 ∈ [-15, -5], x2 ∈ [-3, 3]**.

# Global Minimum:

**f(X\\*) = 0**, at **X\\* = (-10,1)**

# Link

*https://www.sfu.ca/~ssurjano/bukin6.html*
"""
Bukin6(x,y) = 100*sqrt(abs(y - 0.01*x^2)) + 0.01*abs(x + 10)

# CROSS IN TRAY FUNCTION
function Cross_in_tray(x,y)
    s = sqrt(x^2 + y^2)/pi
    return -0.0001*(abs(sin(x)*sin(y)*exp(abs(100-s)))+1)^0.1
end

# DROP WAVE FUNCTION
function Drop_wave(X...)
    s = sum(X.^2)
    return -(1 + cos(12*sqrt(s)))/(0.5*(s) + 2)
end