"""
    dither_mutation(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)
Given a matrix `A` with `lower_limits` and `upper_limits` arrays, crates a new mutated
matrix with a sacale factor `F`.

# Parameters
**`A::AbstractMatrix`**:

Input matrix

## **`lower_limits::AbstractArray`**:

Array that contains all the minimun values
that each individue could obtain.

## **`upper_limits::AbstractArray`**:

Array that contains all the maximum values
that each individue could obtain.

## **`F::Real`**:

The scale factor, F ∈ (0,1+), is a positive real number that controls
the rate at which the population evolves. While there is no upper limit on F,
effective values are seldom greater than 1.0.
(*Price, K., Storn, R. M., & Lampinen, J. A. (2006).
Differential evolution: a practical approach to global optimization.
Springer Science & Business Media.*)

# Examples
```julia-repl
julia> A = 5 .*rand(5,3)
5×3 Matrix{Float64}:
 2.48847   2.80436   1.94457
 4.04176   2.1029    2.66386
 4.05949   4.0124    1.02225
 0.428298  4.07771   4.83569
 2.92697   0.497032  4.71538
```
`A` is a matrix with minimun and maximum possible values equal to 0 and 5 respectively.
```julia-repl
julia> dither_mutation(A,[1,1,1],[5,5,5],1.3)
5×3 Matrix{Float64}:
 3.49753  2.99707  1.3635
 7.51826  3.87943  3.53243
 1.73524  3.0894   3.78479
 4.55726  3.40921  1.96368
 3.89918  3.86246  1.69276
```
"""
function dither_mutation(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)
    
    #Error checks
    if length(lower_limits) != length(upper_limits)
        throw(ErrorException("The limits must have the same length: lower and upper limits of length $(length(lower_limits)) and $(length(upper_limits))."))
    end

    if length(lower_limits) != size(A,2)
        throw(ErrorException("The number of colums of the matrix ($(size(A,2))) must be equal to the length of the vector ($(length(lower_limits)))" ))
    end

    if false in [lower_limits .< upper_limits]
        throw(ErrorException("Lower limits must be less than upper limits. This happens at indexes $(collect(1:length(lower_limits))[lower_limits.>upper_limits]...)"))
    end
    
    mutated = similar(A)
    population_size = size(A,1)
    number_of_parameters = size(A,2)

    @simd for individue in 1:population_size
        
        #Obtaining 3 differents individues without taking into account the current individue
        rnd_individues = rand_exclusive(deleteat!(collect(1:population_size),individue),3)       
        
        @simd for parameter in 1:number_of_parameters
            @inbounds mutated[individue,parameter] = A[rnd_individues[1],parameter] + F*(A[rnd_individues[2],parameter]-A[rnd_individues[3],parameter])
            if (mutated[individue,parameter] < lower_limits[parameter]) | (mutated[parameter] > upper_limits[parameter])
                @inbounds mutated[individue,parameter] = (upper_limits[parameter] - lower_limits[parameter])*rand() + lower_limits[parameter]
            end
        end
    end

    return mutated


end