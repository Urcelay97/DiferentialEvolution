"""
    crossover_classic(population::AbstractMatrix,mutated::AbstractMatrix,Cr::Real)

Perform the crossover operation between a **`population`** matrix and **`mutated`** matrix.
The crossover probability, **`Cr`** ∈ [0,1], is a user-defined value that controls
the fraction of parameter values that are copied from the mutant.


To determine which source contributes a given parameter, uniform crossover
compares **`Cr`** to the output of a uniform random number generator, **`rand()`**.
If the random number is less than or equal to **`Cr`**, the trial parameter
is inherited from the **`mutated`** matrix; otherwise, the parameter is copied
from the **`population`** matrix. In addition, the trial parameter with randomly
chosen index, **`j`**, is taken from the mutant to ensure that the trial vector
does not duplicate. 

(*Price, K., Storn, R. M., & Lampinen, J. A. (2006).
Differential evolution: a practical approach to global optimization.
Springer Science & Business Media.*)

# Parameters

## **`population`**:

Population matrix, each row is an individue; each individue has some parameters.

## **`mutated`**:

The mutated population matrix.

## **`Cr`**:

Crossover probability.

# Examples

```julia-repl
julia> A = rand(5,3)
5×3 Matrix{Float64}:
 0.896444   0.493249   0.51672 
 0.705792   0.0381215  0.48208 
 0.298154   0.491878   0.658559
 0.635572   0.342714   0.443062
 0.0671616  0.553378   0.903124

julia> B = rand(5,3)
5×3 Matrix{Float64}:
 0.946939  0.684977   0.848902
 0.390751  0.0296715  0.789489
 0.772502  0.351496   0.889189
 0.187194  0.451319   0.943257
 0.264035  0.635716   0.477502

julia> crossover_classic(A,B,0.5)
5×3 Matrix{Float64}:
 0.896444   0.684977   0.51672
 0.390751   0.0296715  0.789489
 0.298154   0.351496   0.889189
 0.187194   0.451319   0.943257
 0.0671616  0.635716   0.477502
```
"""
function crossover_classic(population::AbstractMatrix,mutated::AbstractMatrix,Cr::Real)
#Error check
    if size(population) != size(mutated)
        throw(ErrorException("Size of population matrix $(size(population)) must be the same as the mutated matrix $(size(mutated))"))
    end

    if (Cr < 0) | (Cr > 1)
        throw(ErrorException("Cr must be in the interval [0,1]"))
    end

    crossed = similar(population)
    population_size = size(population,1)
    number_of_parameters = size(population,2)
    
    @simd for individue in 1:population_size
        @simd for parameter in 1:number_of_parameters
            if (rand() <= Cr) | (rand(1:number_of_parameters) == parameter)
                @inbounds crossed[individue,parameter] = mutated[individue,parameter]
            else
                @inbounds crossed[individue,parameter] = population[individue,parameter]
            end
        end
    end

    return crossed
end
