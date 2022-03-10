"""
    mutation:classic(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)
Given a matrix `A` with `lower_limits` and `upper_limits` arrays, crates a new mutated
matrix with a constant sacale factor `F`.

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
A = 5 .*rand(5,3)
5×3 Matrix{Float64}:
 1.18384  0.587246  4.46021
 4.58842  3.52444   1.39847
 1.08242  4.8591    4.9314 
 3.54709  3.17776   4.42353
 4.62508  4.84452   4.59266
```
`A` is a matrix with minimun and maximum possible values equal to 0 and 5 respectively.
```julia-repl
julia> mutation_classic(A,[0,0,0],[5,5,5],1.3)
5×3 Matrix{Float64}:        
 4.53165   4.40841   2.46325
 1.55285   1.47686   4.64034
 3.72591   0.302261  1.29996
 1.2315    2.30335   0.0760031
 0.797949  1.3387    0.738242
```
"""
function mutation_classic(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)

    mutated = similar(A)
    population_size = size(A,1)
    number_of_parameters = size(A,2)

    @simd for individue in 1:population_size
        
        #Obtaining 3 differents individues without taking into account the current individue

        rnd_individues = rand_exclusive(deleteat!(collect(1:population_size),individue),3)       
        
        @simd for parameter in 1:number_of_parameters
            @inbounds mutated[individue,parameter] = A[rnd_individues[1],parameter] + F*(A[rnd_individues[2],parameter]-A[rnd_individues[3],parameter])
            if (mutated[individue,parameter] < lower_limits[parameter]) | (mutated[individue,parameter] > upper_limits[parameter])
                @inbounds mutated[individue,parameter] = (upper_limits[parameter] - lower_limits[parameter])*rand() + lower_limits[parameter]
            end            
        end
    end

    return mutated


end

"""
    mutation_jither(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)
Similar to **`mutation_classic`**. Given a population matrix `A`, perform the mutation with a random scale factor
between (0,1]. Each scale factor is diferent for each parameter.
"""
function mutation_jither(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)

    mutated = similar(A)
    population_size = size(A,1)
    number_of_parameters = size(A,2)

    @simd for individue in 1:population_size
        
        #Obtaining 3 differents individues without taking into account the current individue

        rnd_individues = rand_exclusive(deleteat!(collect(1:population_size),individue),3)       
        
        @simd for parameter in 1:number_of_parameters
            @inbounds mutated[individue,parameter] = A[rnd_individues[1],parameter] + rand()*F*(A[rnd_individues[2],parameter]-A[rnd_individues[3],parameter])
            if (mutated[individue,parameter] < lower_limits[parameter]) | (mutated[individue,parameter] > upper_limits[parameter])
                @inbounds mutated[individue,parameter] = (upper_limits[parameter] - lower_limits[parameter])*rand() + lower_limits[parameter]
            end            
        end
    end

    return mutated
end



"""
    mutation_jither(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)
Similar to **`mutation_classic`**. Given a population matrix `A`, perforn the mutation with a random scale factor
in the interval (0,F). Each scale factor is diferent for each individue but the same for all parameters.
"""
function mutation_dither(A::AbstractMatrix, lower_limits::AbstractArray, upper_limits::AbstractArray, F::Real)

    mutated = similar(A)
    population_size = size(A,1)
    number_of_parameters = size(A,2)

    @simd for individue in 1:population_size
        
        #Obtaining 3 differents individues without taking into account the current individue

        rnd_individues = rand_exclusive(deleteat!(collect(1:population_size),individue),3)       
        rnd = rand()
        @simd for parameter in 1:number_of_parameters
            @inbounds mutated[individue,parameter] = A[rnd_individues[1],parameter] + rnd*F*(A[rnd_individues[2],parameter]-A[rnd_individues[3],parameter])
            if (mutated[individue,parameter] < lower_limits[parameter]) | (mutated[individue,parameter] > upper_limits[parameter])
                @inbounds mutated[individue,parameter] = (upper_limits[parameter] - lower_limits[parameter])*rand() + lower_limits[parameter]
            end            
        end
    end

    return mutated
end
