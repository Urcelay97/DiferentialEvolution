"""

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