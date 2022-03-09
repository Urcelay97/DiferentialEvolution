"""

"""
function crossover_classic(population::AbstractMatrix,mutated::AbstractMatrix,Cr::Real)
#Error check
    if size(population) != size(mutated)
        throw(ErrorException("Size of population matrix $(size(population)) must be the same as the mutated matrix $(size(mutated))"))
    end

    if (Cr < 0) | (Cr > 1)
        throw(ErrorException("Cr must be in the interval [0,1]"))
    end
######################
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
