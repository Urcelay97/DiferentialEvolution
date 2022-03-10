"""
    function selection_classic(fob::Function, trial::AbstractMatrix, population::AbstractMatrix)
        
"""
function selection_classic(fob::Function, trial::AbstractMatrix, population::AbstractMatrix)
    #error checks
    if size(trial) != size(population)
        throw(ErrorException("Size of trial matrix $(size(trial)) must be the same as the population matrix $(size(population))"))
    end

    population_size = size(population,1)
    selected = similar(population)

    @simd for individue in 1:population_size
        @inbounds if fob(trial[individue,:]...) < fob(population[individue,:]...)
            @inbounds selected[individue,:] = trial[individue,:]
        else
            @inbounds selected[individue,:] = population[individue,:]
        end
    end

    return selected
end        



