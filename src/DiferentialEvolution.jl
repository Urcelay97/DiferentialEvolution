module DiferentialEvolution

export DE_classic
export DE_jither

include("SupportFunctions.jl")
include("PopulationFunctions.jl")
include("MutationFunctions.jl")
include("CrossoverFunctions.jl")
include("SelectionFunctions.jl")

### Classic Diferential Evolution ###
"""
    classic_DE(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)

"""
function DE_classic(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
    #Error checks

    if length(lower_limits) != length(upper_limits)
        throw(ErrorException("The limits must have the same length: lower and upper limits of length $(length(lower_limits)) and $(length(upper_limits))."))
    end

    if false in [lower_limits .< upper_limits]
        throw(ErrorException("Lower limits must be less than upper limits. This happens at indexes $(collect(1:length(lower_limits))[lower_limits.>upper_limits]...)"))
    end

    if F == 0
        throw(ErrorException("F cannot be zero"))
    end

    if (Cr < 0) | (Cr > 1)
        throw(ErrorException("Cr must be in the interval [0,1]"))
    end

    #DE Algorithm
    
    #First generation
    init_population = initial_population(individues,lower_limits,upper_limits)
    mutated = mutation_classic(init_population,lower_limits,upper_limits,F)
    crossed = crossover_classic(init_population,mutated,Cr)
    selected = selection_classic(fob,crossed,init_population)
    #Futher generations
    for generation in 2:generations
        init_population = selected
        mutated = mutation_classic(init_population,lower_limits,upper_limits,F)
        crossed = crossover_classic(init_population,mutated,Cr)
        selected = selection_classic(fob,crossed,init_population)
    end

    return best_individue(fob,selected)
end  


### Jither Diferential Evolution ###
"""
DE_jither(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
"""
function DE_jither(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
    #Error checks

    if length(lower_limits) != length(upper_limits)
        throw(ErrorException("The limits must have the same length: lower and upper limits of length $(length(lower_limits)) and $(length(upper_limits))."))
    end

    if false in [lower_limits .< upper_limits]
        throw(ErrorException("Lower limits must be less than upper limits. This happens at indexes $(collect(1:length(lower_limits))[lower_limits.>upper_limits]...)"))
    end

    if F == 0
        throw(ErrorException("F cannot be zero"))
    end

    if (Cr < 0) | (Cr > 1)
        throw(ErrorException("Cr must be in the interval [0,1]"))
    end

    #DE Algorithm
    
    #First generation
    init_population = initial_population(individues,lower_limits,upper_limits)
    mutated = mutation_jither(init_population,lower_limits,upper_limits,F)
    crossed = crossover_classic(init_population,mutated,Cr)
    selected = selection_classic(fob,crossed,init_population)
    #Futher generations
    for generation in 2:generations
        init_population = selected
        mutated = mutation_jither(init_population,lower_limits,upper_limits,F)
        crossed = crossover_classic(init_population,mutated,Cr)
        selected = selection_classic(fob,crossed,init_population)
    end

    return best_individue(fob,selected)
end  
end
