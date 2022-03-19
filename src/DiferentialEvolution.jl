module DiferentialEvolution

using Random

export DE_classic
export DE_jither
export DE_dither
export DE_JADE

include("SupportFunctions.jl")
include("PopulationFunctions.jl")
include("MutationFunctions.jl")
include("CrossoverFunctions.jl")
include("SelectionFunctions.jl")

### Classic Diferential Evolution ###
"""
    classic_DE(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
Also known as **DE/rand/1/bin** is the simplest differential evolution algorithm. It has a constant value of **F** and **Cr**.
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

"""
DE_dither(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
"""
function DE_dither(fob::Function,generations::Integer,individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray,F::Real,Cr::Real)
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
    mutated = mutation_dither(init_population,lower_limits,upper_limits,F)
    crossed = crossover_classic(init_population,mutated,Cr)
    selected = selection_classic(fob,crossed,init_population)
    #Futher generations
    for generation in 2:generations
        init_population = selected
        mutated = mutation_dither(init_population,lower_limits,upper_limits,F)
        crossed = crossover_classic(init_population,mutated,Cr)
        selected = selection_classic(fob,crossed,init_population)
    end

    return best_individue(fob,selected)
end 



function DE_JADE(fob::Function, generations::Integer, individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray, p::Real, c::Real, μ_F=0.9, μ_Cr=0.5)

    init_population = initial_population(individues,lower_limits,upper_limits)
    number_of_parameters = length(lower_limits)

    #Initializating archive
    A = Array{Real}(undef, 0, number_of_parameters)
    
    for generation in 1:generations
        S_F = Vector{Real}(undef, 0)
        S_Cr = Vector{Real}(undef, 0)
        p_best = pbest(fob, init_population, p)

        for individue in 1:individues

            original_individue = init_population[individue,:]
            
            # Cr and F
            Cr = rand_normal_trunc(μ_Cr,0.1)
            F = rand_cauchy_trunc(μ_F,0.1)

            # Index of random individue (rnd_best)
            rnd_best = rand(1:size(p_best,1))

            # Union of population and archive
            PuA = vcat(init_population, A)
            PuA_size = size(PuA,1)

            # Random individues indexes of init_population
            rnd1 = rand_exclusive(deleteat!(collect(1:individues),individue),1)[1]
            rnd2 = rand_exclusive(deleteat!(collect(1:PuA_size),sort([individue,rnd1])),1)[1]

            # MUTATION
            mut_individue = mutation_current2pbest(original_individue, p_best[rnd_best,:], init_population[rnd1,:], PuA[rnd2,:], F)
            # CROSSOVER
            crossed_individue = crossover_JADE(original_individue, mut_individue, Cr)
            # SELECTION
            if fob(mut_individue...) <= fob(original_individue...)
                init_population[individue,:] = mut_individue
                A = vcat(A,reshape(original_individue,1,:))
                push!(S_F, F)
                push!(S_Cr, Cr)
            end
        end

        # Randomly removing solutions from A so that size(A,1) ≤ individues
        A_old_size = size(A,1)
        if A_old_size > individues
            A_new_size = rand(1:individues)
            A = A[rand_exclusive(collect(1:A_old_size), A_new_size ),:]
        end
        # Updating μ_F and μ_Cr
        μ_Cr = (1-c)*μ_Cr + c*(sum(S_Cr)/length(S_Cr))
        
        if isempty(S_F)
            continue
        end

        μ_F = (1-c)*μ_F + c*Lehmer_mean(S_F,2)
    end

    return best_individue(fob,init_population)
end



end
