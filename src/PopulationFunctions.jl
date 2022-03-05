module PopulationFunctions

"""
    initial_population(individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray)
Return a `m*n` matrix which `m` is the number of `individues` and `n` the number of parameters. Each individue (row) has `n` parameters
with upper and lower limits given by `lower_limits, upper_limits` arrays.

# Examples
```julia-repl
julia> initial_population(4,[-1,-5],[1,5])
4×3 Array{Float64,2}:
 -0.997462  -3.461379
  0.056423   4.564237
  0.167943  -4.467215
  0.695743   0.561379
```
"""
function initial_population(individues::Integer, lower_limits::AbstractArray, upper_limits::AbstractArray)
    
    #Error checks
    if length(lower_limits) != length(upper_limits)
        throw(ErrorException("The limits must have the same length: lower and upper limits of length $(length(lower_limits)) and $(length(upper_limits))."))
    end

    return SupportFunctions.row_mul( (upper_limits .- lower_limits), rand(individues,length(lower_limits)) ) .+ lower_limits'
    
end

end