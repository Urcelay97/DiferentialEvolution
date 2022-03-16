"""
    row_mul(v::AbstractArray, A::AbstractMatrix)
Multiply all the rows of a matrix `A` by some vector `v` element by element.

"""
function row_mul(v::AbstractArray, A::AbstractMatrix)

    #Error
    if length(v) != size(A,2)
        return ErrorException("The number of colums of the matrix ($(size(A,2))) must be equal to the length of the vector ($(length(v)))" )
    end

    B = similar(A)

    @simd for i in 1:size(A,1)
        @simd for j in 1:size(A,2)
            @inbounds B[i,j] = v[j]*A[i,j]
        end
    end

    return B
end

"""
    rand_exclusive(v::AbstractArray,n::Integer)
Returns an array with `n` different random elements of the array `v`.
"""
function rand_exclusive(v::AbstractArray,n::Integer)
    
    #Error checks
    if n>length(v)
        return ErrorException("Length of the array `v` ($(length(v))) must be less or equal than `n` ($n)")
    end

    a = Array{Number}(undef,n)
    a[1] = rand(v)
    
    @simd for i in 2:n
        tmp = rand(v)
        while tmp in a[1:i-1]
            tmp = rand(v)
        end
        a[i] = tmp
    end
    return a
end

function best_individue(fob::Function,A::AbstractMatrix)
    vals = findmin(@inbounds [fob(A[i,:]...) for i in 1:size(A,1)]) 
    return (vals[1],A[vals[2],:])
end

"""
    rand_cauchy_trunc(μ::Real, c::Real)
Returns a random number with a Cauchy's distribution with **mean `μ`** and truncated between (0,1]. If the value calculated exceds 1, the returned value is 1.
If the value is less or equal to zero, the value es calculated again.
"""
function rand_cauchy_trunc(μ::Real, c::Real)
    r = c*tan(pi*(rand()-0.5)) + μ
    if r > 1
        return 1
    elseif r <= 0
        return rand_cauchy_trunc(μ::Real, c::Real)
    end
    return r
end

"""
    Lehmer_mean(X::AbstractArray,p::Real)
Return the **Lehmer mean** of degree `p` of an array `X`.
"""
function Lehmer_mean(X::AbstractArray,p::Real)
    return sum(X.^p)/sum(X.^(p-1))
end

