function thomas(A::Tridiagonal{Float64}, c::Array{Float64})
# Thomas' algorithm for solving a tridiagonal matrix of equations.
# Given a system of equations on matrix form Au = c, where A is Tridiagonal, this algorithm will eliminate the subdiagonal and solve for u and return that vector.

    # M is the number of equations in the system.
    M = size(A)[1]
    # storage for our modified main diagonal coefficients (d), the "solution" vector (c) and the vector we are solving for (u)
    c′ = Vector{Float64}(M)
    d′ = Vector{Float64}(M)
    u  = Vector{Float64}(M)

    # The first two coefficients of the vectors c′ and d′ are actually not modified at all, they are therefore set equal to their counterparts in the  input vector c, and the diagonal row in the input matrix A. The d denotes diagonal coefficients in A
    c′[1] = c[1]
    d′[1] = A[1, 1]

    for i in 2:M
        c′[i] = c[i] - (c′[i - 1] * A[i, i - 1]) / d′[i - 1]
        d′[i] = A[i, i] - (A[i, i - 1] * A[i - 1, i]) / d′[i - 1]
    end

    # backwards substitution
    u[M] = c′[M] / d′[M]
    N = M - 1
    for i in N:-1:1
        u[i] = (c′[i] - A[i, i + 1] * u[i + 1]) / d′[i]
    end
    return u
end
