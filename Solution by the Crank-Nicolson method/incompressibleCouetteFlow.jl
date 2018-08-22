function incompressibleCouetteFlow(timesteps::Int64)

    # read in IVCs and BCs
    N, Δy, Re, E, Δt, u = IVCs()
    # we have N grid points, there's one equation to be solved for each grid point, since we have Dirichlet BCs, we know the value of the first, and the last u (the variable of interest) for all timesteps. However, the value of the second u to the second to last value of u needs to be solved for each timestep, i.e. we are only really solving M = N - 2 equations
    M = (N - 2)::Int64

    # calculation of the coefficients in the tridiagonal matrix, typical of linear diffusion equations
    A = (-E/2.0)::Float64
    B = (1.0 + E)::Float64
    A = A.*ones(Float64, M - 1)
    B = B.*ones(Float64, M)
    K = ones(Float64, M)

    H = Tridiagonal(A, B, A)

    for n in 1:timesteps
        # calculate the uⁿ values for all grid points except the first and the last grid point (for which we already know the value of uⁿ)
        for j in 2:M+1
            K[j - 1] = (1.0 - E) * u[j] + (E/2.0) * (u[j + 1] + u[j - 1])
        end
        K[M] += - A[M - 1]*u[N]

        # We have Dirichlet type BCs -> we are only interested in the values of the flow variable u between these BCs
        u[2:end-1] = thomas(H, K)
    end

    return u
end
