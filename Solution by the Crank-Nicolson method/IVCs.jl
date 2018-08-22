function IVCs()
    # number of of equidistant segments of lengt Δt spanning the vertical direction of the flow
    N  = 20::Int64

    # step size in the vertical direction
    Δy = (1/N)::Float64

    # Reynolds number
    Re = 5000.0::Float64

    E  = 1.0::Float64
    # timestep value
    Δt = (E*Re*Δy^2.0)::Float64

    # Initial velocity profile
    u  = Vector{Float64}(N + 1)
    u[1:N] = 0.0
    u[end] = 1.0

    # number of equations to solved
    M = N + 1

    return M, Δy, Re, E, Δt, u
end
