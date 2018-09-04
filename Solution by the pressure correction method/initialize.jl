function initialize()
# Variables:
# uₑ: speed of the lid driving the flows
# L: length of the section we are interested in
# D: the heigt of the section we are interested in
# p: the pressure, p = p′ + p_guess, where the primed variable denote a correction, in this case, a pressure correction, the asterisk denote a guessed value of the variabel in question

    # speed of the lid in m/s
    uₑ = 0.3048::Float64

    # dimensions of the simulation domain in m
    L = 0.1524::Float64
    D = 0.003048::Float64

    # density of air at sea level, in kg/m³
    ρ = 1.225055::Float64

    Re =63.6::Float64

    # dynamic viscosity of air
    # μ = (18.6*10.0^-6.0)::Float64
    #
    # # Reynolds number
    # Re = ((ρ * uₑ * D) / μ))::Float64

    # size of computational domain
    M = 5::Int64
    N = 4::Int64
    # check if the size of our matrix A is big enough to contain an interior node of both u,v, and p
    if 2M - 1 >= 9 && 2N - 1 >= 7
        # all is fine
    else
        error("M must be equal to or greater than 5, and N must be equal to or greater than 4")
    end

    # step sizes
    Δx = (L / (N - 3))::Float64
    Δy = (D / (M - 2))::Float64
    Δt = 0.001::Float64 # Δt is somewhat arbitrarily chosen, Δt works as a relaxation factor, the larger the value of Δt, the larger the change in (ρu) and (ρv) from one iteration to the next, if this change becomes to large, instabilities could arise

    v_guess  = zeros(Float64, M, N)
    u_guess  = zeros(Float64, M - 1, N - 1)
    p  = zeros(Float64, M - 1, N - 2)
    p′ = zeros(Float64, M - 1, N - 2)
    p_guess = zeros(Float64, M - 1, N - 2)
    C_guess  = zeros(Float64, 2 * M - 1, 2 * N - 1)

    # filling A with the elements of u,v and p
    C_guess[1:2:end, 1:2:end]     = v_guess[1:end, 1:end]
    C_guess[2:2:end-1, 2:2:end-1] = u_guess[1:end, 1:end]
    C_guess[2:2:end-1, 3:2:end-2] = p_guess[1:end, 1:end]

    # dirichlet boundary conditions for this specific problem
    u_guess[1, :] .= uₑ

    # artificial velocity spike to produce 2D flow during the iteration process
    v_guess[ceil(Int, M / 1.6), ceil(Int, N / 2.5)] = (uₑ/2.0)::Float64

    return u_guess, v_guess, p, p_guess, p′, ρ, Δx, Δy, Δt, C_guess
end
