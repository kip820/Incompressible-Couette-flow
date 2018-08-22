function initialize()
# Variables:
# uₑ: speed of the lid driving the flows
# L: length of the section we are interested in
# D: the heigt of the section we are interested in
# p: the pressure, p = p′ + p✳, where the primed variable denote a correction, in this case, a pressure correction, the asterisk denote a guessed value of the variabel in question

    # in m/s
    uₑ = 0.3048::Float64

    # in m
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
    M = 23::Int64
    N = 12::Int64

    if M*N ≥ 5000
        v = zeros(Float64, M, N)
        sizehint!(v, M*N)
        u = zeros(Float64, M - 1, N - 1)
        sizehint!(u, M*N)
        p = zeros(Float64, M - 1, N - 2)
        sizehint!(p, M*N)
        p′ = zeros(Float64, M - 1, N - 2)
        sizehint!(p′, M*N)
        p✳ = zeros(Float64, M - 1, N - 2)
        sizehint!(p✳, M*N)
    else
        v = zeros(Float64, (M, N)
        u = zeros(Float64, (M - 1, N - 1)
        p = zeros(Float64, (M - 1, N - 2)
        p′ = zeros(Float64, M - 1, N - 2)
        p✳ = zeros(Float64, M - 1, N - 2)
    end

    # dirichlet boundary conditions for this specific problem
    u[1, :] = uₑ

    # artificial velocity spike to produce 2D flow during the iteration process
    v[15, 5] = (uₑ/2.0)::Float64

    return M, N, u, v, p, p✳, p′
end
