u_guess, v_guess, p, p_guess, p′, ρ, Δx, Δy, Δt, C_guess = initialize()

timesteps = 1
test = pressure_correction(p_guess, u_guess, v_guess, ρ, timesteps, Δx, Δy, Δt, C_guess)
