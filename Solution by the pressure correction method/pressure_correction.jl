 function pressure_correction(p_guess::Array{Float64}, u_guess::Array{Float64}, v_guess::Array{Float64}, ρ::Float64, timesteps::Int, Δx::Float64, Δy::Float64, Δt::Float64, C_guess::Array{Float64})

    ρu_guess = ρ .* u_guess
    ρv_guess = ρ .* v_guess
    K, L = size(C_guess)

    # size variables for the ρu_guess and the ρv_guess matrices
    Mᵤ, Nᵤ = size(u_guess)
    Mᵥ, Nᵥ = size(v_guess)
    # starts at 2 because qe are only looking at the interior nodes
    Pᵤ = 2::Int
    Qᵤ = 2::Int
    Pᵥ = 2::Int
    Qᵥ = 2::Int

    μ = 1.0::Float64

    # This is step 2 of the pressure correction method: For each timestep, iterate through matrix C_guess, which contains all the nodes for our guess of u,v, and p in the flow field, and solve for ρu_guess and ρv_guess. Step 1 is done in a separate initialization script
    for t in 1:timesteps
        for i in 3:2:K-2
            for j in 3:2:L-2
                ū = 0.5(C_guess[i + 1, j + 1] + C_guess[i - 1, j + 1])
                û = 0.5(C_guess[i + 1, j - 1] + C_guess[i - 1, j - 1])
                B_guess = -( (ρ*C_guess[i, j + 2]*ū - ρ*C_guess[i, j - 2]*û)/(2*Δx) + (ρ*C_guess[i - 2, j]^2 - ρ*C_guess[i+2, j]^2)/(2*Δy) ) + μ*( (C_guess[i, j + 2] - 2C_guess[i, j] - C_guess[i, j-2])/(Δx^2) + (C_guess[i-2, j] - 2C_guess[i, j] - C_guess[i+2, j])/(Δy^2) )

                ρv_guess[Pᵥ, Qᵥ] = ρv_guess[Pᵥ, Qᵥ] + B_guess*Δt - (Δt/Δy)*(C_guess[i-1, j] - C_guess[i+1, j])

                Qᵥ == Nᵥ-1 ? Qᵥ = 2 : Qᵥ += 1
            end
            Pᵥ == Mᵥ-1 ? Pᵥ = 2 : Pᵥ += 1
        end
        for i in 4:2:K-3
            for j in 4:2:L-3
                v̄ = 0.5(C_guess[i-1, j-1] + C_guess[i-1, j+1])
                v̂ = 0.5(C_guess[i+1, j-1] + C_guess[i+1, j+1])
                A_guess = -( (ρ*C_guess[i, j+2]^2 - ρ*C_guess[i, j-2]^2)/(2*Δx) + (ρ*C_guess[i-2, j]*v̄ - ρ*C_guess[i+2, j]*v̂)/(2*Δy) ) + μ*( (C_guess[i, j+2] - 2C_guess[i, j] + C_guess[i, j-2])/(Δx^2) + (C_guess[i-2, j] - 2C_guess[i, j] + C_guess[i+2, j])/(Δy^2) )

                ρu_guess[Pᵤ, Qᵤ] = ρu_guess[Pᵤ, Qᵤ] + A_guess*Δt - (Δt/Δx)*(C_guess[i, j+1] - C_guess[i, j-1])

                Qᵤ == Nᵤ-1 ? Qᵤ = 2 : Qᵤ += 1
            end
            Pᵤ == Mᵤ-1 ? Pᵤ = 2 : Pᵤ += 1
        end

        # step 3: Solve for p′ from the pressure correction formula
        for i in
            for j in
                b = (-Δt/Δx^2)::Float64
                c = (-Δt/Δy^2)::Float64
                a = (2(-b - c))::Float64
                d = ((1/Δx)*(ρu_guess[i, j+1] - ρu_guess[1, j-1]) +         ((1/Δy))*(ρv_guess[i-1, j] - ρv_guess[i+1, j]))::Float64
            end
        end
    end
    # return v, u, p
    return ρv_guess
end
