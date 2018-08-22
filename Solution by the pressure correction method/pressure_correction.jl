function pressure_correction(M::Int, N::Int )

    for i in 1:M
        for j in 1:N
            # calc v
            if i ≤ M - 1 && j ≤ N - 1
                # calc u
                if j ≤ N - 2
                    # calc p
                end
            end
        end
    end
    return v, u, p
end
