function available_q_points(Lx::Int, Ly::Int)
    ex = 2π / Lx
    ey = 4π / (√3*Ly)
    q_points = [[ex*n, ey*m] for n in 0:Lx-1, m in 0:Ly-1]
    return q_points
end

function smooth_and_fourier_corr(corr::Matrix{ComplexF64}, site0::Int, Lx::Int, Ly::Int, dt, tf, σt; order::Int=1)
    Ns = Lx * Ly
    corr_full = [reverse(conj(corr[2:end, :]), dims=1); corr]
    x0, y0 = inverse_siteindex(Lx, Ly, site0; order=order)
    rx0, ry0 = cartesian_coordinates(x0, y0)

    rs = Vector{Tuple{Float64, Float64}}(undef, Ns)
    for site in 1:Ns
        x, y = inverse_siteindex(Lx, Ly, site; order=order)
        rs[site] = cartesian_coordinates(x, y)
    end

    ts = collect(0:dt:tf)
    ts_full = [reverse(-ts[2:end]); ts]

    # Smooth in time
    for (it, t) in enumerate(ts_full)
        corr_full[it, :] *= exp(-σt * t^2)
    end

    # Fourier transform along the spatial directions
    corr_qt = zero(corr_full)
    q_points = available_q_points(Lx, Ly)
    for iq in 1:Ns
        qx, qy = q_points[iq]
        for site in 1:Ns
            rx, ry = rs[site]
            phase = exp(-1im * (qx*(rx - rx0) + qy*(ry - ry0)))
            corr_qt[:, iq] += corr_full[:, site] * phase
        end
    end
    corr_qω = ifft(corr_qt, 1) * length(ts_full)
    corr_qω = abs.(fftshift(corr_qω, 1))

    return corr_qω
end