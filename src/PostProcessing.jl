function smooth_and_fourier_corr(corr::Matrix{ComplexF64}, dt, tf, L1::Int, L2::Int, σt, σr; order::Int=1)
    a₁ = [3/2, √3/2]
    a₂ = [0, √3]

    Nu = L1 * L2
    corr_full = [reverse(conj(corr[2:end, :]), dims=1); corr]
    corr_reshaped = zeros(ComplexF64, size(corr_full, 1), L1, L2)

    for site in 1:2Nu
        x, y, _ = inverse_siteindex(L2, site; order)
        # Note the sublattice averaging
        corr_reshaped[:, x, y] += corr_full[:, site]
    end

    ts = collect(0:dt:tf)
    ts_full = [reverse(-ts[2:end]); ts]
    # Smooth in time
    for (it, t) in enumerate(ts_full)
        corr_reshaped[it, :, :] *= exp(-σt * t^2)
    end
    # Smooth in space
    r0 = (L1÷2+1, L2÷2+1)
    for i in 1:L1, j in 1:L2
        rvec = (i - r0[1]) * a₁ + (j - r0[2]) * a₂
        r = norm(rvec)
        corr_reshaped[:, i, j] *= exp(-σr * r^2)
    end

    corr_qt = fft(corr_reshaped, (2, 3)) / Nu
    corr_qω = ifft(corr_qt, 1) * length(ts_full)
    corr_qω = abs.(fftshift(corr_qω, 1))

    return corr_qω
end