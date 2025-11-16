"""
    dssf(corr_data::CorrelationData, σt, σr; measure::Symbol=:sperp, order::Int=1)
Compute the dynamical structure factor from time-dependent correlation data. The time smoothing parameter is `σt`, and the spatial smoothing parameter is `σr`. The `measure` argument specifies the type of intensity measure to compute, which can be `:sperp`, `:trace`, or `:component`.
"""
function dssf(corr_data::CorrelationData, σt, σr; measure::Symbol=:sperp, order::Int=1)
    (; r0s, correlations, combiners, dt, tf, Lx, Ly) = corr_data

    ts = collect(0:dt:tf)
    ts_full = [reverse(-ts[2:end]); ts]
    energies_full = fftshift(fftfreq(length(ts_full))) * 2π / dt

    Nsites = Lx * Ly
    intensities = zeros(Float64, length(energies_full), Nsites)
    S_xx = zeros(Float64, length(energies_full), Nsites)
    S_zz = zeros(Float64, length(energies_full), Nsites)

    q_points = available_q_points(Lx, Ly)

    if measure == :sperp
        @assert (1,2) in combiners && (2,1) in combiners && (3,3) in combiners "For :sperp measure, combiners must include (1,2), (2,1), (3,3)"

        for r0_index in eachindex(r0s)
            S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 1, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 2, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            S_zz += smooth_and_fourier_corr(correlations[:, :, 3, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
        end

        for (iq, q_point) in enumerate(q_points)
            q_norm = norm(q_point)
            if q_norm ≈ 0
                intensities[:, iq] += S_xx[:, iq] + S_zz[:, iq]
            else
                intensities[:, iq] += (2-q_point[1]^2/q_norm^2-q_point[2]^2/q_norm^2) * S_xx[:, iq] + S_zz[:, iq]
            end
        end

    elseif measure == :trace
        @assert (1,2) in combiners && (2,1) in combiners && (3,3) in combiners "For :trace measure, combiners must include (1,2), (2,1), and (3,3)."

        for r0_index in eachindex(r0s)
            S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 1, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 2, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            S_zz += smooth_and_fourier_corr(correlations[:, :, 3, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
        end

        @. intensities += S_xx + S_zz

    elseif measure == :component
        if length(combiners) == 2
            for r0_index in eachindex(r0s)
                S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 1, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
                S_zz += smooth_and_fourier_corr(correlations[:, :, 2, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            end
        elseif length(combiners) == 3
            for r0_index in eachindex(r0s)
                S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 1, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
                S_xx += 0.25 * smooth_and_fourier_corr(correlations[:, :, 2, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
                S_zz += smooth_and_fourier_corr(correlations[:, :, 3, r0_index], r0s[r0_index], Lx, Ly, dt, tf, σt, σr; order)
            end
        end
    else
        error("Unsupported measure: $measure. Supported measures are :sperp, :trace, :component.")
    end

    return IntensityData(intensities, S_xx, S_zz, q_points, energies_full, Lx, Ly, measure)
end