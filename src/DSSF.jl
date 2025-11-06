"""
    dssf(corr_data::CorrelationData, σt, σr; measure::Symbol=:sperp, order::Int=1)
Compute the dynamical structure factor from time-dependent correlation data. The time smoothing parameter is `σt`, and the spatial smoothing parameter is `σr`. The `measure` argument specifies the type of intensity measure to compute, which can be `:sperp`, `:trace`, or `:component`.
"""
function dssf(corr_data::CorrelationData, σt, σr; measure::Symbol=:sperp, order::Int=1)
    (; r0s, correlations, combiners, dt, tf, L1, L2) = corr_data

    len_r0s = length(r0s)
    ts = collect(0:dt:tf)
    ts_full = [reverse(-ts[2:end]); ts]
    energies_full = fftshift(fftfreq(length(ts_full))) * 2π / dt

    if measure == :sperp
        # For 2D scattering plane, S^{yz}, and S^{zx} does not contribute to the INS intensity
        b₁ = 4π/3 * [1, 0]
        b₂ = 2π * [-1/3, √3/3]

        @assert (1,1) in combiners && (2,2) in combiners && (3,3) in combiners && (1,2) in combiners "For :sperp measure, combiners must include (1,1), (2,2), (3,3), (1,2)"

        intensities = zeros(Float64, length(energies_full), L1, L2)
        tmp_intensities = zeros(Float64, length(energies_full), L1, L2)

        for (i, combiner) in enumerate(combiners)
            α, β = combiner
            for r0_index in eachindex(r0s)
                tmp_intensities += smooth_and_fourier_corr(correlations[:, :, i, r0_index], dt, tf, L1, L2, σt, σr; order)
            end
            for iqx in 1:L1, iqy in 1:L2
                q_2D = (iqx - 1) * b₁ / L1 + (iqy - 1) * b₂ / L2
                q = [q_2D; 0]
                q_norm = norm(q)
                if q_norm ≈ 0
                    if α == β
                        intensities[:, iqx, iqy] += tmp_intensities[:, iqx, iqy]
                    end
                else
                    if α == β
                        intensities[:, iqx, iqy] += (1-(q[α]^2)/q_norm^2) * tmp_intensities[:, iqx, iqy]
                    else
                        # Contributions from xy and yx components
                        intensities[:, iqx, iqy] += -2*(q[α]*q[β])/q_norm^2 * tmp_intensities[:, iqx, iqy]
                    end
                end
            end
        end

    elseif measure == :trace
        @assert (1,1) in combiners && (2,2) in combiners && (3,3) in combiners "For :trace measure, combiners must include (1,1), (2,2), and (3,3)."
        intensities = zeros(Float64, length(energies_full), L1, L2)
        tmp_intensities = zeros(Float64, length(energies_full), L1, L2)

        for (i, combiner) in enumerate(combiners)
            α, β = combiner
            if α == β
                for r0_index in eachindex(r0s)
                    tmp_intensities += smooth_and_fourier_corr(correlations[:, :, i, r0_index], dt, tf, L1, L2, σt, σr; order)
                end
                @. intensities += tmp_intensities
            end
        end

    elseif measure == :component
        intensities = zeros(Float64, length(combiners), length(energies_full), L1, L2)
        tmp_intensities = zeros(Float64, length(energies_full), L1, L2)
        for i in eachindex(combiners)
            tmp_intensities = zeros(Float64, length(energies_full), L1, L2)
            for r0_index in eachindex(r0s)
                tmp_intensities += smooth_and_fourier_corr(correlations[:, :, i, r0_index], dt, tf, L1, L2, σt, σr; order)
            end
            @. intensities[i, :, :, :] = tmp_intensities
        end
    else
        error("Unsupported measure: $measure. Supported measures are :sperp, :trace, :component.")
    end

    return IntensityData(intensities, energies_full, L1, L2, measure)
end