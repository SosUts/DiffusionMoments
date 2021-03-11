function maximal_excursion(r::AbstractMatrix, τ::Integer, k::Integer)
    @argcheck (k >= 1) && (τ >= 0)
    @argcheck size(r, 1) >= 2
    T = size(r, 1)
    m = displacement(r, 1)
    @inbounds for δt = 1:τ
        mₜ = displacement(r, δt)
        if mₜ > m
            m = mₜ
        end
    end
    m^k
end

function mme(
    df::DataFrame,
    id::Symbol = :TrackID,
    x::Symbol = :POSITION_X,
    y::Symbol = :POSITION_Y,
)
    result = DataFrame(n = Int[], τ = Int[], m1 = Float64[])
    @inbounds for n = 1:maximum(df[!, id])
        m = extract(df, n, id, [x, y])
        @simd for τ = 1:size(m, 1)-1
            m1 = maximal_excursion(m, τ, 1)
            push!(result, [n, τ, m1])
        end
    end
    result
end

function StatsBase.moment(
    df::DataFrame,
    id::Symbol = :TrackID,
    x::Symbol = :POSITION_X,
    y::Symbol = :POSITION_Y,
)
    result = DataFrame(
        TrackID = Int64[],
        fourth_moment = Float64[],
        second_moment = Float64[],
        delta_t = Int64[],
        n = Int64[],
    )
    @inbounds for n in sort(collect(Set(df[!, id])))
        m = extract(df, n, id, [x, y])
        T = size(m, 1)
        for δ = 1:T-1
            c⁴ = 0.0
            c² = 0.0
            @simd for t = 1:(T-δ)
                c⁴ += abs2(abs2(displacement(m, t, δ)))
                c² += squared_displacement(m, t, δ)
            end
            c⁴ /= (T - δ)
            c² /= (T - δ)
            push!(result, [n, c⁴, c², δ, T - δ])
        end
    end
    result
end
