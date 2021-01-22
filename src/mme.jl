function maximal_excursion(r::AbstractMatrix, Δt::Int, k::Int)
    @argcheck (k >= 1) && (Δt >= 0)
    @argcheck size(r, 1) >= 2
    T = size(r, 1)
    m = displacement(r, 1)
    @inbounds for δt = 1:Δt
        mₜ = displacement(r, δt)
        if mₜ > m
            m = mₜ
        end
    end
    m^k
end

function mme(df::DataFrame; idlabel::Symbol, xlabel::Symbol, ylabel::Symbol)
    mme = DataFrame(TrackID = Int[], Δt = Int[], m = Float64[])
    @inbounds for n = 1:maximum(df[!, idlabel])
        _data = extract(df, n, id, [xlabel, ylabel])
        @simd for Δt = 1:size(_data, 1)-1
            m = maximal_excursion(_data, Δt, 1)
            push!(result, [TrackID, Δt, m])
        end
    end
    mme
end

function moment(df::DataFrame; idlabel::Symbol, xlabel::Symbol, ylabel::Symbol)
    result = DataFrame(
        TrackID = Int[],
        Δt = Int[],
        fourth_moment = Float64[],
        second_moment = Float64[],
        n = Int[],
    )
    @inbounds for n in sort(collect(Set(df[!, idlabel])))
        m = extract(df, n, idlabel, [xlabel, ylabel])
        T = size(m, 1)
        for Δt = 1:T-1
            c⁴::Float64 = 0.0
            c²::Float64 = 0.0
            @simd for t = 1:(T-Δt)
                c⁴ += abs2(squared_displacement(m, t, Δt))
                c² += squared_displacement(m, t, Δt)
            end
            c⁴ /= (T - Δt)
            c² /= (T - Δt)
            push!(result, [n, Δt, c⁴, c², T - Δt])
        end
    end
    result
end
