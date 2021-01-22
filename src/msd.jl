# TA-MSD
function time_average(df::DataFrame; idlabel::Symbol, xlabel::Symbol, ylabel::Symbol)
    tamsd = DataFrame(TrackID = Int[], Δt = Int[], msd = Float64[], n = Int[])
    @inbounds @threads for n in sort(collect(Set(df[!, id])))
        data = extract(df, Int(n), id, [xlabel, ylabel])
        T = size(data, 1)
        for Δt = 1:T-1
            r² = 0.0
            @simd for t = 1:(T-Δt)
                r² += squared_displacement(data, t, Δt)
            end
            r² /= (T - Δt)
            push!(tamsd, [n, Δt, r², T - Δt])
        end
    end
    tamsd
end

# EA-MSD
function ensemble_average(df::DataFrame; idlabel::Symbol, xlabel::Symbol, ylabel::Symbol)
    eamsd = DataFrame(TrackID = Int[], Δt = Int[], msd = Float64[], n = Int[])
    @inbounds @threads for n in sort(collect(Set(df[!, id])))
        data = extract(df, Int(n), idlabel, [xlabel, ylabel])
        T = size(data, 1)
        @simd for Δt = 1:T-1
            push!(eamsd, [n, Δt, squared_displacement(data, 1, Δt), 1])
        end
    end
    eamsd
end

# Ensmble average of TA-MSD
function ensemble_msd(df::DataFrame)
    eatamsd =
        DataFrame(Δt = Int[], msd = Float64[], n = Int[], std = Float64[], sem = Float64[])
    @inbounds @threads for Δt = 1:maximum(df.Δt)
        data = df[df.Δt.==Δt, :]
        push!(eatamsd, [Δt, mean(data.msd), sum(data.n), std(data.msd), sem(data.msd)])
    end
    eatamsd
end

function msd(df::DataFrame;
    idlabel::Symbol,
    xlabel::Symbol,
    ylabel::Symbol;
    return_tamsd::Bool = false,
    method,
)
    @argcheck String(method) in ["tamsd", "eamsd"]
    msd = method(df, idlabel, xlabel, ylabel)
    return_tamsd ? (ensemble_msd(msd), msd) : ensemble_msd(msd)
end
