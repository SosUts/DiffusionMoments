"""
    mean_square_disaplcement(df; method=[:ensemble_average, time_average], non_averaged::Bool) -> DataFrame, [DataFrame]

Compute mean_square_disaplcement over files.

**Example**
```julia
using SPT

```
"""
function _ensemble_time_average_msd(
    df::DataFrame,
    id::Symbol,
    x::Symbol,
    y::Symbol;
    return_tamsd::Bool = false,
)
    ta_msd = time_average_msd(df, id, x, y)
    return_tamsd ? (ensemble_tamsd(ta_msd), ta_msd) : ensemble_tamsd(ta_msd)
end

function _time_average_msd(
    df::DataFrame,
    id::Symbol = :TrackID,
    x::Symbol = :POSITION_X,
    y::Symbol = :POSITION_Y;
)
    tamsd = DataFrame(TrackID = Int64[], cell_name = String[], msd = Float64[], delta_t = Int64[], n = Int64[])
    @inbounds Threads.@threads for cell = sort(unique(df.cell_name))
        _data = df[df.cell_name .== cell, :]
        for n = sort(unique(_data[!, id]))
            data = extract(_data, Int(n), id, [x, y])
            T = size(data, 1)
            for δ = 1:T-1
                r² = 0.0
                @simd for t = 1:(T-δ)
                    r² += squared_displacement(data, t, δ)
                end
                r² /= (T - δ)
                push!(tamsd, [n, cell, r², δ, T - δ])
            end
        end
    end
    tamsd
end


function _ensemble_msd(df::DataFrame, id::Symbol, x::Symbol, y::Symbol)
    eamsd = DataFrame(TrackID = Int64[], cell_name = String[], msd = Float64[], n = Int[], delta_t = Int64[])
    @inbounds Threads.@threads for cell = sort(unique(df.cell_name))
        _data = df[df.cell_name .== cell, :]
        for n = sort(unique(_data[!, id]))
            data = extract(_data, Int(n), id, [x, y])
            T = size(data, 1)
            @simd for δ = 1:T-1
                push!(eamsd, [n, cell, squared_displacement(data, 1, δ), 1, δ])
            end
        end
    end
    eamsd
end

function _ensemble_tamsd(df::DataFrame)
    eatamsd = DataFrame(
        cell_name = String[],
        delta_t = Float64[],
        msd = Float64[],
        n = Int64[],
        std = Float64[],
        sem = Float64[],
        ci = Float64[],
    )
    @inbounds Threads.@threads for cell = sort(unique(df.cell_name))
        _data = df[df.cell_name .== cell, :]
        for i = 1:maximum(_data.delta_t)
            data = _data[_data.delta_t .== i, :]
            push!(
                eatamsd,
                [
                    cell,
                    i,
                    mean(data.msd),
                    sum(data.n),
                    std(data.msd),
                    sem(data.msd),
                    1.96 * sem(data.msd),
                ],
            )
        end
    end
    eatamsd
end