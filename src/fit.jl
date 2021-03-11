function fit(
    df::DataFrame;
    max_time::Int = 10,
    fps = 0.022,
    loc_error = 0.03,
    p0::AbstractVector = [0.5, 0.5],
)
    @argcheck loc_error >= 0.0

    @. model(x, p) = 4 * p[1] * x^p[2] + 4 * loc_error^2
    result = curve_fit(
        model,
        df.delta_t[1:max_time] .* fps,
        df.msd[1:max_time],
        p0,
        lower = [0.0, 0.0],
        upper = [10.0, 2.0],
    )
    result
end

function estimate_α(df::DataFrame, δ1::Int, δ2::Int)
    @argcheck 1 <= δ2 <= δ1
    log(mean(df[df.Δt.==δ1, :msd]) / mean(df[df.Δt.==δ2, :msd])) / log(δ1 / δ2)
end
