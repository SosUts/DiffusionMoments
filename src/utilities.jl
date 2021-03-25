function displacement(r::AbstractMatrix, δ::Int)
    sqrt((r[1+δ, 2] - r[1, 2])^2 + (r[1+δ, 1] - r[1, 1])^2)
end

function displacement(r::AbstractMatrix, t::Int, δ::Int)
    sqrt((r[t+δ, 2] - r[t, 2])^2 + (r[t+δ, 1] - r[t, 1])^2)
end

squared_displacement(r::AbstractMatrix, t::Int, δ::Int) = displacement(r, t, δ)^2

function extract(df::DataFrame, n::Int, id::Symbol, x::Symbol)
    Vector(@where(df, cols(id) .== n)[!, x])
    # Vector{Float64}(df[df[:, id] .== n, x])
end

function extract(df::DataFrame, n::Int, id::Symbol, label::Vector{Symbol})
    Matrix(@where(df, cols(id) .== n)[!, label])
    # df[df[!, id] .== n, label]
end
