function extract(df::DataFrame, n::Int, idlabel::Symbol, datalabel::Symbol)
    Vector(@where(df, cols(idlabel) .== n)[!, datalabel])
end

function extract(df::DataFrame, n::Int, idlabel::Symbol, datalabels::Vector{Symbol})
    Matrix(@where(df, cols(idlabel) .== n)[!, datalabels])
end

function displacement(r::AbstractMatrix, δ::Int)
    sqrt((r[1+δ, 2] - r[1, 2])^2 + (r[1+δ, 1] - r[1, 1])^2)
end

function displacement(r::AbstractMatrix, t::Int, δ::Int)
    sqrt((r[t+δ, 2] - r[t, 2])^2 + (r[t+δ, 1] - r[t, 1])^2)
end

squared_displacement(r::AbstractMatrix, t::Int, δ::Int) = displacement(r, t, δ)^2