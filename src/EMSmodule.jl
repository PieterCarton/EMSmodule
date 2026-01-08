module EMSmodule

using JuMP, InfiniteOpt, LinearAlgebra, Distributions, Interpolations
using Parameters, Random
using JSON3, MAT, CSV, StatsBase, DataFrames, Dates, LiiBRA, GaussianMixtures, Serialization
using LaTeXStrings, ColorSchemes
using XLSX, StateSpaceModels, Statistics, TimeSeries
using Makie, CairoMakie, GLMakie

include("makeEMSobjs.jl")
include("EMSfns.jl")
include("ESSfns.jl")
include("EMSrfofns.jl")
include("thermfns.jl")
include("makeEMSplots.jl")
include("simEMS.jl")
include("testEMS.jl")
include("makeForecasts.jl")

Random.seed!(1234);

# automatically export all functions and types
# export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end
# # manual export
# export makeInputsplot, makeEBplot, makeEMSplots, compareEB, compareEMSplots, compareTESS, compareEMSaging, hist2d, countourSP, plotWeightedCost, compareCgrid, compareDiffCgrid, optStatusPlots
# export ....

end # module EMSmodule
