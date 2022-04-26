module GeneralAnalysis
using Reexport
@reexport using CSV, DataFrames, MixedModels,
                CategoricalArrays, Random, BrowseTables, StatsBase,
                MixedModelsExtras, MixedModelsSim,
                EffectSizes, StatsPlots

include(joinpath("Analysis", "PowerAnalysis.jl"))

export mm_effectsize

# Write your package code here.

end
