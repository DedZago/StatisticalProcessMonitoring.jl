module SPM

using Parameters

include("phase1/phase1-interface.jl")
include("nominal/nominal-interface.jl")
include("statistics/stats-interface.jl")
include("limits/limits-interface.jl")
include("charts/charts-interface.jl")
include("charts/retrospective.jl")
include("algorithms/sacl.jl")
include("algorithms/bisectioncl.jl")
include("algorithms/grid-search.jl")
end
