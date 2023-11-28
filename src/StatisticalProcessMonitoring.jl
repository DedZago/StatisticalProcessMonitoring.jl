module StatisticalProcessMonitoring

using Parameters
using LinearAlgebra

include("phase2/phase2-interface.jl")
include("nominal/nominal-interface.jl")
include("statistics/stats-interface.jl")
include("limits/limits-interface.jl")
include("charts/charts-interface.jl")
include("charts/retrospective.jl")
include("algorithms/optsettings.jl")
include("algorithms/interface.jl")
include("algorithms/control-limit/sacl.jl")
include("algorithms/control-limit/bisectioncl.jl")
include("algorithms/control-limit/double-bootstrap.jl")
include("algorithms/design/grid-search.jl")
include("algorithms/design/nlopt.jl")
end
