module ODE

# Libraries this package relies upon
using FFTW, ProgressMeter

# The objects 'exported' here become available for the user
# after importing this package
export integrate_ode, derivative, rungekutta4, rungekutta10_4, euler

# Its awkward to have everything in the same file.
# Lets split the package into multiple files!
include("include/tools.jl")
include("include/integrators.jl")
include("include/integrate_ode.jl")

end # module
