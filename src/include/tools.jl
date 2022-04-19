"""
    derivative(fx, x; order=1)

Computes the `order`-th derivative of the array `fx` with respect to the space
array `x`. The derivative is computed by means of the Fast Fourier Transform provided by FFTW.jl.

Notes
=====
* Assumes periodic boundary conditions.
"""
function derivative(fx, x; order=1)
    Nx = length(x) # Number of space nodes
    dx = x[2]-x[1] # Space step

    # Frequency axis multiplied by i=√(-1)
    ik = im .* rfftfreq(Nx, 2pi/dx)

    return irfft(rfft(fx) .* ik.^order, Nx)
end
