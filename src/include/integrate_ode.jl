function integrate_ode(dydt, y0, Nt, dt; integrator=euler)
    # Definimos el eje de tiempo
    time = range(0, step=dt, length=Nt)

    # Creamos una variable en donde guargar los valores de y para cada instante
    # trajectory[i, j] is the i-th variable at the j-th time
    trajectory = Array{Float64}(undef, length(y0), Nt)

    # Ciclo principal. Aqui integramos el valor de y desde t=0 hasta t=tf.
    # El macro @showprogress es parte de la librería ProgressMeter, y
    # nos geneera una barra de progreso informando el avance de este ciclo

    y = copy(y0) # Condicion inicial
    @showprogress for i in eachindex(time)
        # Obtenemos y(t+dt) en base a F(x, t), y(t), t y dt
        y = integrator(dydt, y, time[i], dt)

        # Guardamos el valor de y en el i-esimo instante
        trajectory[:, i] .= y
    end

    return time, trajectory
end
