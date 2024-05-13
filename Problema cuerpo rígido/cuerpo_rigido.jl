using Plots, TaylorIntegration, DifferentialEquations


#Resolviendo el sistema a lo ingenuo, hacemos:

params = (9.81, 1.0, sqrt(9.81))
tspan = (0.0, 10.0)

function cuerpo_rigido(ddu, du, u, p, t)
    g, l, Ω = p


    ddu[1] = (g * cos(t * Ω) / 3) + (sin(u[2]) * du[2]^2 / 2) - (cos(u[2]) * ddu[2] / 2)
    ddu[2] = -3/2 * (g * sin(u[2]) + cos(u[2]) * ddu[1]) 
end

u0 = [0.0, 0.0]
du0 = [0.0, 0.0]

prob_CR = SecondOrderODEProblem(cuerpo_rigido, du0, u0, tspan, params)
sol_CR= solve(prob_CR, KenCarp47(), maxiters=1e7)

plot(sol_CR)


# Ahora resolvemos convirtiendo a primer orden y desacoplando


x0 = [0.0, 0.0, 0.0, 0.0]

function cuerpo_rigido_primer_orden(dx, x, p, t)
    g, l, Ω = p
    
    # Desacoplando las ecuaciones
    dx[1] = x[2]
    dx[3] = x[4]
    
    # dx[2] calculado de forma desacoplada
    dx[2] = -(9.0*g*sin(2.0*x[3]) + 8.0*g*cos(t*Ω) + 12.0*x[4]^2*sin(x[3])) / (18.0*cos(x[3])^2 - 24.0)
    
    # dx[4] calculado usando el resultado de dx[2]
    dx[4] = -1.5*g*sin(x[3]) + 1.5 * ((9.0*g*sin(2.0*x[3]) + 8.0*g*cos(t*Ω) + 12.0*x[4]^2*sin(x[3])) / (18.0*cos(x[3])^2 - 24.0)) * cos(x[3])
end

prob3 = ODEProblem(cuerpo_rigido_primer_orden, x0, tspan, params)
sol3 = solve(prob3, Vern9(), maxiters=1e7)

plot(sol3)


# Dio lo mimso que resolviendo de forma desacoplada, pero se tuvo que cambiar el algoritmo utilizado


# Ahora resolvemos sin desacoplar pero iterando, como se acordó en la última reunión

function cuerpo_rigido_acoplado2(dx, x, p, t)
    g, l, Ω = p
    
    # Inicializamos dx[2] y dx[4] con un valor inicial, por ejemplo cero
    dx_2 = 0.0
    dx_4 = 0.0

    # Iteramos hasta alcanzar la convergencia en dx[2] y dx[4]
    for i in 1:50  # Un número fijo de iteraciones
        new_dx_2 = (g * cos(t * Ω) / 3) + (sin(x[3]) * x[4]^2 / 2) - (cos(x[3]) * dx_4 / 2)
        new_dx_4 = -3/2 * (g * sin(x[3]) + cos(x[3]) * new_dx_2)

        # Comprobamos la convergencia (podría ser un criterio basado en la diferencia absoluta, por ejemplo)
        if abs(new_dx_2 - dx_2) < 1e-6 && abs(new_dx_4 - dx_4) < 1e-6
            break
        end

        # Actualizamos los valores para la siguiente iteración
        dx_2 = new_dx_2
        dx_4 = new_dx_4
    end

    # Asignamos los valores calculados a dx
    dx[1] = x[2]
    dx[2] = dx_2
    dx[3] = x[4]
    dx[4] = dx_4
end

prob2 = ODEProblem(cuerpo_rigido_acoplado2, x0, tspan, params)
sol2 = solve(prob2, Vern9(), maxiters=1e7)

plot(sol2)

# Volvió a dar lo mismo, ahora usemos TaylorIntegration

t0 = 0.0
tmax = 10.0
order = 28  
abstol = 1e-20 

t, solution = taylorinteg(cuerpo_rigido_acoplado2, x0, t0, tmax, order, abstol, params)


plot(t, solution[:,1])
plot!(t, solution[:,2])
plot!(t, solution[:,3])
plot!(t, solution[:,4], xlabel="t")

# Volvió a dar lo mismo, aunque el tiempo de ejecución, por razones obivas de las iteraciones, aumentó significativamente. 

# Podemos comparar resolviendo con TaylorIntegration el sistema desacoplado.


t2, solution2 = taylorinteg(cuerpo_rigido_acoplado2, x0, t0, tmax, order, abstol, params)

plot(t2, solution2[:,1])
plot!(t2, solution2[:,2])
plot!(t2, solution2[:,3])
plot!(t2, solution2[:,4], xlabel="t")

