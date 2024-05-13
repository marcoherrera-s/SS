using TaylorIntegration, Plots

function harmonic(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
end

u₀ = [0, 1.0]
varorder = 8
ξ = set_variables("ξ", numvars=2, order=varorder)
q0TN = u₀ .+ ξ

t0 = 0.0  
tf = 15.0  
step = 0.01  
time_vector = t0:step:tf
order = 8  
abstol = 1e-10 

result = taylorinteg(harmonic, q0TN, time_vector, order, abstol)

n_rows, n_cols = size(result)

# Lista para guardar los vectores de resultados
results_by_time = []

for time in time_vector
    # Vector para guardar los resultados de todos los polinomios evaluados en 'time'
    current_time_results = Vector{Float64}(undef, n_rows * n_cols)

    # Índice para colocar los resultados en el vector
    index = 1
    for i in 1:n_rows
        for j in 1:n_cols
            current_time_results[index] = evaluate(result[i, j], [time, time])
            index += 1
        end
    end
    # Agregar el vector de resultados a la lista
    push!(results_by_time, current_time_results)
end

anim = @animate for i in 1:500
    plot(time_vector, results_by_time[i][1502:3002], alpha=0.6, label=false, ylim=(-10,10), ylabel="Posición", xlabel="tiempo") #Selecciono los valores asociados a la posición, del 1-1501 son los valores asociados a la velocidad
end
gif(anim, "anim_oscilador.gif", fps = 60)  