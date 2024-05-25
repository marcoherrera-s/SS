using TaylorIntegration, Plots

function artificial_ode(du, u, p, t)
    du[1] = 2.0*u[1]*u[2]
    du[2] = -u[2]^2 + u[1]
end

t0 = 0.0  
tf = 1.7
step = 0.2  
time_vector = t0:step:tf
order = 15
abstol = 1e-16

solution1 = taylorinteg(artificial_ode, [0.7, -sqrt(0.7/2)], time_vector, order, abstol, maxsteps=5000)


varorder = 16
ξ = set_variables("ξ", numvars=2, order=varorder)
q0TN = [0.7, -sqrt(0.7/2)] .+ ξ


result = taylorinteg(artificial_ode, q0TN, time_vector, order, abstol, maxsteps=5000)

polar2cart(r, ϕ) = [r*cos(ϕ), r*sin(ϕ)] # convert radius r and angle ϕ to cartesian coordinates
r = 0.1 #the radius of the neighborhood
ϕ = 0.0:0.1:(2π+0.1) #the values of the angle
ξv = polar2cart.(r, ϕ)

xjet_plot2 = map(λ->λ.(ξv), result[:,1])
vjet_plot2 = map(λ->λ.(ξv), result[:,2])

plot!(xjet_plot2, vjet_plot2, legend=false)
scatter!(solution1[:,1], solution1[:,2], color=:black, alpha = 0.7)

gr()


# Crear rangos para q y p
q = range(-1, stop=1, length=100)'
p = range(-1, stop=1, length=100)



HHH = @. q * p^2 - 0.5 * q^2


# Crea el gráfico de contorno
contour(q', p, HHH, fill=true, levels=10, title="Curvas de nivel de H(q, p)", xlabel="q", ylabel="p")