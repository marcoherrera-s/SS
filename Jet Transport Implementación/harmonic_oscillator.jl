using TaylorIntegration, Plots

u₀ = [0.0, 1.0]

t0 = 0.0  
tf = 2π
step = 0.1  
time_vector = t0:step:tf
order = 15
abstol = 1e-20

function harmonic(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
end

solution = taylorinteg(harmonic, u₀, time_vector, order, abstol, maxsteps=5000)

varorder = 13
ξ = set_variables("ξ", numvars=2, order=varorder)
q0TN = u₀ .+ ξ

result = taylorinteg(harmonic, q0TN, time_vector, order, abstol, maxsteps=5000)

polar2cart(r, ϕ) = [r*cos(ϕ), r*sin(ϕ)] # convert radius r and angle ϕ to cartesian coordinates
r = 0.05 #the radius of the neighborhood
ϕ = 0.0:0.1:(2π+0.1) #the values of the angle
ξv = polar2cart.(r, ϕ)

xjet_plot2 = map(λ->λ.(ξv), result[:,1])
vjet_plot2 = map(λ->λ.(ξv), result[:,2])

begin
    plot!(xjet_plot2, vjet_plot2, legend=false)
    scatter!(solution[:,1], solution[:,2], color=:black, alpha = 0.7, aspect_ratio=:equal)
end

# Crear rangos para q y p
q = range(-1.5, stop=1.5, length=100)'
p = range(-1.2, stop=1.2, length=100)

H3 = @. p^2/2 + q^2/2

contour(q', p, H3, levels=50, title="Curvas de nivel de H(q, p)", xlabel="q", ylabel="p", aspect_ratio=:equal)
plotly()