function forced_oscillator(du, u, p, t)
    epsilon = 0.1  
    Omega = 6.5  
    F0 = 1.0

    du[1] = u[2]
    du[2] = -Omega^2 * u[1] + F0 * cos(Omega * (1 + epsilon) * t)
end

u₀ = [0.0, 1.0]

t0 = 0.0  
tf = 15
step = 0.1  
time_vector = t0:step:tf
order = 15
abstol = 1e-20

solution = taylorinteg(forced_oscillator, u₀, time_vector, order, abstol, maxsteps=5000)

varorder = 13
ξ = set_variables("ξ", numvars=2, order=varorder)
q0TN = u₀ .+ ξ

result = taylorinteg(forced_oscillator, q0TN, time_vector, order, abstol, maxsteps=5000)

polar2cart(r, ϕ) = [r*cos(ϕ), r*sin(ϕ)] # convert radius r and angle ϕ to cartesian coordinates
r = 0.03 #the radius of the neighborhood
ϕ = 0.0:0.1:(2π+0.1) #the values of the angle
ξv = polar2cart.(r, ϕ)

xjet_plot2 = map(λ->λ.(ξv), result[:,1])
vjet_plot2 = map(λ->λ.(ξv), result[:,2])

begin
    plot(xjet_plot2, vjet_plot2, legend=false)
    scatter!(solution[:,1], solution[:,2], color=:black, alpha = 0.7)
end

plotly()