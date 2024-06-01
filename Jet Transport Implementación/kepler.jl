using TaylorIntegration, Plots

@taylorize function kepler_eqs!(dq, q, params, t)
    dq[1] = q[3]
    dq[2] = q[4]
    rr = ( q[1]^2 + q[2]^2 )^(3/2)
    dq[3] = - q[1] / rr
    dq[4] = - q[2] / rr
end

t0 = 0.0  
tf = 6.3
step = 0.3
time_vector = t0:step:tf
order = 25
abstol = 1e-20

const mu = 1.0
const mass = 1.0
const aKep = 1.0
const eKep = 0.8

function ini_cond(a, e)
    x0  = a*(one(e)-e)
    vy0 = mass * sqrt( mu * a * (1-e^2) ) / x0
    y0  = zero(vy0)
    vx0 = zero(vy0)
    return [x0, y0, vx0, vy0]
end
q0 = ini_cond(aKep, eKep)

t, q = taylorinteg(kepler_eqs!, [0.18500000000000005, 0.0, 0.0, 3.1322213859832457], 0.0, 10, 25, 1.0e-20);

begin
    x = view(q, :, 1)
    y = view(q, :, 2)
    vx = view(q, :, 3)
    vy = view(q, :, 4)
    plot(x, y, legend=false)
    # scatter!([0], [0], shape=:circle, ms=5)
    # xaxis!("x", (-2.0, 0.5))
    # yaxis!("y", (-1.0, 1.0))
    # title!("Fig. 1")
end

solution1 = taylorinteg(kepler_eqs!, q0, time_vector, order, abstol, maxsteps=5000)

scatter(solution1[:,1], solution1[:,2])

varorder = 5
ξ = set_variables("ξ", numvars=4, order=varorder)
q0TN = q0 .+ ξ

result = taylorinteg(kepler_eqs!, q0TN, time_vector, order, abstol, maxsteps=5000)

exs = -0.007:0.001:0.007

q0s =[ini_cond.(aKep, eKep .+ exs)[i] .- q0 for i in 1:length(exs)]


xjet_plot2 = map(λ->λ.(q0s), result[:,1])
vjet_plot2 = map(λ->λ.(q0s), result[:,2])

plot(xjet_plot2, vjet_plot2, legend=false)
scatter!(solution1[:,1], solution1[:,2], color=:black, alpha = 0.7)


function polar2cart(r, ϕ)
    return [r*cos(ϕ), r*sin(ϕ)]
end

r_pos = 0.00001
r_vel = 0.00001

ϕ_pos = 0.0:0.5:(2π+0.1)
ϕ_vel = 0.0:0.5:(2π+0.1)

ξ_pos = [polar2cart(r_pos, ϕ) for ϕ in ϕ_pos]
ξ_vel = [polar2cart(r_vel, ϕ) for ϕ in ϕ_vel]

ξ = [(pos[1], pos[2], vel[1], vel[2]) for pos in ξ_pos for vel in ξ_vel]

q0s2 =[ξ[i] .- q0 for i in 1:length(ξ)]

xjet_plot22 = map(λ->λ.(q0s2), result[:,1])
vjet_plot22 = map(λ->λ.(q0s2), result[:,2])

plot(xjet_plot22, vjet_plot22, legend=false)
scatter!(solution1[:,1], solution1[:,2], color=:black, alpha = 0.7)


