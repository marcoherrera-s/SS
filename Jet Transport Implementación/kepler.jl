using TaylorIntegration, Plots

@taylorize function kepler_eqs!(dq, q, params, t)
    dq[1] = q[3]
    dq[2] = q[4]
    rr = ( q[1]^2 + q[2]^2 )^(3/2)
    dq[3] = - q[1] / rr
    dq[4] = - q[2] / rr
end

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
t[end], q[end,:]

begin
x = view(q, :, 1)
y = view(q, :, 2)
vx = view(q, :, 3)
vy = view(q, :, 4)
plot(x, y, legend=false)
scatter!([0], [0], shape=:circle, ms=5)
xaxis!("x", (-2.0, 0.5))
yaxis!("y", (-1.0, 1.0))
title!("Fig. 1")
end


solution1 = taylorinteg(kepler_eqs!, [0.18500000000000005, 0.0, 0.0, 3.1322213859832457], time_vector, order, abstol, maxsteps=5000)

scatter(solution1[:,1], solution1[:,2])

exs = 0.81:0.001:0.9

q0s = ini_cond.(aKep, exs)

t0 = 0.0  
tf = 6.3
step = 0.4 
time_vector = t0:step:tf
order = 25
abstol = 1e-20

varorder = 9
ξ = set_variables("ξ", numvars=4, order=varorder)
q0TN = q0 .+ ξ

result = taylorinteg(kepler_eqs!, q0TN, time_vector, order, abstol, maxsteps=5000)

xjet_plot2 = map(λ->λ.(q0s), result[:,1])
vjet_plot2 = map(λ->λ.(q0s), result[:,2])

plot(xjet_plot2, vjet_plot2, legend=false)
scatter!(solution1[:,1], solution1[:,2], color=:black, alpha = 0.7)


res = []
vel = []
for i in 1:16
    xv = evaluate(result[i, :],  [0.18500000000000005, 0.0, 0.0, 3.1322213859832457])
    push!(res, xv[1])
    push!(vel, xv[2])
end

plot(res, vel)  