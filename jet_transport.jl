using TaylorIntegration, Plots

@taylorize function cuerpo_rigido_primer_orden(dx, x, p, t)
    g, l, Ω = p
    
    # Desacoplando las ecuaciones
    dx[1] = x[2]
    dx[3] = x[4]
    
    # dx[2] calculado de forma desacoplada
    dx[2] = -(9.0*g*sin(2.0*x[3]) + 8.0*g*cos(t*Ω) + 12.0*x[4]^2*sin(x[3])) / (18.0*cos(x[3])^2 - 24.0)
    
    # dx[4] calculado usando el resultado de dx[2]
    dx[4] = -1.5*g*sin(x[3]) + 1.5 * ((9.0*g*sin(2.0*x[3]) + 8.0*g*cos(t*Ω) + 12.0*x[4]^2*sin(x[3])) / (18.0*cos(x[3])^2 - 24.0)) * cos(x[3])
end

g = 9.81  

l = 1.0   

Ω = sqrt(g)  

p = (g, l, Ω)

x0 = [0.0, 0.0, 0.0, 0.0] 


varorder = 3
ξ = set_variables("ξ", numvars=4, order=varorder)
q0TN = x0 .+ ξ

t0 = 0.0  
tf = 10.0  
step = 0.1  
time_vector = t0:step:tf
order = 5  
abstol = 1e-8 

result = taylorinteg(cuerpo_rigido_primer_orden, q0TN, time_vector, order, abstol, p)