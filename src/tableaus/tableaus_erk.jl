
"Tableau for explicit Euler method"
function getTableauExplicitEuler()
    a = zeros(Float64, 1, 1)
    b = [1.0]
    c = [0.0]
    o = 1

    TableauERK(:explicit_euler, o, a, b, c)
end

"Tableau for explicit midpoint method"
function getTableauExplicitMidpoint()
    a = [[0.0 0.0]
         [0.5 0.0]]
    b = [0.0, 1.0]
    c = [0.0, 0.5]
    o = 2

    TableauERK(:explicit_midpoint, o, a, b, c)
end

"Tableau for Runge's method"
function getTableauRunge()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    TableauERK(:runge, o, a, b, c)
end

"Tableau for Heun's method"
function getTableauHeun()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.0, 1.0]
    o = 2

    TableauERK(:heun, o, a, b, c)
end

"Tableau for Kutta's method of order three"
function getTableauKutta()
    a = [[ 0.0 0.0 0.0]
         [ 0.5 0.0 0.0]
         [-1.0 2.0 0.0]]
    b = [1/6, 4/6, 1/6]
    c = [0.0, 0.5, 1.0]
    o = 3

    TableauERK(:kutta, o, a, b, c)
end

"Tableau for explicit Runge-Kutta method of order four (1/6 rule)"
function getTableauERK4()
    a = [[0.0 0.0 0.0 0.0]
         [0.5 0.0 0.0 0.0]
         [0.0 0.5 0.0 0.0]
         [0.0 0.0 1.0 0.0]]
    b = [1/6, 1/3, 1/3, 1/6]
    c = [0.0, 0.5, 0.5, 1.0]
    o = 4

    TableauERK(:erk4, o, a, b, c)
end

"Tableau for explicit Runge-Kutta method of order four (3/8 rule)"
function getTableauERK438()
    a = [[ 0.0  0.0  0.0  0.0]
         [ 1/3  0.0  0.0  0.0]
         [-1/3  1.0  0.0  0.0]
         [ 1.0 -1.0  1.0  0.0]]
    b = [1/8, 3/8, 3/8, 1/8]
    c = [0.0, 1/3, 2/3, 1.0]
    o = 4

    TableauERK(:erk438, o, a, b, c)
end

"Tableau for explicit Verner's method of order six"
function getTableauVerner()
    a = [[ 0.0          0.0     0.0           0.0       0.0          0.0  0.0         0.0]
         [ 1/6          0.0     0.0           0.0       0.0          0.0  0.0         0.0]
         [ 4/75         16/75   0.0           0.0       0.0          0.0  0.0         0.0]
         [ 5/6         -8/3     5/2           0.0       0.0          0.0  0.0         0.0]
         [-165/64       55/6   -425/64        85/96     0.0          0.0  0.0         0.0]
         [ 12/5        -8.0     4015/612     -11/36     88/255       0.0  0.0         0.0]
         [-8263/15000   124/75 -643/680      -81/250    2484/10625   0.0  0.0         0.0]
         [ 3501/1720   -300/43  297275/52632 -319/2322  24068/84065  0.0  3850/26703  0.0]]
    b = [3/40, 0.0, 875/2244, 23/72, 264/1955, 0.0, 125/11592, 43/616]
    c = [0.0, 1/6, 4/15, 2/3, 5/6, 1.0, 1/15, 1.0]
    o = 6

    TableauERK(:verner, o, a, b, c)
end
