var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#GeometricIntegrators.jl-1",
    "page": "Home",
    "title": "GeometricIntegrators.jl",
    "category": "section",
    "text": "Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.(Image: Build Status) (Image: Coverage Status) (Image: codecov)GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles."
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\"tutorial.md\"]"
},

{
    "location": "index.html#Modules-1",
    "page": "Home",
    "title": "Modules",
    "category": "section",
    "text": "Pages = [\"modules/basis_functions.md\",\n         \"modules/equations.md\",\n         \"modules/integrators.md\",\n         \"modules/interpolation.md\",\n         \"modules/solvers_linear.md\",\n         \"modules/solvers_nonlinear.md\",\n         \"modules/solutions.md\",\n         \"modules/tableaus.md\"\n]"
},

{
    "location": "index.html#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "Copyright (c) 2016 Michael Kraus <michael.kraus@ipp.mpg.de>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

{
    "location": "tutorial.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "modules/basis_functions.html#",
    "page": "Basis Functions",
    "title": "Basis Functions",
    "category": "page",
    "text": ""
},

{
    "location": "modules/basis_functions.html#Basis-Functions-1",
    "page": "Basis Functions",
    "title": "Basis Functions",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.BasisFunctions]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/equations.html#",
    "page": "Equations",
    "title": "Equations",
    "category": "page",
    "text": ""
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.DAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.DAE",
    "category": "Type",
    "text": "DAE: Differential Algebraic Equation\n\nDefines a differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t)) + u(t q(t) lambda(t))   q(t_0) = q_0  \n0 = phi (t q(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector field f, projection u, algebraic constraint phi=0, initial conditions q_0 and lambda_0, the dynamical variable q taking values in mathbbR^m and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nm: dimension of dynamical variable q and the vector field f\nn: dimension of algebraic variable lambda and the constraint phi\nv: function computing the vector field\nu: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\nλ₀: initial condition for algebraic variable lambda\n\nThe function v, providing the vector field, takes three arguments, v(t, q, v), the functions u and ϕ, providing the projection and the algebraic constraint take four arguments, u(t, q, λ, u) and ϕ(t, q, λ, ϕ), where t is the current time, q and λ are the current solution vectors, and v, u and ϕ are the vectors which hold the result of evaluating the vector field v, the projection u and the algebraic constraint phi on t, q and λ.\n\nExample\n\n    function v(t, q, v)\n        v[1] = q[1]\n        v[2] = q[2]\n    end\n\n    function u(t, q, λ, u)\n        u[1] = +λ[1]\n        u[2] = -λ[1]\n    end\n\n    function ϕ(t, q, λ, ϕ)\n        ϕ[1] = q[2] - q[1]\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ₀ = [0.]\n\n    dae = DAE(v, u, ϕ, t₀, q₀, λ₀)\n\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.IDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.IDAE",
    "category": "Type",
    "text": "IDAE: Implicit Differential Algebraic Equation\n\nDefines a partitioned differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t) + u(t q(t) p(t) lambda(t))   q(t_0) = q_0  \ndotp (t) = f(t q(t) v(t)) + r(t q(t) p(t) lambda(t))   p(t_0) = p_0  \np(t) = p(t q(t) v(t))   \n0 = phi (t q(t) p(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector field f, the momentum defined by p, projection u and r, algebraic constraint phi=0, conditions (q_0 p_0) and lambda_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nm: dimension of dynamical variables q and p as well as the vector fields f and p\nn: dimension of algebraic variable lambda and the constraint phi\nf: function computing the vector field f\np: function computing p\nu: function computing the projection\ng: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\np₀: initial condition for dynamical variable p\nλ₀: initial condition for algebraic variable lambda\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.IODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.IODE",
    "category": "Type",
    "text": "IODE: Implicit Ordinary Differential Equation\n\nDefines an implicit initial value problem\n\nbeginalign*\ndotq (t) = v(t)  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) v(t))  \np(t_0) = p_0  \np(t) = (t q(t) v(t))\nendalign*\n\nwith vector field f, the momentum defined by p, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d. This is a special case of a differential algebraic equation with dynamical variables (qp) and algebraic variable v.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nα: function determining the momentum\nf: function computing the vector field\ng: function determining the projection, given by ∇α(q)λ\nv: function computing an initial guess for the velocity field (optional)\nt₀: initial time (optional)\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions α and f must have the interface\n\n    function α(t, q, v, p)\n        p[1] = ...\n        p[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, v, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, v is the current velocity and f and p are the vectors which hold the result of evaluating the functions f and  on t, q and v. The funtions g and v are specified by\n\n    function g(t, q, λ, g)\n        g[1] = ...\n        g[2] = ...\n        ...\n    end\n\nand\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.ODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.ODE",
    "category": "Type",
    "text": "ODE: Ordinary Differential Equation\n\nDefines an initial value problem\n\ndotq (t) = v(t q(t))  qquad q(t_0) = q_0 \n\nwith vector field v, initial condition q_0 and the solution q taking values in mathbbR^d.\n\nFields\n\nd: dimension of dynamical variable q and the vector field v\nv: function computing the vector field\nt₀: initial time\nq₀: initial condition\n\nThe function v providing the vector field must have the interface\n\n    function v(t, q, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, and v is the vector which holds the result of evaluating the vector field v on t and q.\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.PDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PDAE",
    "category": "Type",
    "text": "PDAE: Partitioned Differential Algebraic Equation\n\nDefines a partitioned differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t) p(t)) + u(t q(t) p(t) lambda(t))   q(t_0) = q_0  \ndotp (t) = f(t q(t) p(t)) + r(t q(t) p(t) lambda(t))   p(t_0) = p_0  \n0 = phi (t q(t) p(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector fields v and f, projection u and r, algebraic constraint phi=0, conditions (q_0 p_0) and lambda_0, the dynamical variables (qp) taking values in mathbbR^d times mathbbR^d and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nm: dimension of dynamical variables q and p as well as the vector fields v and f\nn: dimension of algebraic variable lambda and the constraint phi\nv: function computing the vector field v\nf: function computing the vector field f\nu: function computing the projection\ng: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\np₀: initial condition for dynamical variable p\nλ₀: initial condition for algebraic variable lambda\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.PODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PODE",
    "category": "Type",
    "text": "IODE: Partitioned Ordinary Differential Equation\n\nDefines a partitioned initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t) p(t))  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) p(t))  \np(t_0) = p_0 \nendalign*\n\nwith vector fields v and f, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields v and f\nv: function computing the vector field v\nf: function computing the vector field f\nt₀: initial time\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions v and f must have the interface\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, p, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q and p are the current solution vectors and v and f are the vectors which hold the result of evaluating the vector fields v and f on t, q and p.\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.VODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.VODE",
    "category": "Type",
    "text": "VODE: Variational Ordinary Differential Equation\n\nDefines an implicit initial value problem\n\nbeginalign*\ndotq (t) = v(t)  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) v(t))  \np(t_0) = p_0  \np(t) = (t q(t) v(t))\nendalign*\n\nwith vector field f, the momentum defined by p, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d. This is a special case of a differential algebraic equation with dynamical variables (qp) and algebraic variable v.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and p\nα: function determining the momentum\nf: function computing the vector field\ng: function determining the projection, given by ∇α(q)λ\nv: function computing an initial guess for the velocity field (optional)\nt₀: initial time (optional)\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions α and f must have the interface\n\n    function α(t, q, v, p)\n        p[1] = ...\n        p[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, v, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, v is the current velocity and f and p are the vectors which hold the result of evaluating the functions f and  on t, q and v. The funtions g and v are specified by\n\n    function g(t, q, λ, g)\n        g[1] = ...\n        g[2] = ...\n        ...\n    end\n\nand\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\n\n\n"
},

{
    "location": "modules/equations.html#Equations-1",
    "page": "Equations",
    "title": "Equations",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Equations]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/integrators.html#",
    "page": "Integrators",
    "title": "Integrators",
    "category": "page",
    "text": ""
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.DAE,GeometricIntegrators.Tableaus.TableauARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.DAE,GeometricIntegrators.Tableaus.TableauSARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for special additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.Equation,GeometricIntegrators.Tableaus.AbstractTableau,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Print error for integrators not implemented, yet.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.IDAE,GeometricIntegrators.Tableaus.TableauVPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for variational partitioned additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.IDAE,GeometricIntegrators.Tableaus.TableauVSPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for variational special partitioned additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.IODE,GeometricIntegrators.Tableaus.TableauVPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for variational partitioned Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Tableaus.TableauDIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for diagonally implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Tableaus.TableauERK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for explicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Tableaus.TableauFIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for fully implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Tableaus.TableauSIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for singly implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PDAE,GeometricIntegrators.Tableaus.TableauPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for partitioned additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PDAE,GeometricIntegrators.Tableaus.TableauSPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for special partitioned additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PODE,GeometricIntegrators.Tableaus.TableauEPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for explicit partitioned Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PODE,GeometricIntegrators.Tableaus.TableauIPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for implicit partitioned Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDIRK",
    "category": "Type",
    "text": "Diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorEPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorEPRK",
    "category": "Type",
    "text": "Explicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorERK",
    "category": "Type",
    "text": "Explicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorFIRK",
    "category": "Type",
    "text": "Fully implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorIPRK",
    "category": "Type",
    "text": "Implicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSARK",
    "category": "Type",
    "text": "Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSIRK",
    "category": "Type",
    "text": "Singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSPARK",
    "category": "Type",
    "text": "Special Partitioned Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPARK",
    "category": "Type",
    "text": "Variational partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRK",
    "category": "Type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVPRKpMidpoint",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpMidpoint",
    "category": "Type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVPRKpStandard",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpStandard",
    "category": "Type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVPRKpSymmetric",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVPRKpSymmetric",
    "category": "Type",
    "text": "Variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorVSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorVSPARK",
    "category": "Type",
    "text": "Variational special partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.initialize!-Tuple{GeometricIntegrators.Integrators.Integrator,GeometricIntegrators.Solutions.Solution,Int64,Int64}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.initialize!",
    "category": "Method",
    "text": "Initialize integrator for initial conditions m with m₁ ≤ m ≤ m₂ and time step 0.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Apply integrator for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Integrate given equation with given tableau for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-NTuple{4,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{Any,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE for all initial conditions.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorARK,GeometricIntegrators.Solutions.SolutionDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate DAE with Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorDIRK,GeometricIntegrators.Solutions.SolutionODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorDIRK,GeometricIntegrators.Solutions.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSARK,GeometricIntegrators.Solutions.SolutionDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate DAE with Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSIRK,GeometricIntegrators.Solutions.SolutionODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSIRK,GeometricIntegrators.Solutions.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSPARK,GeometricIntegrators.Solutions.SolutionPDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned DAE with Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVPARK{DT,TT,FT,PT,UT,GT,ϕT,VT,SPT,ST,IT} where IT where ST where SPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{VT}, Tuple{ϕT}} where N where VT where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate DAE with variational partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Union{Tuple{DT}, Tuple{GeometricIntegrators.Integrators.Integrator{DT,TT},GeometricIntegrators.Solutions.Solution{DT,TT,N},Int64,Int64,Int64,Int64}, Tuple{N}, Tuple{TT}} where N where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE for initial conditions m with m₁ ≤ m ≤ m₂ for time steps n with n₁ ≤ n ≤ n₂.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{DT,1},Array{DT,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersPARK{DT,TT,FT,PT,UT,GT,ϕT}}, Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{ϕT}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{DT,1},Array{DT,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVPARK{DT,TT,FT,PT,UT,GT,ϕT}}, Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{ϕT}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of variational partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{DT,1},Array{DT,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVSPARK{DT,TT,FT,PT,UT,GT,ϕT}}, Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{ϕT}} where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of variational special partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK{DT,TT,VT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{VT}} where S where D where VT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of fully implicit Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK{DT,TT,VT,FT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{FT}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{VT}} where S where D where FT where VT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of implicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpMidpoint{DT,TT,ΑT,FT,GT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{FT}, Tuple{GT}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{ΑT}} where S where D where GT where FT where ΑT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpStandard{DT,TT,ΑT,GT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{GT}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{ΑT}} where S where D where GT where ΑT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of projected variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpSymmetric{DT,TT,ΑT,FT,GT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{FT}, Tuple{GT}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{ΑT}} where S where D where GT where FT where ΑT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Union{Tuple{Array{ST,1},Array{ST,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S}}, Tuple{DT}, Tuple{D}, Tuple{FT}, Tuple{ST}, Tuple{S}, Tuple{TT}, Tuple{ΑT}} where S where D where FT where ΑT where TT where DT where ST",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorARK",
    "category": "Type",
    "text": "Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorPARK",
    "category": "Type",
    "text": "Implicit partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK",
    "category": "Type",
    "text": "Parameters for right-hand side function of fully implicit Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK",
    "category": "Type",
    "text": "Parameters for right-hand side function of implicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersPARK",
    "category": "Type",
    "text": "Parameters for right-hand side function of partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVPARK",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRK",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpMidpoint",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpMidpoint",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpStandard",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpStandard",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpSymmetric",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVPRKpSymmetric",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersVSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersVSPARK",
    "category": "Type",
    "text": "Parameters for right-hand side function of variational special partitioned additive Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.aitken_neville-Union{Tuple{Array{TT,1},Array{DT,2},TT,Array{DT,1}}, Tuple{DT}, Tuple{TT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.aitken_neville",
    "category": "Method",
    "text": "Compute p(x) where p is the unique polynomial of degree length(xi), such that p(x[i]) = y[i]) for all i.\n\nti: interpolation nodes\nxi: interpolation values\nt:  evaluation point\nx:  evaluation value\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.computeStageP!-Tuple{GeometricIntegrators.Integrators.IntegratorEPRK,Int64,Int64,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageP!",
    "category": "Method",
    "text": "Compute P stages of explicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.computeStageQ!-Tuple{GeometricIntegrators.Integrators.IntegratorEPRK,Int64,Int64,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageQ!",
    "category": "Method",
    "text": "Compute Q stages of explicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.euler_extrapolation-Union{Tuple{DT}, Tuple{Function,TT,TT,Array{DT,1},Array{DT,1},Int64}, Tuple{TT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.euler_extrapolation",
    "category": "Method",
    "text": "Euler extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nt₀: initial time\nt₁: final   time\nx₀: initial value\nx₁: final   value\ns:  number of interpolations (order p=s+1)\n\nTODO This is probably broken!\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorPARK{DT,TT,FT,PT,UT,GT,ϕT,ST} where ST,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{ϕT}} where N where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate DAE with partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVPRKpMidpoint{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} where IT where ST where FPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}, Tuple{ΑT}} where N where VT where GT where FT where ΑT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVPRKpStandard{DT,TT,ΑT,FT,GT,VT,SPT,PPT,SST,STP,IT} where IT where STP where SST where PPT where SPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}, Tuple{ΑT}} where N where VT where GT where FT where ΑT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVPRKpSymmetric{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} where IT where ST where FPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}, Tuple{ΑT}} where N where VT where GT where FT where ΑT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVPRK{DT,TT,ΑT,FT,GT,VT,FPT,ST,IT} where IT where ST where FPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}, Tuple{ΑT}} where N where VT where GT where FT where ΑT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with variational partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GT}, Tuple{GeometricIntegrators.Integrators.IntegratorVSPARK{DT,TT,FT,PT,UT,GT,ϕT,VT,SPT,ST,IT} where IT where ST where SPT,GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{PT}, Tuple{TT}, Tuple{UT}, Tuple{VT}, Tuple{ϕT}} where N where VT where ϕT where GT where UT where PT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate DAE with variational special partitioned additive Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GeometricIntegrators.Integrators.IntegratorEPRK{DT,TT,VT,FT},GeometricIntegrators.Solutions.SolutionPODE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}} where N where FT where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GeometricIntegrators.Integrators.IntegratorERK{DT,TT,FT},GeometricIntegrators.Solutions.SolutionODE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}} where N where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with explicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GeometricIntegrators.Integrators.IntegratorFIRK{DT,TT,FT,SPT,ST,IT,N} where N,GeometricIntegrators.Solutions.SolutionODE{DT,TT,N},Int64,Int64}, Tuple{IT}, Tuple{N}, Tuple{SPT}, Tuple{ST}, Tuple{TT}} where N where IT where ST where SPT where FT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with fully implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate_step!-Union{Tuple{DT}, Tuple{FT}, Tuple{GeometricIntegrators.Integrators.IntegratorIPRK{DT,TT,VT,FT,SPT,ST,IT} where IT where ST where SPT,GeometricIntegrators.Solutions.SolutionPODE{DT,TT,N},Int64,Int64}, Tuple{N}, Tuple{TT}, Tuple{VT}} where N where FT where VT where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate_step!",
    "category": "Method",
    "text": "Integrate ODE with implicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.midpoint_extrapolation-Union{Tuple{DT}, Tuple{Function,Function,TT,TT,Array{DT,1},Array{DT,1},Array{DT,1},Array{DT,1},Int64}, Tuple{TT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.midpoint_extrapolation",
    "category": "Method",
    "text": "Midpoint extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nf:  function to compute force  field\nt₀: initial time\nt₁: final   time\nq₀: initial positions\np₀: initial momenta\nq₁: final   positions\np₁: final   momenta\ns:  number of interpolations (order p=2s+2)\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.midpoint_extrapolation-Union{Tuple{DT}, Tuple{Function,TT,TT,Array{DT,1},Array{DT,1},Int64}, Tuple{TT}} where TT where DT",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.midpoint_extrapolation",
    "category": "Method",
    "text": "Midpoint extrapolation method with arbitrary order p.\n\nv:  function to compute vector field\nt₀: initial time\nt₁: final   time\nx₀: initial value\nx₁: final   value\ns:  number of interpolations (order p=2s+2)\n\n\n\n"
},

{
    "location": "modules/integrators.html#Integrators-1",
    "page": "Integrators",
    "title": "Integrators",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Integrators]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/interpolation.html#",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "page",
    "text": ""
},

{
    "location": "modules/interpolation.html#GeometricIntegrators.Interpolation.HermiteInterpolation",
    "page": "Interpolation",
    "title": "GeometricIntegrators.Interpolation.HermiteInterpolation",
    "category": "Type",
    "text": "Hermite's Interpolating Polynomials\n\nHere, we implement a two point Hermite interpolation function which passes through the function and its first derivative for the interval 01. The polynomial is determined by four constraint equations, matching the function and its derivative at the points 0 and 1.\n\nStart by defining the 3rd degree polynomial and its derivative by\n\nbeginalign*\ng(x) = a_0 + a_1 x + a_2 x^2 + a_3 x^3  \ng(x) = a_1 + 2 a_2 x + 3 a_3 x^2 \nendalign*\n\nand apply the constraints\n\nbeginalign*\ng(0) = f_0   Rightarrow  a_0 = f_0  \ng(1) = f_1   Rightarrow  a_0 + a_1 + a_2 + a_3 = f_1  \ng(0) = f_0   Rightarrow  a_1 = f_0  \ng(1) = f_1   Rightarrow  a_1 + 2 a_2 + 3 a_3 = f_1  \nendalign*\n\nSolving for a_0 a_1 a_2 a_3 leads to\n\nbeginalign*\na_0 = f_0  \na_1 = f_0  \na_2 = - 3 f_0 + 3 f_1 - 2 f_0 - f_1  \na_3 = 2 f_0 - 2 f_1 + f_0 + f_1 \nendalign*\n\nso that the polynomial g(x) reads\n\ng(x) = f_0 + f_0 x + (- 3 f_0 + 3 f_1 - 2 f_0 - f_1) x^2 + (2 f_0 - 2 f_1 + f_0 + f_1) x^3 \n\nThe function and derivative values can be factored out, so that g(x) can be rewritten as\n\ng(x) = f_0 (1 - 3 x^2 + 2 x^3) + f_1 (3 x^2 - 2 x^3) + f_0 (x - 2 x^2 + x^3) + f_1 (- x^2 + x^3) \n\nor in generic form as\n\ng(x) = f_0 a_0(x) + f_1 a_1(x) + f_0 b_0(x) + f_1 b_1(x) \n\nwith basis functions\n\nbeginalign*\na_0 (x) = 1 - 3 x^2 + 2 x^3  \nb_0 (x) = x - 2 x^2 + x^3  \na_1 (x) = 3 x^2 - 2 x^3  \nb_1 (x) = - x^2 + x^3 \nendalign*\n\nThe derivative g(x) accordingly reads\n\ng(x) = f_0 a_0(x) + f_1 a_1(x) + f_0 b_0(x) + f_1 b_1(x) \n\nwith\n\nbeginalign*\na_0 (x) = - 6 x + 6 x^2  \nb_0 (x) = 1 - 4 x + 3 x^2  \na_1 (x) = 6 x - 6 x^2  \nb_1 (x) = - 2 x + 3 x^2 \nendalign*\n\nThe basis functions a_0and a_1 are associated with the function values at x_0 and x_1, respectively, while the basis functions b_0 and b_1 are associated with the derivative values at x_0 and x_1. The basis functions satisfy the following relations,\n\nbeginalign*\na_i (x_j) = delta_ij  \nb_i (x_j) = 0  \na_i (x_j) = 0  \nb_i (x_j) = delta_ij  \nij = 0 1 \nendalign*\n\nwhere delta_ij denotes the Kronecker-delta, so that\n\nbeginalign*\ng(0) = f_0  \ng(1) = f_1  \ng(0) = f_0  \ng(1) = f_1 \nendalign*\n\n\n\n"
},

{
    "location": "modules/interpolation.html#Interpolation-1",
    "page": "Interpolation",
    "title": "Interpolation",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Interpolation]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_linear.html#",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_linear.html#Linear-Solvers-1",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/linear/linear_solvers.jl\",\n           \"solvers/linear/lu_solver_lapack.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_nonlinear.html#",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_nonlinear.html#Nonlinear-Solvers-1",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/nonlinear/nonlinear_solvers.jl\",\n           \"solvers/nonlinear/jacobian.jl\",\n           \"solvers/nonlinear/abstract_fixed_point_solver.jl\",\n           \"solvers/nonlinear/fixed_point_solver.jl\",\n           \"solvers/nonlinear/abstract_newton_solver.jl\",\n           \"solvers/nonlinear/newton_solver.jl\",\n           \"solvers/nonlinear/quasi_newton_solver.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solutions.html#",
    "page": "Solutions",
    "title": "Solutions",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.PSolutionPDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.PSolutionPDAE",
    "category": "Type",
    "text": "Parallel Solution of a partitioned differential algebraic equation.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.SSolutionPDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SSolutionPDAE",
    "category": "Type",
    "text": "Serial Solution of a partitioned differential algebraic equation.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Print error for solutions of equations not implemented, yet.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for ODE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for partitioned DAE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for DAE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for implicit ODE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for partitioned ODE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.Solution",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.Solution",
    "category": "Type",
    "text": "Create solution for implicit DAE.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.SolutionDAE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionDAE",
    "category": "Type",
    "text": "Solution of a differential algebraic equation.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.SolutionODE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionODE",
    "category": "Type",
    "text": "SolutionODE: Solution of an ordinary differential equation\n\nContains all fields necessary to store the solution of an ODE.\n\nFields\n\nnd: dimension of the dynamical variable q\nnt: number of time steps to store\nni: number of initial conditions\nt:  time steps\nq:  solution q[nd, nt+1, ni] with q[:,0,:] the initial conditions\nntime: number of time steps to compute\nnsave: save every nsave'th time step\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.SolutionPODE",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.SolutionPODE",
    "category": "Type",
    "text": "Solution of a partitioned ordinary differential equation.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "Function",
    "text": "Append solution to HDF5 file.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "Function",
    "text": "Append solution to HDF5 file.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "Function",
    "text": "Append solution to HDF5 file.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "page": "Solutions",
    "title": "GeometricIntegrators.CommonFunctions.write_to_hdf5",
    "category": "Function",
    "text": "Append solution to HDF5 file.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{DT}, Tuple{GeometricIntegrators.Solutions.SolutionDAE{DT,TT,N} where N,AbstractString,Int64}, Tuple{GeometricIntegrators.Solutions.SolutionDAE{DT,TT,N} where N,AbstractString}, Tuple{TT}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "Method",
    "text": "Creates HDF5 file and initialises datasets for DAE solution object.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{DT}, Tuple{GeometricIntegrators.Solutions.SolutionODE{DT,TT,N} where N,AbstractString,Int64}, Tuple{GeometricIntegrators.Solutions.SolutionODE{DT,TT,N} where N,AbstractString}, Tuple{TT}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "Method",
    "text": "Creates HDF5 file and initialises datasets for ODE solution object.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{DT}, Tuple{GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N} where N,AbstractString,Int64}, Tuple{GeometricIntegrators.Solutions.SolutionPDAE{DT,TT,N} where N,AbstractString}, Tuple{TT}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "Method",
    "text": "Creates HDF5 file and initialises datasets for PDAE solution object.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.create_hdf5-Union{Tuple{DT}, Tuple{GeometricIntegrators.Solutions.SolutionPODE{DT,TT,N} where N,AbstractString,Int64}, Tuple{GeometricIntegrators.Solutions.SolutionPODE{DT,TT,N} where N,AbstractString}, Tuple{TT}} where TT where DT",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.create_hdf5",
    "category": "Method",
    "text": "Creates HDF5 file and initialises datasets for PODE solution object.\n\n\n\n"
},

{
    "location": "modules/solutions.html#GeometricIntegrators.Solutions.writeSolutionToHDF5-Tuple{GeometricIntegrators.Solutions.Solution,AbstractString}",
    "page": "Solutions",
    "title": "GeometricIntegrators.Solutions.writeSolutionToHDF5",
    "category": "Method",
    "text": "Creates HDF5 file, writes solution to file, and closes file.\n\n\n\n"
},

{
    "location": "modules/solutions.html#Solutions-1",
    "page": "Solutions",
    "title": "Solutions",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solutions]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/tableaus.html#",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "page",
    "text": ""
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.AbstractTableau",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.AbstractTableau",
    "category": "Type",
    "text": "Holds the information for the various methods' tableaus.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.AbstractTableauIRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.AbstractTableauIRK",
    "category": "Type",
    "text": "Holds the tableau of an implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.AbstractTableauPRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.AbstractTableauPRK",
    "category": "Type",
    "text": "Holds the tableau of a partitioned Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.AbstractTableauRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.AbstractTableauRK",
    "category": "Type",
    "text": "Holds the tableau of a Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.CoefficientsARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.CoefficientsARK",
    "category": "Type",
    "text": "Holds the coefficients of an additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.CoefficientsMRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.CoefficientsMRK",
    "category": "Type",
    "text": "Holds the multiplier Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.CoefficientsPRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.CoefficientsPRK",
    "category": "Type",
    "text": "Holds the coefficients of a projective Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.CoefficientsRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.CoefficientsRK",
    "category": "Type",
    "text": "Holds the coefficients of a Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauARK",
    "category": "Type",
    "text": "Holds the tableau of a additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauDIRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauDIRK",
    "category": "Type",
    "text": "Holds the tableau of a diagonally implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauEPRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauEPRK",
    "category": "Type",
    "text": "TableauEPRK: Tableau of an Explicit Partitioned Runge-Kutta method\n\nbeginalign*\nV_ni = hphantom- dfracpartial Hpartial p (Q_ni P_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ki = - dfracpartial Hpartial q (Q_ni P_ni)  \nP_ni = p_n + h  sum limits_i=1^s bara_ij  F_nj  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni \nendalign*\n\nusually satisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauERK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauERK",
    "category": "Type",
    "text": "Holds the tableau of an explicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauFIRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauFIRK",
    "category": "Type",
    "text": "Holds the tableau of a fully implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauGLM",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauGLM",
    "category": "Type",
    "text": "Holds the tableau of a general linear method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauIPRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauIPRK",
    "category": "Type",
    "text": "TableauIPRK: Tableau of an Implicit Partitioned Runge-Kutta method\n\nbeginalign*\nP_ni = dfracpartial Lpartial v (Q_ni V_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ki = dfracpartial Lpartial q (Q_ni V_ni)  \nP_ni = p_n + h  sum limits_i=1^s bara_ij  F_nj  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni \nendalign*\n\nusually satisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauPARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauPARK",
    "category": "Type",
    "text": "Holds the tableau of an partitioned additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauSARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauSARK",
    "category": "Type",
    "text": "Holds the tableau of a spezialized additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauSIRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauSIRK",
    "category": "Type",
    "text": "Holds the tableau of a singly implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauSPARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauSPARK",
    "category": "Type",
    "text": "Holds the tableau of a spezialized partitioned additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauVPARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauVPARK",
    "category": "Type",
    "text": "Holds the tableau of an variational partitioned additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauVPRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauVPRK",
    "category": "Type",
    "text": "TableauVPRK: Tableau of a Variational Partitioned Runge-Kutta method\n\nbeginalign*\nP_ni = dfracpartial Lpartial v (Q_ni V_ni)  \nQ_ni = q_n + h sum limits_j=1^s a_ij  V_nj  \nq_n+1 = q_n + h sum limits_i=1^s b_i  V_ni  \nF_ki = dfracpartial Lpartial q (Q_ni V_ni)  \nP_ni = p_n + h  sum limits_i=1^s bara_ij  F_nj - d_i lambda  \np_n+1 = p_n + h sum limits_i=1^s barb_i  F_ni  \n\n0 = sum limits_i=1^s d_i V_i  \nendalign*\n\nsatisfying the symplecticity conditions\n\nbeginalign*\nb_i bara_ij + b_j a_ji = b_i b_j  \nbarb_i = b_i \nendalign*\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.TableauVSPARK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.TableauVSPARK",
    "category": "Type",
    "text": "Holds the tableau of an variational special partitioned additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauERK4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK4",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (1/6 rule)\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauERK438-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK438",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (3/8 rule)\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauExplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitEuler",
    "category": "Method",
    "text": "Tableau for explicit Euler method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauExplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitMidpoint",
    "category": "Method",
    "text": "Tableau for explicit midpoint method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauGLRKpSymmetric-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRKpSymmetric",
    "category": "Method",
    "text": "Tableau for Gauss-Legendre method with s stages and symplectic projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauGLRKpSymplectic-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRKpSymplectic",
    "category": "Method",
    "text": "Tableau for Gauss-Legendre method with s stages and symplectic projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauHeun-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauHeun",
    "category": "Method",
    "text": "Tableau for Heun's method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauImplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitEuler",
    "category": "Method",
    "text": "Implicit Euler\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauImplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitMidpoint",
    "category": "Method",
    "text": "Implicit Midpoint\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauKutta-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauKutta",
    "category": "Method",
    "text": "Tableau for Kutta's method of order three\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIA4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIA4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIA Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymmetric-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymmetric",
    "category": "Method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymplectic-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB2pSymplectic",
    "category": "Method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symplectic projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymmetric-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymmetric",
    "category": "Method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymplectic-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIAIIIB3pSymplectic",
    "category": "Method",
    "text": "Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symplectic projection.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIB2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIB3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIB4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIB4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIB Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIC2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIC3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIC4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIC4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIC Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIID2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIID3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIID4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIID4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIID Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIE2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIE3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIE4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIE4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIE Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIF2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIF3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIF4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIF4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIF Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIG2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG2",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIG3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG3",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauLobIIIG4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauLobIIIG4",
    "category": "Method",
    "text": "Gauss-Lobatto-IIIG Runge-Kutta, s=4\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauRadIIA2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRadIIA2",
    "category": "Method",
    "text": "Gauss-Radau-IIA Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauRadIIA3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauRadIIA3",
    "category": "Method",
    "text": "Gauss-Radau-IIA Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauSRK3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSRK3",
    "category": "Method",
    "text": "Gauss-Legendre Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauSymplecticEulerA-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerA",
    "category": "Method",
    "text": "Tableau for symplectic Euler-A method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauSymplecticEulerB-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerB",
    "category": "Method",
    "text": "Tableau for symplectic Euler-B method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPGLRK-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPGLRK",
    "category": "Method",
    "text": "Tableau for variational Gauss-Legendre method with s stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB2",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with two stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB3",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with three stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIAIIIB4",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIA-IIIB method with four stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIC2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC2",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with two stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIC3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC3",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with three stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIC4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIC4",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIC-III method with four stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIID2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID2",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with two stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIID3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID3",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with three stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIID4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIID4",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIID method with four stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIE2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE2",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with two stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIE3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE3",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with three stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIE4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIE4",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIE method with four stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIF2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF2",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with two stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIF3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF3",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with three stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPLobIIIF4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPLobIIIF4",
    "category": "Method",
    "text": "Tableau for variational Gauss-Lobatto IIIF method with four stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauVPSRK3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauVPSRK3",
    "category": "Method",
    "text": "Tableau for variational symmetric Runge-Kutta method with 3 stages\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.readTableauERKFromFile-Tuple{AbstractString,AbstractString}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.readTableauERKFromFile",
    "category": "Method",
    "text": "Read explicit Runge-Kutta tableau from file.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.writeTableauToFile-Union{Tuple{AbstractString,GeometricIntegrators.Tableaus.AbstractTableauRK{T}}, Tuple{T}} where T",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.writeTableauToFile",
    "category": "Method",
    "text": "Write Runge-Kutta tableau to file.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.CoefficientsPGLRK",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.CoefficientsPGLRK",
    "category": "Type",
    "text": "Holds the coefficients of a projected Gauss-Legendre Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Base.show-Tuple{IO,GeometricIntegrators.Tableaus.CoefficientsARK}",
    "page": "Tableaus",
    "title": "Base.show",
    "category": "Method",
    "text": "Print additive Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Base.show-Tuple{IO,GeometricIntegrators.Tableaus.CoefficientsMRK}",
    "page": "Tableaus",
    "title": "Base.show",
    "category": "Method",
    "text": "Print multiplier Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Base.show-Tuple{IO,GeometricIntegrators.Tableaus.CoefficientsPGLRK}",
    "page": "Tableaus",
    "title": "Base.show",
    "category": "Method",
    "text": "Print Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Base.show-Tuple{IO,GeometricIntegrators.Tableaus.CoefficientsPRK}",
    "page": "Tableaus",
    "title": "Base.show",
    "category": "Method",
    "text": "Print projective Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Base.show-Tuple{IO,GeometricIntegrators.Tableaus.CoefficientsRK}",
    "page": "Tableaus",
    "title": "Base.show",
    "category": "Method",
    "text": "Print Runge-Kutta coefficients.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.readTableauRKHeaderFromFile-Tuple{Any}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.readTableauRKHeaderFromFile",
    "category": "Method",
    "text": "Reads and parses Tableau metadata from file.\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Tableaus-1",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Tableaus]\nOrder   = [:constant, :type, :macro, :function]"
},

]}
