module LotkaVolterra2dProblem

    using GeometricIntegrators.Equations

    export lotka_volterra_2d_ode,  lotka_volterra_2d_pdae,
           lotka_volterra_2d_iode, lotka_volterra_2d_idae,
           lotka_volterra_2d_vode, lotka_volterra_2d_vdae,
           lotka_volterra_2d_dg, lotka_volterra_2d_vdae_slrk


    Δt = 0.01
    nt = 10

    const A1=1.0
    const A2=1.0
    const B1=1.0
    const B2=2.0

    const X0=1.0
    const Y0=1.0


    function ϑ₁(t, q)
        q[2] + log(q[2]) / q[1]
    end

    function ϑ₂(t, q)
        q[1]
    end


    function v₁(t, q)
        q[1] * (A2*q[2] - B2)
    end

    function v₂(t, q)
        q[2] * (B1 - A1*q[1])
    end


    const q₀=[X0, Y0]
    const v₀=[v₁(0, q₀), v₂(0, q₀)]


    function dϑ₁d₁(t, q)
        - log(q[2]) / q[1]^2
    end

    function dϑ₁d₂(t, q)
        1 + 1 / (q[1] * q[2])
    end

    function dϑ₂d₁(t, q)
        one(eltype(q))
    end

    function dϑ₂d₂(t, q)
        zero(eltype(q))
    end


    function d²ϑ₁d₁d₁(t, q)
        + 2 * log(q[2]) / q[1]^3
    end

    function d²ϑ₁d₁d₂(t, q)
        - 1 / (q[1]^2 * q[2])
    end

    function d²ϑ₁d₂d₁(t, q)
        - 1 / (q[1]^2 * q[2])
    end

    function d²ϑ₁d₂d₂(t, q)
        - 1 / (q[1] * q[2]^2)
    end

    function d²ϑ₂d₁d₁(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₁d₂(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₂d₁(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₂d₂(t, q)
        zero(eltype(q))
    end


    function ϑ(t::Number, q::AbstractVector, Θ::AbstractVector)
        Θ[1] = ϑ₁(t,q)
        Θ[2] = ϑ₂(t,q)
        nothing
    end

    function ϑ(t::Number, q::AbstractVector)
        [ϑ₁(t,q), ϑ₂(t,q)]
    end

    function ω(t, q, Ω)
        Ω[1,1] = 0
        Ω[1,2] = dϑ₁d₂(t,q) - dϑ₂d₁(t,q)

        Ω[2,1] = dϑ₂d₁(t,q) - dϑ₁d₂(t,q)
        Ω[2,2] = 0

        nothing
    end


    function f₁(t, q, v)
        dϑ₁d₁(t,q) * v[1] + dϑ₂d₁(t,q) * v[2]
    end

    function f₂(t, q, v)
        dϑ₁d₂(t,q) * v[1] + dϑ₂d₂(t,q) * v[2]
    end


    function g₁(t, q, v)
        dϑ₁d₁(t,q) * v[1] + dϑ₁d₂(t,q) * v[2]
    end

    function g₂(t, q, v)
        dϑ₂d₁(t,q) * v[1] + dϑ₂d₂(t,q) * v[2]
    end


    function g̅₁(t, q, v)
        d²ϑ₁d₁d₁(t,q) * q[1] * v[1] + d²ϑ₁d₂d₁(t,q) * q[1] * v[2] + d²ϑ₂d₁d₁(t,q) * q[2] * v[1] + d²ϑ₂d₂d₁(t,q) * q[2] * v[2]
    end

    function g̅₂(t, q, v)
        d²ϑ₁d₁d₂(t,q) * q[1] * v[1] + d²ϑ₁d₂d₂(t,q) * q[1] * v[2] + d²ϑ₂d₁d₂(t,q) * q[2] * v[1] + d²ϑ₂d₂d₂(t,q) * q[2] * v[2]
    end


    function hamiltonian(t, q)
        A1*q[1] + A2*q[2] - B1*log(q[1]) - B2*log(q[2])
    end


    function dHd₁(t, q)
        A1 - B1 / q[1]
    end

    function dHd₂(t, q)
        A2 - B2 / q[2]
    end

    function lotka_volterra_2d_dH(t, q, dH)
        dH[1] = dHd₁(t, q)
        dH[2] = dHd₂(t, q)
        nothing
    end


    function lotka_volterra_2d_ϑ(t, q, Θ)
        ϑ(t, q, Θ)
    end

    function lotka_volterra_2d_ϑ(t, q, v, Θ)
        ϑ(t, q, Θ)
    end

    function lotka_volterra_2d_ϑ(κ, t, q, v, Θ)
        Θ[1] = (1-κ) * ϑ₁(t,q) - κ * f₁(t,q,q)
        Θ[2] = (1-κ) * ϑ₂(t,q) - κ * f₂(t,q,q)
        nothing
    end


    function lotka_volterra_2d_ω(t, q, Ω)
        ω(t, q, Ω)
    end


    function lotka_volterra_2d_v(t, q, v)
        v[1] = v₁(t, q)
        v[2] = v₂(t, q)
        nothing
    end

    function lotka_volterra_2d_v(t, q, p, v)
        lotka_volterra_2d_v(t, q, v)
    end

    function lotka_volterra_2d_v_ham(t, q, p, v)
        v[1] = 0
        v[2] = 0
        nothing
    end

    function lotka_volterra_2d_f(t::Real, q::Vector, v::Vector, f::Vector)
        f[1] = f₁(t,q,v) - dHd₁(t,q)
        f[2] = f₂(t,q,v) - dHd₂(t,q)
        nothing
    end

    function lotka_volterra_2d_f(κ::Real, t::Real, q::Vector, v::Vector, f::Vector)
        f[1] = (1-κ) * f₁(t,q,v) - κ * (g₁(t,q,v) + g̅₁(t,q,v)) - dHd₁(t,q)
        f[2] = (1-κ) * f₂(t,q,v) - κ * (g₂(t,q,v) + g̅₂(t,q,v)) - dHd₂(t,q)
        nothing
    end

    function lotka_volterra_2d_f_ham(t::Real, q::Vector, f::Vector)
        f[1] = - dHd₁(t,q)
        f[2] = - dHd₂(t,q)
        nothing
    end

    function lotka_volterra_2d_f_ham(t::Real, q::Vector, v::Vector, f::Vector)
        lotka_volterra_2d_f_ham(t, q, f)
    end

    function lotka_volterra_2d_g(t::Real, q::Vector, v::Vector, g::Vector)
        g[1] = f₁(t,q,v)
        g[2] = f₂(t,q,v)
        nothing
    end

    function lotka_volterra_2d_g(t::Real, q::Vector, p::Vector, v::Vector, g::Vector)
        lotka_volterra_2d_g(t, q, v, g)
    end

    function lotka_volterra_2d_g(κ::Real, t::Real, q::Vector, v::Vector, g::Vector)
        g[1] = (1-κ) * f₁(t,q,v) - κ * (g₁(t,q,v) + g̅₁(t,q,v))
        g[2] = (1-κ) * f₂(t,q,v) - κ * (g₂(t,q,v) + g̅₂(t,q,v))
        nothing
    end

    # function lotka_volterra_2d_g(κ::Real, t::Real, q::Vector, v::Vector, g::Vector)
    #     g[1] = (1-κ) * g₁(t,q,v) - κ * g̅₁(t,q,v) - κ * f₁(t,q,v)
    #     g[2] = (1-κ) * g₂(t,q,v) - κ * g̅₂(t,q,v) - κ * f₂(t,q,v)
    #     nothing
    # end

    function lotka_volterra_2d_g̅(t::Real, q::Vector, v::Vector, g::Vector)
        g[1] = g₁(t,q,v)
        g[2] = g₂(t,q,v)
        nothing
    end

    function lotka_volterra_2d_g̅(t::Real, q::Vector, p::Vector, v::Vector, g::Vector)
        lotka_volterra_2d_g̅(t, q, v, g)
    end

    function lotka_volterra_2d_u(t, q, v, u)
        u[1] = v[1]
        u[2] = v[2]
        nothing
    end

    function lotka_volterra_2d_u(t, q, p, v, u)
        lotka_volterra_2d_u(t, q, v, u)
    end


    function lotka_volterra_2d_ϕ(t, q, p, ϕ)
        ϕ[1] = p[1] - ϑ₁(t,q)
        ϕ[2] = p[2] - ϑ₂(t,q)
        nothing
    end

    function lotka_volterra_2d_ψ(t, q, p, v, f, ψ)
        ψ[1] = f[1] - g₁(t,q,v)
        ψ[2] = f[2] - g₂(t,q,v)
        nothing
    end


    function lotka_volterra_2d_ode(q₀=q₀)
        ODE(lotka_volterra_2d_v, q₀; h=hamiltonian)
    end


    function lotka_volterra_2d_iode(q₀=q₀, p₀=ϑ(0, q₀))
        IODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
             lotka_volterra_2d_g, q₀, p₀;
             h=hamiltonian, v=lotka_volterra_2d_v)
    end

    function lotka_volterra_2d_vode(q₀=q₀, p₀=ϑ(0, q₀))
        VODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
             lotka_volterra_2d_g, q₀, p₀;
             v=lotka_volterra_2d_v, Ω=lotka_volterra_2d_ω, ∇H=lotka_volterra_2d_dH)
    end

    function lotka_volterra_2d_idae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀))
        IDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
             lotka_volterra_2d_u, lotka_volterra_2d_g,
             lotka_volterra_2d_ϕ, q₀, p₀, λ₀;
             v=lotka_volterra_2d_v)
    end

    function lotka_volterra_2d_pdae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀))
        PDAE(lotka_volterra_2d_v_ham, lotka_volterra_2d_f_ham,
             lotka_volterra_2d_u, lotka_volterra_2d_g,
             lotka_volterra_2d_ϕ, q₀, p₀, λ₀)
    end

    function lotka_volterra_2d_vdae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀))
        VDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
             lotka_volterra_2d_g, lotka_volterra_2d_g̅,
             lotka_volterra_2d_ϕ, lotka_volterra_2d_ψ,
             q₀, p₀, λ₀; v=lotka_volterra_2d_v)
    end

    function lotka_volterra_2d_vdae_slrk(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀))
        VDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
             lotka_volterra_2d_g, lotka_volterra_2d_g̅,
             lotka_volterra_2d_ϕ, lotka_volterra_2d_ψ,
             q₀, p₀, λ₀; v=lotka_volterra_2d_v)
    end

    function lotka_volterra_2d_dg(q₀=q₀; κ=0)
        lotka_volterra_2d_ϑ_κ(t, q, v, p) = lotka_volterra_2d_ϑ(κ, t, q, v, p)
        lotka_volterra_2d_f_κ(t, q, v, f) = lotka_volterra_2d_f(κ, t, q, v, f)
        lotka_volterra_2d_g_κ(t, q, λ, g) = lotka_volterra_2d_g(κ, t, q, λ, g)

        IODE(lotka_volterra_2d_ϑ_κ, lotka_volterra_2d_f_κ,
             lotka_volterra_2d_g_κ, q₀, p₀;
             v=lotka_volterra_2d_v)
    end

end
