
"Holds the tableau of a stochastic explicit Runge-Kutta method."
struct TableauSERK{T} <: AbstractTableauERK{T}
    name::Symbol
    s::Int

    qdrift::CoefficientsRK{T}
    qdiff::CoefficientsRK{T}
    qdiff2::CoefficientsRK{T}

    # Order of the tableau is not included, because unlike in the deterministic
    # setting, it depends on the properties of the noise (e.g., the dimension of
    # the Wiener process and the commutativity properties of the diffusion matrix)
    #
    # Orders stored in qdrift, qdiff and qdiff2 are understood as the classical orders of these methods.

    function TableauSERK{T}(name, qdrift, qdiff, qdiff2) where {T}
        @assert qdrift.s == qdiff.s == qdiff2.s
        @assert qdrift.c[1] == qdiff.c[1] == qdiff2.c[1] == 0.
        @assert istrilstrict(qdrift.a)
        @assert istrilstrict(qdiff.a)
        @assert istrilstrict(qdiff2.a)
        @assert !(qdrift.s==1. && qdrift.a[1,1] ≠ 0.)
        @assert !(qdiff.s==1. && qdiff.a[1,1] ≠ 0.)
        @assert !(qdiff2.s==1. && qdiff2.a[1,1] ≠ 0.)
        new(name, qdrift.s, qdrift, qdiff, qdiff2)
    end
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}, qdiff2::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, qdiff2)
end

function TableauSERK(name, qdrift::CoefficientsRK{T}, qdiff::CoefficientsRK{T}) where {T}
    TableauSERK{T}(name, qdrift, qdiff, CoefficientsRK(T, :NULL, 0, zero(qdrift.a), zero(qdrift.b), zero(qdrift.c)) )
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T},
                                    order_diff2::Int, a_diff2::Matrix{T}, b_diff2::Vector{T}, c_diff2::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(name, order_diff2, a_diff2, b_diff2, c_diff2))
end

function TableauSERK(name::Symbol, order_drift::Int, a_drift::Matrix{T}, b_drift::Vector{T}, c_drift::Vector{T},
                                    order_diff::Int, a_diff::Matrix{T}, b_diff::Vector{T}, c_diff::Vector{T}) where {T}
    TableauSERK{T}(name, CoefficientsRK(name, order_drift, a_drift, b_drift, c_drift), CoefficientsRK(name, order_diff, a_diff, b_diff, c_diff),
                    CoefficientsRK(:NULL, 0, zero(a_drift), zero(b_drift), zero(c_drift)))
end




"Stochastic Explicit Runge-Kutta integrator."
struct IntegratorSERK{DT, TT, ET <: SDE{DT,TT}} <: StochasticIntegrator{DT,TT}
    equation::ET
    tableau::TableauSERK{TT}
    Δt::TT
    ΔW::Vector{DT}
    ΔZ::Vector{DT}
    Δy::Vector{DT}

    q::Vector{Vector{DT}}     # q[m]    - holds the previous time step solution (for m-th sample path)
    Q::Vector{Vector{DT}}     # Q[j][k] - the k-th component of the j-th internal stage
    V::Vector{Vector{DT}}     # V[j][k] - the k-th component of v(Q[j])
    B::Vector{Matrix{DT}}     # B[j]    - the diffusion matrix B(Q[j])


    function IntegratorSERK{DT,TT}(equation::ET, tableau, Δt::TT) where {DT, TT, ET <: SDE{DT,TT}}
        D = equation.d
        M = equation.m
        NS= max(equation.ns,equation.ni)
        S = tableau.s

        # create solution vectors
        q = create_solution_vector(DT, D, NS)

        # create internal stage vectors
        Q = create_internal_stage_vector(DT, D, S)
        V = create_internal_stage_vector(DT, D, S)
        B = create_internal_stage_vector(DT, D, M, S)

        new{DT,TT,ET}(equation, tableau, Δt, zeros(DT,M), zeros(DT,M), zeros(DT,M), q, Q, V, B)
    end
end

function IntegratorSERK(equation::SDE{DT,TT}, tableau::TableauSERK{TT}, Δt::TT) where {DT,TT}
    IntegratorSERK{DT,TT}(equation, tableau, Δt)
end

function initialize!(int::IntegratorSERK, sol::SolutionSDE, m::Int)
    check_solution_dimension_asserts(sol, m)

    # copy initial conditions from solution
    get_initial_conditions!(sol, int.q[m], m)
end

"""
Integrate SDE with explicit Runge-Kutta integrator.
  Calculating the n-th time step of the explicit integrator for the sample path m
"""
function integrate_step!(int::IntegratorSERK{DT,TT,FT}, sol::SolutionSDE{DT,TT,NQ,NW}, m::Int, n::Int) where {DT,TT,FT,NQ,NW}
    local tᵢ::TT
    local ydrift::DT

    # copy the increments of the Brownian Process
    if NW==2
        # Multidimensional Brownian motion, 1 sample path
        for l = 1:sol.nm
            int.ΔW[l] = sol.W.ΔW[l,n]
            int.ΔZ[l] = sol.W.ΔZ[l,n]
        end
    elseif NW==3
        # Multidimensional Brownian motion, m-th sample path
        for l = 1:sol.nm
            int.ΔW[l] = sol.W.ΔW[l,n,m]
            int.ΔZ[l] = sol.W.ΔZ[l,n,m]
        end
    end

    # calculates v(t,tQ) and assigns to the i-th column of V
    int.equation.v(sol.t[0] + (n-1)*int.Δt, int.q[m], int.V[1])
    # calculates B(t,tQ) and assigns to the matrix BQ[1][:,:]
    int.equation.B(sol.t[0] + (n-1)*int.Δt, int.q[m], int.B[1])


    @inbounds for i in 2:int.tableau.s
        for k in eachindex(int.Q[i])
            # contribution from the drift part
            ydrift = 0.
            for j = 1:i-1
                ydrift += int.tableau.qdrift.a[i,j] * int.V[j][k]
            end

            # ΔW contribution from the diffusion part
            int.Δy .= 0.
            for j = 1:i-1
                for l = 1:sol.nm
                    int.Δy[l] += int.tableau.qdiff.a[i,j] * int.B[j][k,l]
                end
            end

            int.Q[i][k] = int.q[m][k] + int.Δt * ydrift + dot(int.Δy,int.ΔW)

            # ΔZ contribution from the diffusion part
            if int.tableau.qdiff2.name ≠ :NULL
                int.Δy .= 0.
                for j = 1:i-1
                    for l = 1:sol.nm
                        int.Δy[l] += int.tableau.qdiff2.a[i,j] * int.B[j][k,l]
                    end
                end

                int.Q[i][k] += dot(int.Δy,int.ΔZ)/int.Δt
            end

        end
        tᵢ = sol.t[0] + (n-1)*int.Δt + int.Δt * int.tableau.qdrift.c[i]
        int.equation.v(tᵢ, int.Q[i], int.V[i])
        int.equation.B(tᵢ, int.Q[i], int.B[i])
    end

    # compute final update
    if int.tableau.qdiff2.name == :NULL
        update_solution!(int.q[m], int.V, int.B, int.tableau.qdrift.b, int.tableau.qdiff.b, int.Δt, int.ΔW, int.Δy)
    else
        update_solution!(int.q[m], int.V, int.B, int.tableau.qdrift.b, int.tableau.qdiff.b, int.tableau.qdiff2.b, int.Δt, int.ΔW, int.ΔZ, int.Δy)
    end

    # take care of periodic solutions
    cut_periodic_solution!(int.q[m], int.equation.periodicity)

    # copy to solution
    set_solution!(sol, int.q[m], n, m)
end
