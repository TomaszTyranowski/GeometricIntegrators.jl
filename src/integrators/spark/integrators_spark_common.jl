
function create_integrator_cache(int::AbstractIntegratorSPARK{DT,TT}) where {DT,TT}
    IntegratorCacheSPARK{DT, TT, ndims(int), nstages(int), pstages(int)}()
end


function initialize!(int::AbstractIntegratorSPARK, cache::IntegratorCacheSPARK)
    cache.t̅ = cache.t - timestep(int)

    equation(int).v(cache.t, cache.q, cache.p, cache.v)
    equation(int).f(cache.t, cache.q, cache.p, cache.f)

    initialize!(int.iguess, cache.t, cache.q, cache.p, cache.v, cache.f,
                            cache.t̅, cache.q̅, cache.p̅, cache.v̅, cache.f̅)
end


function update_solution!(int::AbstractIntegratorSPARK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # compute final update
    update_solution!(sol.q, sol.q̃, int.cache.Vi, int.params.tab.q.b, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Fi, int.params.tab.p.b, timestep(int))

    # compute projection
    update_solution!(sol.q, sol.q̃, int.cache.Up, int.params.tab.q.β, timestep(int))
    update_solution!(sol.p, sol.p̃, int.cache.Gp, int.params.tab.p.β, timestep(int))
    # TODO # update_multiplier!(sol.λ, int.cache.Λp, int.params.tab.λ.b)
end


"Integrate an implicit DAE with a specialised partitioned additive Runge-Kutta integrator."
function integrate_step!(int::AbstractIntegratorSPARK{DT,TT}, sol::AtomicSolutionPDAE{DT,TT}) where {DT,TT}
    # update nonlinear solver parameters from cache
    update_params!(int.params, sol)

    # compute initial guess
    initial_guess!(int, sol)

    # reset cache
    reset!(sol, timestep(int))

    # call nonlinear solver
    solve!(int.solver)

    # check_jacobian(int.solver)
    # print_jacobian(int.solver)

    # print solver status
    print_solver_status(int.solver.status, int.solver.params)

    # check if solution contains NaNs or error bounds are violated
    check_solver_status(int.solver.status, int.solver.params)

    # compute vector fields at internal stages
    compute_stages!(int.solver.x, int.cache, int.params)

    # compute final update
    update_solution!(int, sol)

    # copy solution to initial guess
    update!(int.iguess, sol.t, sol.q, sol.p, sol.v, sol.f)
end
