
using Polynomials: polyint
using GeometricIntegrators.BasisFunctions: legendre_polynomial

@testset "$(rpad("Legendre Basis Tests",80))" begin

    b = LegendreBasis(Float64, 5)

    @test evaluate(b, 1, 0.0) == 1.0
    @test evaluate(b, 1, 0.5) == 1.0
    @test evaluate(b, 1, 1.0) == 1.0

    @test evaluate(b, 2, 0.5) == 0.0
    @test evaluate(b, 3, 0.5) != 0.0
    @test evaluate(b, 4, 0.5) == 0.0
    @test evaluate(b, 5, 0.5) != 0.0

    @test derivative(b, 1, 0.0) == 0.0
    @test derivative(b, 1, 0.5) == 0.0
    @test derivative(b, 1, 1.0) == 0.0

    @test derivative(b, 2, 0.0) == derivative(b, 2, 0.5) == derivative(b, 2, 1.0)

    @test integral(b, 1, 0.0) == 0.0
    @test integral(b, 1, 0.5) == 0.5
    @test integral(b, 1, 1.0) == 1.0

    @test integral(b, 2, 0.0) == 0.0
    @test integral(b, 2, 1.0) == 0.0


    for q in 1:10, p in 1:10
        f1, p1 = legendre_polynomial(q)
        f2, p2 = legendre_polynomial(p)
        int = polyint(p1 * p2)
        @test f1 * f2 * (int(1.0) - int(0.0)) ≈ (q == p ? 1.0 : 0.0) atol=1e-3
    end

end
