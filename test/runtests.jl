using Rayden
using Base.Test

@testset "Components" begin 
    @testset "Intersection" begin 
        @testset "Sphere" begin 

            r = Ray(Vec(1., 2., 10.), Vec(.01, -0.12, -1.))
            s = Dome(Sphere(Vec(-1., 1., 2.), 5.), Vec(0., 0., 1.), deg2rad(90))
            intersect!(s, r)
            d = s.round.center - r.orig
            @test Rayden.norm(d) ≈ s.round.radius

            r = Ray(Vec(1., 2., 10.), Vec(.01, -0.12, -1.))
            s = Dome(Sphere(Vec(-1., 1., 2.), 5.), Vec(0., 0., 1.), deg2rad(10))
            intersect!(s, r)
            @test isinf(r.orig.x)

        end
        @testset "Ellipsoid" begin 

            r = Ray(Vec(1.1, 1.98, 10.11), Vec(-0.022, 0.2, -1.1))
            s = Dome(Ellipsoid(Vec(-1.2, 1.1, 2.4), 5.2, 4., 2.), Vec(0., 0., 1.), deg2rad(90))
            intersect!(s, r)
            d = s.round.center - r.orig
            @test Rayden.norm(Rayden.unEllipsoidise(d, s.round)) ≈ s.round.radiusx

            r = Ray(Vec(1.1, 1.98, 10.11), Vec(-0.022, 0.2, -1.1))
            s = Dome(Ellipsoid(Vec(-1.2, 1.1, 2.4), 5.2, 4., 2.), Vec(0., 0., 1.), deg2rad(40))
            intersect!(s, r)
            @test isinf(r.orig.x)

        end
    end
    @testset "Refraction" begin

            r = Ray(Vec(1., 2., 10.), Vec(.01, -0.12, -1.))
            dir_org = deepcopy(r.dir)
            s = Dome(Sphere(Vec(-1., 1., 2.), 5.), Vec(0., 0., 1.), deg2rad(90))
            intersect!(s, r)
            i = Interface(s, 1, 1.)
            rotate!(i, r)
            @test dir_org == r.dir

    end
end


