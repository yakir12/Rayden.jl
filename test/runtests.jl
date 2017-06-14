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
            i = Interface(s, false, 1.)
            rotate!(i, r)
            @test dir_org == r.dir

            r = Ray(Vec(1., 2., 10.), Vec(.01, -0.12, -1.))
            dir_org = deepcopy(r.dir)
            s = Dome(Ellipsoid(Vec(-1.2, 1.1, 2.4), 5.2, 4., 2.), Vec(0., 0., 1.), deg2rad(90))
            intersect!(s, r)
            i = Interface(s, false, 1.)
            rotate!(i, r)
            @test dir_org == r.dir

    end
    @testset "Eye" begin


        function coordinates(d::Dome{Sphere})
    s = d.round
    n = 25
    x = zeros(n, n)
    y = zeros(n, n)
    z = zeros(n, n)
    for (i, θ) in enumerate(linspace(0, 2pi, 25)), (j, ψ) in enumerate(linspace(0, d.open, 25))
        sinψ = sin(ψ)
        x[i, j] = s.center.x + s.radius*cos(θ)*sinψ
        y[i, j] = s.center.y + s.radius*sin(θ)*sinψ
        z[i, j] = s.center.z + d.dir.z*s.radius*cos(ψ)
    end
    return (x, y, z)
end
        function coordinates(d::Dome{Ellipsoid})
    s = d.round
    n = 25
    x = zeros(n, n)
    y = zeros(n, n)
    z = zeros(n, n)
    for (i, θ) in enumerate(linspace(0, 2pi, 25)), (j, ψ) in enumerate(linspace(0, d.open, 25))
        sinψ = sin(ψ)
        x[i, j] = s.center.x + s.radiusx*cos(θ)*sinψ
        y[i, j] = s.center.y + s.radiusy*sin(θ)*sinψ
        z[i, j] = s.center.z + d.dir.z*s.radiusz*cos(ψ)
    end
    return (x, y, z)
end



        using Rayden, Plots
        plotlyjs()
        x = Float64[]
        y = Float64[]
        z = Float64[]
        ray = Ray(Vec(0., 0., 100.), Vec(.0, 0., -1.))
        ss = scallop()
        map(x -> surface!(coordinates(x.dome)), ss)
        for i in ss[1:5]
            push!(x, ray.orig.x)
            push!(y, ray.orig.y)
            push!(z, ray.orig.z)
            intersect!(i, ray)
            if i.n != 1
                Rayden.rotate!(i, ray)
            end
        end
        plot!(x,y,z)



    end
end


