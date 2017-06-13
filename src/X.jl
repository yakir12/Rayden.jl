module X

using Z, Plots
#plotlyjs()

import Plots.plot

export fwhm_aperture_distance, fwhm_aperture, plot

getfwhm2(a::Dict{String,Interface}, aperture::Float64, L::Float64) = getfwhm(a, Light(L, aperture, a))

function fwhm_aperture(a::Dict{String,Interface}, aperturelim::Tuple{Float64, Float64}, L::Float64)
    na = 10
    apertures = linspace(aperturelim...,na)
    y = Dict("proximal_retina" => zeros(na), "distal_retina" => zeros(na))
    for (i, ap) in enumerate(apertures), (k,v) in getfwhm2(a, ap, L)
        y[k][i] = rad2deg(v)
    end
    p = plot(xlabel="Aperture (microns)", ylabel = "FWHM (degrees)")
    for (k,v) in y
        plot!(apertures, v, label = k)
    end
    plot!()
    return p
end

function fwhm_aperture_distance(a::Dict{String,Interface}, aperturelim::Tuple{Float64, Float64}, distancelim::Tuple{Float64, Float64})
    na = 10
    apertures = linspace(aperturelim...,na)
    distances = linspace(distancelim...,na)
    y = Dict("proximal_retina" => zeros(na, na), "distal_retina" => zeros(na, na))
    for (i, ds) in enumerate(distances), (j, ap) in enumerate(apertures), (k,v) in getfwhm2(a, ap, ds)
        y[k][i,j] = rad2deg(v)
    end
    p = Any[]
    for (k,v) in y
        push!(p, contour(apertures, distances*1e-3, log(v), fill=true, yscale = :log10, title = k, xlabel="Aperture (microns)", ylabel = "Distance (mm)"))
    end
    pp = plot(p...)
    return pp
end

##Plotting

function plot(c::Circle)
    n = 50
    x = zeros(n)
    y = zeros(n)
    for (i, theta) in enumerate(linspace(c.cc.interval..., n))
        xy = c.cc.r*exp(1im*theta)
        x[i] = real(xy)
        y[i] = imag(xy) + c.h
    end
    return (x,y)
end

function plot(el::Ellipse)
    x, y = plot(el.c)
    return (x, y/el.ratio)
end

function plot(s::Side)
    p1 = s.o
    p2 = s.o + s.l*s.d
    x = [p1[1], p2[1]]
    y = [p1[2], p2[2]]
    return (x, y)
end

plot(i::Interface) = plot(i.l)

function plot(a::Dict{String, Interface}, l::Light, n::Int)
    r = Ray(l, rand())
    b = collect(linspace(0,1,n))
    pop!(b)
    p = [(Float64[], Float64[])]
    for i in b
        x = Float64[]
        y = Float64[]
        Ray(r, l, i)
        push!(x, r.o[1])
        push!(y, r.o[2])
        trace!(r, a, x, y)
        push!(p, (x, y))
        push!(p, (-x, y))
    end
    return p
end


function plot(a::Dict{String,Interface}, l::Light)
    buff = 10
    ylim = (a["mirror"].l.c.h/a["mirror"].l.ratio - a["mirror"].l.r2 - buff,buff)
    xlim = (-a["mirror"].l.c.cc.r - buff, a["mirror"].l.c.cc.r + buff)

    p = plot(aspect_ratio = :equal, legend = false, xlims = xlim, ylim = ylim)
    for i in values(a)
        x, y = plot(i)
        plot!(x, y, linecolor=:black)
    end
    ps = plot(a, l, 10)
    for (x,y) in ps
        plot!(x, y, linecolor=:red)
    end
    plot!()
    return p
end

end
