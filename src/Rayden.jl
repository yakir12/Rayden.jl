module Rayden

import Base:intersect!, ⋅

export Vec, Ray, Hit, Sphere, Ellipsoid, Dome, intersect!, Interface, rotate!

#const FWHM = 2*sqrt(2*log(2))

immutable Vec
    x::Float64
    y::Float64
    z::Float64
end
import Base: +, -, *, /
+(a::Vec, b::Vec) = Vec(a.x+b.x, a.y+b.y, a.z+b.z)
-(a::Vec) = Vec(-a.x, -a.y, -a.z)
-(a::Vec, b::Vec) = Vec(a.x-b.x, a.y-b.y, a.z-b.z)
*(a::Float64, b::Vec) = Vec(a*b.x, a*b.y, a*b.z)
*(a::Int, b::Vec) = Vec(a*b.x, a*b.y, a*b.z)
*(a::Vec, b::Float64) = *(b,a)
/(a::Vec, b::Vec) = Vec(a.x/b.x, a.y/b.y, a.z/b.z)
⋅(a::Vec, b::Vec) = (a.x*b.x + a.y*b.y + a.z*b.z)
norm(a::Vec) = sqrt(a⋅a)
unitize(a::Vec) = (1. / norm(a) * a)

type Ray
    orig::Vec
    dir::Vec
    Ray(o, d) = new(o, unitize(d))
end

type Hit
    lambda::Float64
    normal::Vec
end

abstract Round

immutable Sphere <: Round
    center::Vec
    radius::Float64
    radius2::Float64
    Sphere(c, r) = new(c, r, r*r)
end

function distance(s::Sphere, ray::Ray)
    v = s.center - ray.orig
    b = v⋅ray.dir
    disc = b^2 - v⋅v + s.radius2
    if disc >= 0
        d = sqrt(disc)
        t2 = b + d
        if t2 >= 0
            t1 = b - d
            return t1 > 0 ? t1 : t2
        end
    end
    return Inf
end

immutable Ellipsoid <: Round
    center::Vec
    radiusx::Float64
    radiusy::Float64
    radiusz::Float64
    radiusx2::Float64
    radiusy2::Float64
    radiusz2::Float64
    yfactor::Float64
    zfactor::Float64
    sphere::Sphere
    function Ellipsoid(c, rx, ry, rz)
        yf = rx/ry
        zf = rx/rz
        cc = Vec(c.x, yf*c.y, zf*c.z)
        new(c, rx, ry, rz, rx^2, ry^2, rz^2, yf, zf, Sphere(cc, rx))
    end
end

Ellipsoidise(v::Vec, s::Ellipsoid) = Vec(v.x, v.y/s.yfactor, v.z/s.zfactor)
unEllipsoidise(v::Vec, s::Ellipsoid) = Vec(v.x, v.y*s.yfactor, v.z*s.zfactor)
unEllipsoidise(ray::Ray, s::Ellipsoid) = Ray(unEllipsoidise(ray.orig, s), unEllipsoidise(ray.dir, s))

function advance!(ray::Ray, l::Float64) 
    ray.orig += l*ray.dir 
end

function distance(s::Ellipsoid, ray::Ray)
    r = unEllipsoidise(ray, s)
    l = distance(s.sphere, r)
    isinf(l) && return Inf
    advance!(r, l)
    return norm(ray.orig - Ellipsoidise(r.orig, s))
end

#=function intersect{R <: Round}(s::R, ray::Ray)
    l = distance(s,  ray)
    return advance(ray, l)
end=#

immutable Dome{R <: Round}
    round::R
    dir::Vec # the center (normalized) direction of the dome 
    open::Float64 # acos(α) where α is the opening angle of the dome
    Dome(r, d, o) = new(r, unitize(d), o)
end
Dome{R <: Round}(r::R, d::Vec, o::Float64) = Dome{R}(r, d, o)

function infRay!(r::Ray)
    r.orig = Vec(Inf, Inf, Inf)
end

function intersect!(d::Dome, ray::Ray)
    l = distance(d.round, ray)
    isinf(l) && infRay!(ray)
    advance!(ray, l)
    acosα = d.dir⋅unitize(ray.orig)
    acosα > d.open && infRay!(ray)
end

immutable Interface{R <: Round}
    dome::Dome{R}
    inout::Int
    n::Float64
    n2::Float64
    Interface(d, io, n) = new(d, io, n, n*n)
end
Interface{R <: Round}(d::Dome{R}, i::Int, n::Float64) = Interface{R}(d, i, n)

normal(i::Interface{Sphere}, ray::Ray) = i.inout*unitize(ray.orig - i.dome.round.center)
normal(i::Interface{Ellipsoid}, ray::Ray) = i.inout*unitize((ray.orig - i.dome.round.center)/Vec(i.dome.round.radiusx2, i.dome.round.radiusy2, i.dome.round.radiusz2))

function rotate!(i::Interface, ray::Ray)
    N = normal(i, ray)
    a = -ray.dir⋅N
    b = i.n2*(1 - a^2)
    if b <= 1
        c = 1 - b 
        ray.dir = unitize(i.n*ray.dir + (i.n*a - sqrt(c))*N)
    else
        ray.dir = unitize(ray.dir + 2a*N)
    end
end

end # module
#=
function intersect(s::Sphere, i::Hit, ray::Ray)
    l = ray_scene(s, ray)
    if l >= i.lambda
        return i
    else
        n = ray.orig + l * ray.dir - s.center
        return Hit(l, unitize(n))
    end
end



function coordinates(s::Sphere)
    n = 25
    x = zeros(n, n)
    y = zeros(n, n)
    z = zeros(n, n)
    for (i, θ) in enumerate(linspace(0, 2pi, 25)), (j, ψ) in enumerate(linspace(0, pi, 25))
        sinψ = sin(ψ)
        x[i, j] = s.center.x + s.radius*cos(θ)*sinψ
        y[i, j] = s.center.y + s.radius*sin(θ)*sinψ
        z[i, j] = s.center.z + s.radius*cos(ψ)
    end
    return (x, y, z)
end







immutable Dome <: Scene
    sphere::Sphere
    dir::Vec # the center (normalized) direction of the dome 
    open::Float64 # acos(α) where α is the opening angle of the dome
end






function intersect(d::Dome, i::Hit, ray::Ray)
    l_ = ray_scene(d.sphere, ray)
    p2 = ray.orig + l_ * ray.dir
    acosα = dot(d.dir, unitize(p2))
    l = acosα <= d.open ? l_ : Inf
    if l >= i.lambda
        return i
    else
        n = ray.orig + l * ray.dir - s.center
        return Hit(l, unitize(n))
    end
end








immutable Ray
    o::Vector{Float64}
    d::Vector{Float64}
    Ray(o,d) = new(o, normalize(d))
end

##Shape types

abstract Line

immutable CenterCircle <: Line
    r::Float64
    interval::Tuple{Float64, Float64}
    r2::Float64
    function CenterCircle(r, t)
        interval = t > 0 ? (t, pi - t) : t < 0 ? (-pi - t, t) : (0, 2pi)
        new(r, interval, r^2)
    end
end

immutable Circle <: Line
    h::Float64
    cc::CenterCircle
    Circle(h, r, t) = new(h, CenterCircle(r, t))
end

immutable Ellipse <: Line
    h::Float64
    r1::Float64
    r2::Float64
    ratio::Float64
    c::Circle
    function Ellipse(h, r1, r2, t)
        ratio = r1/r2
        r = r1*r2/sqrt((r2*cos(t))^2 + (r1*sin(t))^2)
        y = r*sin(t)
        t2 = asin(y*ratio/r1)
        new(h, r1, r2, ratio, Circle(h*ratio, r1, t2))
    end
end

immutable Side <: Line
    o::Vector{Float64}
    d::Vector{Float64}
    l::Float64
    Side(o, d, l) = new(o, normalize(d), l)
end

immutable Direction
    n::Float64
    n2::Float64
    normal::Function
end

immutable Interface{T <: Line}
    l::T
    d::Tuple{Direction, Direction}
end

const cw = [0 1; -1 0]
const ccw = [0 -1; 1 0]

function Interface(l::Side, ri1, ri2)
    n = ri1/ri2 #here, ri1 is above the line and ri2 is below
    wc, wcc = l.d[1] > 0 ? (ccw, cw) : (cw, ccw)
    normal = x::Vector{Float64} -> wc*l.d #counter-clockwise
    d1 = Direction(n, n^2, normal)
    normal = x::Vector{Float64} -> wcc*l.d #clockwise
    d2 = Direction(1/n, 1/n^2, normal)
    Interface{Side}(l, (d1, d2))
end

function Interface(l::Circle, ri1, ri2)
    n = ri1/ri2 #here, ri1 is inside the circle and ri2 is outside
    normal = x::Vector{Float64} -> normalize([x[1], x[2] - l.h])
    d1 = Direction(n, n^2, normal)
    normal = x::Vector{Float64} -> normalize([-x[1], l.h - x[2]])
    d2 = Direction(1/n, 1/n^2, normal)
    Interface{Circle}(l, (d1, d2))
end

function Interface(l::Ellipse, ri1, ri2)
    n = ri1/ri2 #here, ri1 is inside the circle and ri2 is outside
    normal = x::Vector{Float64} -> normalize([x[1]/l.ratio, (x[2] - l.h)*l.ratio])
    d1 = Direction(n, n^2, normal)
    normal = x::Vector{Float64} -> normalize([-x[1]/l.ratio, (l.h - x[2])*l.ratio])
    d2 = Direction(1/n, 1/n^2, normal)
    Interface{Ellipse}(l, (d1, d2))
end

##Create the light

gety(x::Float64, s::Circle) = sqrt(s.cc.r2 - x^2) + s.h
gety(x::Float64, s::Ellipse) = s.c.cc.r2*sqrt(1 - (x/s.r1)^2) + s.h
gety(x::Float64, s::Side) = x*s.d[2]/s.d[1]
gety(x::Float64, a::Dict{String, Interface}) = x > a["cornea_right"].l.o[1] ? gety(x, a["cornea_right"].l) : gety(x, a["cornea_round"].l)

immutable Light
    alpha::Float64
    cosa::Float64
    o::Vector{Float64}
    function Light(l, a, s)
        r = a/2
        h = abs(gety(r, s))
        alpha = atan(r/(l + h))
        cosa = cos(alpha)
        L = (l + h)/cosa
        orig = [0, l]
        new(alpha, cosa, orig)
    end
end

function Ray(r::Ray, l::Light, b::Float64)
    z = b*(l.cosa - 1) - l.cosa
    theta = -pi/2 - atan(sqrt(-z^2 + 1)/z)
    r.o .= l.o
    r.d[1] = cos(theta)
    r.d[2] = sin(theta)
end

function Ray(l::Light, b::Float64)
    r = Ray([0.,0.],[0.,0.])
    Ray(r, l, b)
    return r
end

##Distance functions

angleinterval(p::Vector{Float64}, interval::Tuple{Float64, Float64}) = interval[1] <= atan2(p[2], p[1]) <= interval[2]

function distance(r::Ray, cc::CenterCircle)
    v = -r.o
    b =  v⋅r.d
    disc = b^2 - v⋅v + cc.r2
    if disc == 0
        p = r.o + b*r.d
        angleinterval(p, cc.interval) && return b::Float64
    elseif disc > 0
        d = sqrt(disc)
        if b > d + 1e-6
            m = Inf
            for l in [b - d, b + d]
                p = r.o + l*r.d
                if angleinterval(p, cc.interval)
                    m = min(m, l)
                end
            end
            return m
        elseif b > -d + 1e-6
            l = b + d
            p = r.o + l*r.d
            angleinterval(p, cc.interval) && return l::Float64
        end
    end
    return Inf
end

function distance(r::Ray, c::Circle)
    r.o[2] -= c.h
    l = distance(r, c.cc)
    r.o[2] += c.h
    return l
end

function distance(r::Ray, el::Ellipse)
    r1 = deepcopy(r)
    r1.o[2] *= el.ratio
    r1.d[2] *= el.ratio
    normalize!(r1.d)
    l = distance(r1, el.c)
    p2 = r1.o + l*r1.d
    p2[2] /= el.ratio
    return norm(r.o - p2)
end

cross2d(x::Vector{Float64}, y::Vector{Float64}) = x[1]*y[2] - x[2]*y[1]

function distance(r::Ray, s::Side)
    rs = cross2d(r.d, s.d)
    rs == 0 && return Inf
    qp = s.o - r.o
    f = cross2d(qp, r.d/rs)
    if 0 < f <= s.l
        d = cross2d(qp, s.d/rs)
        d > 1e-6 && return d
    end
    return Inf
end

function distance(r::Ray, a::Dict{String, Interface}, excludeme::String)
    k = ""
    m = Inf
    for (ki, face) in a
        d = distance(r, face.l)
        ki == excludeme && d < 1 && continue
        if d < m
            m = d
            k = ki
        end
    end
    return (m, k)
end

function advance!(r::Ray, l::Float64)
    r.o .+= l.*r.d
end

function rotate!{T <: Line}(r::Ray, face::Interface{T})
    i = 1
    N = face.d[i].normal(r.o)
    a = -r.d⋅N
    if a < 0
        i = 2
        N = face.d[i].normal(r.o)
        a = -r.d⋅N
    end
    b = face.d[i].n2*(1 - a*a)
    if b <= 1
        c = 1 - b 
        r.d .= face.d[i].n*r.d + (face.d[i].n*a - sqrt(c))*N
    else
        r.d .+= 2a*N
    end
    normalize!(r.d)
end

function hit!(r::Ray, a::Dict{String, Interface}, excludeme::String)
    l, k = distance(r, a, excludeme)
    advance!(r, l)
    rotate!(r, a[k])
    return k
end

function getfwhm(a::Dict{String, Interface}, l::Light)
    srand(1)
    r = Ray(l, rand())
    n = 1000
    x = Dict("distal_retina" => Dict("s" => 0.0, "weight" => 0.0),"proximal_retina" => Dict("s" => 0.0, "weight" => 0.0))
    for i = 1:n
        Ray(r, l, rand())
        trace!(r, a, x)
    end
    return Dict(k => sqrt(-2*log(abs(v["s"])/v["weight"])) for (k,v) in x)
    end

    getr(x::Circle) = x.cc.r
    getr(x::Ellipse) = getr(x.c)
    getmuy(x::Circle) = x.h - x.r
    getmuy(x::Ellipse) = x.h - x.r2

    function trace!(r::Ray, a::Dict{String, Interface}, x::Dict{String,Dict{String,Float64}})
        k = "nothing"
        for i = 1:8
            k = Z.hit!(r, a, k)
            r"box"(k) && break
            if r"retina"(k)
                yy = r.o[2] - getmuy(a[k].l) - a["mirror"].l.h
                m = sqrt(r.o[1]^2 + yy^2)
                yy /= m
                xradius = getr(a[k].l)
                w = xradius/abs(r.o[1])
                x[k]["s"] += 2w*yy
                x[k]["weight"] += 2w
            end
        end
    end

    function trace!(r::Ray, a::Dict{String, Interface}, x::Vector{Float64}, y::Vector{Float64})
        k = "nothing"
        for i = 1:8
            k = Z.hit!(r, a, k)
            push!(x, r.o[1])
            push!(y, r.o[2])
            r"box"(k) && break
        end
    end

end


c = 13.6
p = MvNormal([0,0], c)
xs = linspace(-5, 5, 100)
ys = linspace(-5, 5, 100)
zs = [pdf(p, [x, y]) for x in xs, y in ys]
    pdf(p, [0,0])/pdf(p, [sqrt(2*log(2))*c, 0])

    surface(xs, ys, zs)
    plot(xs, [pdf(p, [x, 0]) for x in xs])

        n = 10000
        xy = rand(p, n)
        std(xy)

        r = vec(sqrt(sum(xy.^2, 1)))
        std(r)
        =#
