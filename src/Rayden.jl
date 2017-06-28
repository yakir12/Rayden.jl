__precompile__()
module Rayden

using CoordinateTransformations, StaticArrays

export Vec, Ray, Ellipsoid, Interface, OpticUnit, raytrace!, advance!, bend!

const Vec = SVector{3,Float64}

mutable struct Ray
    orig::Vec # origin of the ray
    dir::Vec # (normalized) direction of the ray
    Ray(o::Vec, d::Vec) = new(o, normalize(d))
end

Ray() = Ray(Vec(0,0,0), Vec(1,0,0))

struct Ellipsoid
    c::Vec # the center of the ellipsoid
    r::Vec # the radii of the ellipsoid
    dir::Vec # the center (normalized) direction of the dome 
    open::Float64 # cos(α), where α is half the opening angle of the dome
    # all of the following are transformations
    center::Translation{Vec} # translate the ellipsoid to zero
    scale::LinearMap{SArray{Tuple{3,3},Float64,2,9}} # scale the ellipsoid to a unit sphere
    center_scale::AffineMap{SArray{Tuple{3,3},Float64,2,9},SVector{3,Float64}}
    uncenter::Translation{Vec}
    unscale::LinearMap{SArray{Tuple{3,3},Float64,2,9}}
    uncenter_scale::AffineMap{SArray{Tuple{3,3},Float64,2,9},SVector{3,Float64}}
    function Ellipsoid(c::Vec, r::Vec, dir::Vec, open::Float64)
        uncenter = Translation(c)
        unscale = LinearMap(@SMatrix([r[1] 0 0; 0 r[2] 0; 0 0 r[3]])) # fix with Diagonal
        center = inv(uncenter)
        scale = inv(unscale)
        center_scale = scale∘center
        uncenter_scale = inv(center_scale)
        new(c, r, normalize(dir), open, center, scale, center_scale, uncenter, unscale, uncenter_scale)
    end
end

function distance(orig::Vec, dir::Vec)::Tuple{Float64, Float64}
    # distance between ray and the two (possibly identical) intersection points with a centered unit sphere
    b = -orig⋅dir
    disc = b^2 - orig⋅orig + 1
    if disc >= 0
        d = sqrt(disc)
        t2 = b + d
        if t2 >= 0
            t1 = b - d
            return t1 > 0 ? (t1, t2) : (Inf, t2)
        end
    end
    return (Inf, Inf)
end

function advance!(r::Ray, s::Ellipsoid)::Bool
    # move the ray's origin to the intersection point that is within the ellipsoid's dome, and return failiure
    orig = s.center_scale(r.orig)
    dir = normalize(s.scale(r.dir))
    ls = distance(orig, dir)
    for l in ls
        if !isinf(l)
            o = orig + l*dir
            p = s.unscale(o)
            cosα = s.dir⋅normalize(p)
            if cosα > s.open
                r.orig = s.uncenter(p)
                return false
            end
        end
    end
    return false
end

struct Interface
    normal::AffineMap{SArray{Tuple{3,3},Float64,2,9},SVector{3,Float64}}
    n::Float64
    n2::Float64
    Interface(normal::AffineMap{SArray{Tuple{3,3},Float64,2,9},SVector{3,Float64}}, n::Float64) = new(normal, n, n^2)
end

struct OpticUnit
    body::Ellipsoid
    interface::Interface
    register::Bool
    name::String
end


function OpticUnit(body::Ellipsoid, pointin::Bool, n::Float64, register::Bool, name::String)
    i = Float64((-1)^pointin)
    dir = LinearMap(@SMatrix([i 0 0; 0 i 0; 0 0 i]))
    normal = dir∘body.scale∘body.scale∘body.center
    interface = Interface(normal, n)
    return OpticUnit(body, interface, register, name)
end


function bend!(r::Ray, i::Interface)::Bool
    i.n == 1 && return false
    N = normalize(i.normal(r.orig))
    a = -r.dir⋅N
    b = i.n2*(1 - a^2)
    dir = b <= 1 ? i.n*r.dir + (i.n*a - sqrt(1 - b))*N : r.dir + 2a*N
    r.dir = normalize(dir)
    return false
end

raytrace!(r::Ray, c::OpticUnit) = advance!(r, c.body) || bend!(r, c.interface)

end #module

