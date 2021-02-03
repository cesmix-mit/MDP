struct Cartesian{T}
    x::T
    y::T
    z::T
end

struct Spherical{T}
    r::T
    θ::T
    ϕ::T
end

Base.:+(f::Cartesian, g::Cartesian) = Cartesian(f.x + g.x, f.y + g.y, f.z + g.z)
Base.:-(f::Cartesian, g::Cartesian) = Cartesian(f.x - g.x, f.y - g.y, f.z - g.z)
Base.:norm(f::Cartesian) = sqrt(f.x^2 + f.y^2 + f.z^2)
Base.:normsq(f::Cartesian) = (f.x^2 + f.y^2 + f.z^2)

function convert(::Type{Spherical}, coord::Cartesian)
    x, y, z = coord.x, coord.y, coord.z
    r = sqrt(x^2 + y^2 + z^2)
    θ = acos(z / r)
    ϕ = atan(y, x)
    return Spherical(r, θ, ϕ)
end

function convert(::Type{Cartesian}, coord::Spherical)
    r, θ, ϕ = coord.r, coord.θ, coord.ϕ
    x = r * sin(θ) * cos(ϕ)
    y = r * sin(θ) * sin(ϕ)
    z = r * cos(θ)
    return Cartesian(x, y, z)
end

function rotc(coord::Cartesian{T}, m::SMatrix{3, 3, T}) where T
    c = @SArray [coord.x, coord.y, coord.z]
    c′ = m*c
    Cartesian(c'...)
end

function euler2rotm(α, β, γ)
    r11 =  cos(α) * cos(β) * cos(γ) - sin(α) * sin(γ)
    r12 = -cos(α) * cos(β) * sin(γ) - sin(α) * cos(γ)
    r13 =  cos(α) * sin(β)
    r21 =  sin(α) * cos(β) * cos(γ) + cos(α) * sin(γ) 
    r22 = -sin(α) * cos(β) * sin(γ) + cos(α) * cos(γ)
    r23 =  sin(α) * sin(β)
    r31 = -sin(β) * cos(γ)
    r32 =  sin(β) * cos(γ)
    r33 =  cos(β)

    @SMatrix [
        r11 r12 r13
        r21 r22 r23 
        r31 r32 r33
    ] 
end
