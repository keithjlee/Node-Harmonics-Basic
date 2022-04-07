using FastSphericalHarmonics, FastTransforms, Clustering, BoundingSphere

###### We are using the physicist's standard of spherical coordinates
# θ = angle from positive Z pole ∈ [0, π]
# ϕ = angle from positive X in XY plane ∈ [0, 2π]
#Coordinates to a point
function spherical2cartesian(r,θ,ϕ)
    x = @.r*sin(θ)*cos(ϕ)
    y = @.r*sin(θ)*sin(ϕ)
    z = @.r*cos(θ)
    return x, y, z
end

function spherical2cartesian(θ,ϕ)
    x = @.sin(θ)*cos(ϕ)
    y = @.sin(θ)*sin(ϕ)
    z = @.cos(θ)
    return x, y, z
end

#Coordinates to a vector (from origin)
function spherical2vector(θ,ϕ; scale = 1.0)
    # θ is angle from positive Z 
    # ϕ is the angle from positive X
    x = @.scale * sin(θ) * cos(ϕ)
    y = @.scale * sin(θ) * sin(ϕ)
    z = @.scale * cos(θ)

    return Vec3.(x,y,z)
end

#Coordinates to a vector (from origin)
function spherical2point(θ,ϕ; scale = 1.0)
    # θ is angle from positive Z 
    # ϕ is the angle from positive X
    x = @.scale * sin(θ) * cos(ϕ)
    y = @.scale * sin(θ) * sin(ϕ)
    z = @.scale * cos(θ)

    return Point3.(x,y,z)
end

function nodeforce2spherical(nf::nodeForce)
    vectors = nf.startPoints
    value = nf.TC .* norm.(nf.forces)
    phi = [atan(vector[2], vector[1]) for vector in vectors]
    phi[phi .< 0] .+= 2pi
    theta = [acos(vector[3] / norm(vector)) for vector in vectors]

    return [theta phi value]
end

function nodeforce2cartesian(nf::nodeForce; radius = 2)
    x = [sp[1] for sp in nf.startPoints] ./ radius
    y = [sp[2] for sp in nf.startPoints] ./ radius
    z = [sp[3] for sp in nf.startPoints] ./ radius
    value = nf.TC .* norm.(nf.forces)

    return [x y z value]
end

function rotationMatrix(θx, θy, θz)    
    return [1 0 0; 0 cos(θx) -sin(θx); 0 sin(θx) cos(θx)] * [cos(θy) 0 sin(θy); 0 1 0; -sin(θy) 0 cos(θy)] * [cos(θz) -sin(θz) 0; sin(θz) cos(θz) 0; 0 0 1]
end

function rotationMatrixd(θx, θy, θz)
    return [1 0 0; 0 cosd(θx) -sind(θx); 0 sind(θx) cosd(θx)] * [cosd(θy) 0 sind(θy); 0 1 0; -sind(θy) 0 cosd(θy)] * [cosd(θz) -sind(θz) 0; sind(θz) cosd(θz) 0; 0 0 1]
end

# Spherical harmonics analysis
# DISCRETIZATION TO BE CALIBRATED

sampleRate = 181 #number of discretizations of π

θ, ϕ = sph_points(sampleRate) # rows = θ, cols = ϕ

thet = collect(θ)
thet[1] = 0.0
thet[end] = pi

ph = collect(ϕ)
ph[end] = 2pi

frequencyRange = 30 # number of frequencies to sample

const Y = [[[sphevaluate(theta, phi, l, m) for theta in θ, phi in ϕ] for m = -l:l] for l = 0:frequencyRange]

const xsphere = [cos(phi)sin(theta) for theta in thet, phi in ph]
const ysphere = [sin(phi)sin(theta) for theta in thet, phi in ph]
const zsphere = [cos(theta) for theta in thet, phi in ph]

realx = Vector{Vector{Matrix{Float64}}}()
realy = Vector{Vector{Matrix{Float64}}}()
realz = Vector{Vector{Matrix{Float64}}}()

for i = 1:8
    x2 = Vector{Matrix{Float64}}()
    y2 = Vector{Matrix{Float64}}()
    z2 = Vector{Matrix{Float64}}()
    for j = 1:2i-1
        Y2 = abs.(Y[i][j])
        push!(x2, xsphere .* Y2)
        push!(y2, ysphere .* Y2)
        push!(z2, zsphere .* Y2)
    end
    push!(realx, x2)
    push!(realy, y2)
    push!(realz, z2)
end


# function harmonicAnalysis(theta, phi, value; l = 16)
#     idx1 = Int.(round.(theta ./ (pi/sampleRate)))
#     idx1[idx1 .== 0] .= 1
#     idx2 = Int.(round.(phi ./ (2pi/(2sampleRate - 1))))
#     idx2[idx2 .== 0] .= 1

#     sphericalFunction = sparse(idx1, idx2, value, sampleRate, 2sampleRate - 1)
#     coefficients = sph_transform(sphericalFunction)

#     C = Vector{Vector{Float64}}()
#     for i = 0:l
#         push!(C, [coefficients[sph_mode(i,m)] for m = -i:i])
#     end

#     featureVector = [norm(sum(C[i] .* Y[i])) for i = 1:l+1]

#     return featureVector
# end


# closest distance in phi range
function phiDifference(phi1, phi2)
    return min(abs(phi1-phi2), abs(2pi-phi1+phi2))
end



function sphericalGaussian(func::Matrix{Float64}, θ, ϕ; δ = 1)
    val = 0.0

    px = cos(ϕ)sin(θ)
    py = sin(ϕ)sin(θ)
    pz = cos(θ)

    for i = 1:size(func)[1] #for each force
        val += func[i,4] * exp(-δ*((px-func[i,1])^2 + (py-func[i,2])^2 + (pz-func[i,3])^2))
    end
    return val
end


function harmonicAnalysis(nodeforce::nodeForce; l = 16, δ = 5, sampleRate = 181)
    
    #extract matrix of [θ ϕ val]; nRows = number of forces acting at node
    # sph = nodeforce2spherical(nodeforce)
    sph = nodeforce2cartesian(nodeforce)

    #discretization of spherical coordinates
    θ, ϕ = sph_points(sampleRate)

    #convert to a gaussian distribution
    # sphericalFunction = [sphericalGaussian(sph, thet, phi; δ = δ) for thet in θ, phi in ϕ]
    sphericalFunction = [sphericalGaussian(sph, thet, phi; δ = δ) for thet in θ, phi in ϕ]
    #extract expansion coefficients for frequency/harmonic
    coefficients = sph_transform(sphericalFunction)

    #extract scalar coefficient values for each frequency/harmonic
    C = Vector{Vector{Float64}}()
    for i = 0:l
        push!(C, [coefficients[sph_mode(i,m)] for m = -i:i])
    end

    #calculate rotation-invariant energy function
    featureVector = [norm(sum(C[i] .* Y[i])) for i = 1:l+1]

    return featureVector
end

function sphericalFunc(nodeforce::nodeForce, δ; sampleRate = 181)
    # sph = nodeforce2spherical(nodeforce)
    sph = nodeforce2cartesian(nodeforce)

    θ, ϕ = sph_points(sampleRate)

    func = [sphericalGaussian(sph, thet, phi; δ = δ) for thet in θ, phi in ϕ]

    func[abs.(func) .< 1e-10] .= 0.0

    return func
end

function sphericalFunc(nodeforce::Matrix{Float64}, δ; sampleRate = 181)
    θ, ϕ = sph_points(sampleRate)

    func = [sphericalGaussian(nodeforce, thet, phi; δ = δ) for thet in θ, phi in ϕ]

    func[abs.(func) .< 1e-10] .= 0.0

    return func
end

function harmonicAnalysis(sphFunction::Matrix{Float64}; l = 16)

    coefficients = sph_transform(sphFunction)

    C = Vector{Vector{Float64}}()
    for i = 0:l
        push!(C, [coefficients[sph_mode(i,m)] for m = -i:i])
    end

    featureVector = [norm(sum(C[i] .* Y[i])) for i = 1:l+1]

    return featureVector
end

# hold on to just coefficients for visualization purposes
function harmonicCoefficients(nodeforce::nodeForce; l = 16)
    sph = nodeforce2spherical(nodeforce)

    theta = sph[:,1]
    phi = sph[:,2]
    value = sph[:,3]

    idx1 = Int.(round.(theta ./ (pi/sampleRate)))
    idx1[idx1 .== 0] .= 1
    idx2 = Int.(round.(phi ./ (2pi/(2sampleRate - 1))))
    idx2[idx2 .== 0] .= 1

    sphericalFunction = sparse(idx1, idx2, value, sampleRate, 2sampleRate - 1)
    coefficients = sph_transform(sphericalFunction)

    C = Vector{Vector{Float64}}()
    for i = 0:l
        push!(C, [coefficients[sph_mode(i,m)] for m = -i:i])
    end
    return C
end

function harmonicCoefficients(sphFunction::Matrix{Float64}; l = 16)
    coefficients = sph_transform(sphFunction)

    C = Vector{Vector{Float64}}()
    for i = 0:l
        push!(C, [coefficients[sph_mode(i,m)] for m = -i:i])
    end
    return C
end

function comparator(func1::Vector{Matrix{Float64}}, func2::Vector{Matrix{Float64}}, index::Int; cm = :viridis)
    f = Figure(resolution = (1000,500),)
    heatmap(f[1,1], func1[index], axis = (title = "Function", xlabel = L"\theta", ylabel = L"\phi"), colormap = cm)
    heatmap(f[1,2], func2[index], axis = (title = "Harmonic Estimate",xlabel = L"\theta", ylabel = L"\phi"), colormap = cm)
    heatmap(f[1,3], func1[index] .- func2[index], axis = (title = "Error",xlabel = L"\theta", ylabel = L"\phi"), colormap = cm)
    display(f)
end

function forceDemand(nf::nodeForce)
    return sum(abs.(nf.axialForces))
end

function estimateError(nf::nodeForce, l::Int; δ = 20)
    sph1 = sphericalFunc(nf, δ)
    c = harmonicCoefficients(sph1; l = l)
    sph2 = sum([sum(c[i] .* Y[i]) for i = 1:l+1])

    return norm(sph1 .- sph2) / norm(sph1)
end

function estimateError(sph1::Matrix{Float64}, l::Int; δ = 20)
    c = harmonicCoefficients(sph1; l = l)
    sph2 = sum([sum(c[i] .* Y[i]) for i = 1:l+1])

    return norm(sph1 .- sph2) / norm(sph1)
end

function estimatedFunction(sph::Matrix{Float64}, l::Int; δ = 20)
    c = harmonicCoefficients(sph; l = l)
    sph2 = sum([sum(c[i] .* Y[i]) for i = 1:l+1])

    return sph2
end

function estimatedFunction2(sph::Matrix{Float64}, l::Int; δ = 20)
    coefficients = sph_transform(sph)

    C = [coefficients[sph_mode(l-1,m)] for m = -(l-1):(l-1)]

    return C .* Y[l]
end

