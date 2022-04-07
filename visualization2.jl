function axo(geo::Geometry; 
    cmforce = redWhiteBlue,
    mode = :undisplaced,
    forces = true,
    scale = false,
    elev = pi/6,
    az = pi/4,
    lw = 2,
    lcoverride = nothing,
    bgc = :transparent,
    arrowSize = 0.15,
    res = (3000,3000))

    fig = Figure(resolution = res, backgroundcolor = bgc)

    if scale == false
        ax1 = Axis3(fig[1,1],
            backgroundcolor = bgc,
            elevation = elev,
            azimuth = az,
            aspect = :data)
        hidedecorations!(ax1)
        hidespines!(ax1)
    else
        ax1 = Axis3(fig[1,1],
            backgroundcolor = bgc,
            elevation = elev,
            azimuth = az,
            aspect = :data)
    end

    if mode == :undisplaced
        els = geo.elements
        normalizedAreas = [e.A for e in geo.structure.elements]
        normalizedAreas ./= maximum(normalizedAreas)
        thickness = lw .* normalizedAreas
        lc = - ones(geo.structure.nElements)
        cm = :grays
        cr = (-1,1)
    else
        els = geo.displacedElements
        peakforce = maximum(abs.(geo.axialForce))
        thickness = 5 .* lw .* abs.(geo.axialForce) ./ peakforce
        lc = geo.axialForce
        cm = cmforce
        cr = (-peakforce, peakforce)
    end

    # override linecolor if desired
    if !isnothing(lcoverride)
        lc = lcoverride
    end

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    if mode == :undisplaced && forces
        arrows!(geo.loads...,
            color = (niceRed, 0.5),
            arrowsize = arrowSize,
            linewidth = arrowSize / 3)
    end

    ax = Axis3(fig[2,1],
        backgroundcolor = bgc,
        elevation = 0,
        azimuth = 0,
        aspect = :data)

    hidespines!(ax)
    hidedecorations!(ax)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    ax3 = Axis3(fig[2,2],
        backgroundcolor = bgc,
        elevation = 0,
        azimuth = pi/2,
        aspect = :data)

    hidedecorations!(ax3)
    hidespines!(ax3)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    ax4 = Axis3(fig[1,2],
        backgroundcolor = bgc,
        elevation = pi/2,
        azimuth = pi/2,
        aspect = :data)
        
    hidedecorations!(ax4)
    hidespines!(ax4)

    linesegments!(vcat(els...),
        color = lc,
        colormap = cm,
        colorrange = cr,
        linewidth = thickness)

    # arrows!(geo.pinDOFS...,
    #     arrowhead = pin,
    #     arrowcolor = (:black,0.6),
    #     arrowsize = 2 * arrowSize,
    #     linewidth = 0)

    return fig
end

function internalForcePlot(geo::Geometry;
    res = (3000,3000),
    bgc = :transparent,
    cm = red2blue,
    elev = pi/6,
    az = pi/4,
    arrowSize = 0.1,
    lw = 1,
    lc = :black,
    showElements = true,
    showNodes = true)

    fig = Figure(resolution = res, backgroundcolor = bgc)

    ax = Axis3(fig[1,1],
        backgroundcolor = :transparent,
        elevation = elev,
        azimuth = az,
        aspect = :data)

    hidedecorations!(ax)
    hidespines!(ax)

    smallestLength = minimum([e.length for e in geo.structure.elements])
    reductionFactor = smallestLength / 4

    if showElements
        linesegments!(vcat(geo.elements...),
            color = (lc, 0.5),
            linewidth = lw)
    end

    if showNodes
        meshscatter!(geo.nodes,
            markersize = arrowSize / 2,
            color = lc)
    end

    

    colorlimit = maximum([maximum(abs.(n.axialForces)) for n in geo.nodalForces])

  
    for n in geo.nodalForces
        s = copy(n.startPoints) .* reductionFactor
        s .+= n.position

        arrows!(s, normalize.(n.forces) .* reductionFactor,
            color = n.TC,
            colormap = cm,
            arrowsize = arrowSize,
            linewidth = arrowSize / 2,
            colorrange = (-colorlimit, colorlimit)
            )
    end

    return fig
end

function allForcePlot(nodeforces::Vector{nodeForce}; 
    res = (10000,10000), 
    elev = pi/6,
    az = pi/4,
    dim = :auto, 
    summed = true,
    cm = redWhiteBlue,
    mc = :black)
    forcefig = Figure(resolution = res,backgroundcolor = :transparent)

    if dim == :auto
        l = Int(floor(sqrt(length(nodeforces))))
        dif = length(nodeforces) - l^2
        if dif != 0
            println("$dif nodes omitted.")
        end
        dims = l
    elseif dim^2 > length(nodeforces)
        error("dim exceeds total number of nodes.")
    else
        dims = dim
    end

    nf_reshaped = reshape(nodeforces[1:dims^2], dims, dims)
    if summed
        for i = 1:dims
            for j = 1:dims
                nf = nf_reshaped[i, j]
                forces = nf.forces ./ nf.normalizingFactor

                ax = Axis3(forcefig[i,j],
                    backgroundcolor = :transparent,
                    elevation = elev,
                    azimuth = az,
                    aspect = :data)
                hidedecorations!(ax)
                hidespines!(ax)

                sumforce = arrows!(vcat([Point3([0.0,0.0,0.0])], cumsum(Point3.([f.data for f in forces]))[1:end-1]), 
                    forces,
                    arrowsize = 0.08,
                    linewidth = 0.03,
                    color = nf.TC,
                    colormap = cm)
            end
        end
    else
        for i = 1:dims
            for j = 1:dims
                nf = nf_reshaped[i, j]
                forces = nf.forces ./ nf.normalizingFactor .* 1.5
                ax = Axis3(forcefig[i,j],
                    elevation = elev,
                    azimuth = az,
                    backgroundcolor = :transparent,
                    aspect = :data)
                hidedecorations!(ax)
                hidespines!(ax)

                sumforce = arrows!(nf.startPoints, 
                    forces,
                    arrowsize = 0.4,
                    linewidth = 0.1,
                    show_axis = false,
                    color = nf.TC,
                    colormap = cm)
                
                meshscatter!(Point3([0.0,0.0,0.0]), color = mc, markersize = 0.1)
            end
        end

    end
    return forcefig
end

function clusterPlot(geo::Geometry,
        clusters::Vector{Vector{Int64}};
        elev = pi/6,
        az = pi/4,
        lc = :black,
        lw = 2,
        res1 = (3000,3000),
        res2 = (6000,3000),
        shade = true,
        cm = :tab20,
        ms = .2,
        FV = nothing)

    if !isnothing(FV)
        res = res2
    else
        res = res1
    end

    fig = Figure(resolution = res,
        backgroundcolor = :transparent)
    
    if !isnothing(FV)
        ax = Axis3(fig[1,1:2],
            elevation = elev,
            azimuth = az,
            backgroundcolor = :transparent,
            aspect = :data)
        hidedecorations!(ax)
        hidespines!(ax)
    else
        ax = Axis3(fig[1,1],
            elevation = elev,
            azimuth = az,
            backgroundcolor = :transparent,
            aspect = :data)
        hidedecorations!(ax)
        hidespines!(ax)
    end

    # plot elements
    linesegments!(vcat(geo.elements...),
        color = (lc, 0.6),
        linewidth = lw)
    
    nclusters = length(clusters)
    colorspace = [cgrad(cm, [0.0,1.0])[z] for z ∈ range(0.0, 1.0, length = nclusters)]

    # plot clusters
    for i = 1:nclusters
        meshscatter!(geo.nodes[clusters[i]],
            color = colorspace[i],
            shading = shade,
            markersize = ms)
    end

    if !isnothing(FV)

        ax2 = Axis(fig[1,3],
            backgroundcolor = :transparent,
            aspect = 1)

        hidedecorations!(ax2)
        hidespines!(ax2)

        allradius = boundingsphere(FV)[2]
        clusterradius = [boundingsphere(FV[clust])[2] for clust in clusters]

        thet = collect(0:0.01:2pi)
        x = cos.(thet)
        y = sin.(thet)

        lines!(x .* allradius, y .* allradius,
            color = :black,
            linewidth = lw * 2,
            linestyle = :dash)

        for i = 1:nclusters
            lines!(x .* clusterradius[i], y .* clusterradius[i],
                color = colorspace[i],
                linewidth = lw)
        end
    end

    return fig
end

