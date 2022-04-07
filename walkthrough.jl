## parsing JSON files for example structures
begin
    include("types.jl")
    include("structuralAnalysis.jl")
    include("geometricAnalysis.jl")
    include("visualization.jl")
    include("visualization2.jl")
    include("parsing.jl")
end

# Look in parsing.jl for how a external structure is parsed, defined, and analyzed
structure = renaud2structure("00053");

# converting to Geometry types for visualization
geo = Geometry(structure; nodeAnalysis = true)

# interactive visualization
fig = structurePlot(geo; colorMap = red2blue, scaleToForce = true)

```
This is the brunt of our analysis.
Converts a sets of points to a Gaussian spherical function using sphericalGaussian()
Then expands that function into its spherical harmonic weights using sph_transform()
Then norms within each frequency to extract feature vectors.
```
featureVecs = [harmonicAnalysis(nf; l = 16, δ = 20) for nf in geo.nodalForces];

# cost matrix
costs = [norm(fv1 - fv2) for fv1 in featureVecs, fv2 in featureVecs];

# visualize distance matrix
costmap = heatmap(reverse(costs, dims = 1), 
    axis = (aspect = 1,
        xlabel = "Node ID",
        ylabel = "Node ID"))


# Clustering 
nclusters = 10

#concatenate all feature vectors
fvMatrix = hcat(featureVecs...)

# results of k-means clustering
kresults =  kmeans(fvMatrix, nclusters)

# assign each node to a cluster
kclusters = [findall(kresults.assignments .== i) for i = 1:nclusters]

#visualize
kclusterplot = clusterPlot(geo, kclusters)


# feature vectors in each cluster
FV_k_clusters = [featureVecs[set] for set in kclusters]

# initial complexity score
radii_all = boundingsphere(featureVecs)[2]

# Complexity score for each cluster
kClusterRadii = [boundingsphere(vecs)[2] for vecs in FV_k_clusters]
