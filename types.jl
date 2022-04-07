using Base: Float64, Symbol
using Statistics, LinearAlgebra, SparseArrays, GLMakie, GeometryBasics

mutable struct Node
    position :: Vector{Float64} # [x, y, z], z is optional
    DOFS :: Vector{Bool} # [dof1, dof2, dof3, ...] type of analysis deduced from length
    x :: Float64 # individual cartesian coordinates of Node
    y :: Float64 #
    z :: Float64 #
    load :: Vector{Float64} # external load, if applicable
    reaction :: Vector{Float64} # external reaction, if applicable
    disp :: Vector{Float64} # displacements of dofs (after analysis)
    elements :: Vector{Tuple{Int64,Int64}} # [(element index, 0 = start / 1 = end), ...]
    globalIndex :: Vector{Int64} #positions in global DOF order

    # create node with assigned DOFs
    function Node(position, DOFS::Vector{Bool})
        node = new(position, DOFS)
        if length(position) == 3
            node.x, node.y, node.z = node.position
        else
            node.x, node.y = node.position
        end
        node.elements = Vector{Tuple{Int64,Int64}}([])
        return node
    end

    # node with predefined DOFs for typical scenarios
    function Node(position, type::Symbol, fixity::Symbol)
        dims = length(position)
        node = new(position, dofMaker(dims, type, fixity))
        if length(position) == 3
            node.x, node.y, node.z = node.position
        else
            node.x, node.y = node.position
        end
        node.elements = Vector{Tuple{Int64,Int64}}([])
        return node
    end

    # basic node generation (position only)
    function Node(position)
        node = new(position)
        if length(position) == 3
            node.x, node.y, node.z = node.position
        else
            node.x, node.y = node.position
        end
        node.elements = Vector{Tuple{Int64,Int64}}([])
        return node
    end
end


mutable struct Element
    nodeIndex :: Vector{Int64} # index of start/end nodes
    type :: Symbol # :truss or :frame
    posStart :: Vector{Float64} # coordinates of start node
    posEnd :: Vector{Float64} # coordinates of end node
    length :: Float64 # element length
    k :: Symmetric{Float64, Matrix{Float64}} # stiffness matrix
    dofIndex ::Vector{Int64} # DOF indices in global expanded coordinates
    localForce :: Vector{Float64}  #list of internal loads in local coordinates
    axialForce :: Float64 #-ve compression, +ve tension
    globalForce :: Vector{Float64}
    A :: Float64 # section area
    E :: Float64 # modulus of elasticity
    G :: Float64 # shear modulus, OPTIONAL
    Iz :: Float64 # strong axis moment of inertia, OPTIONAL
    Iy :: Float64 # weak axis moment of inertia, OPTIONAL
    J :: Float64 # torsional constant, OPTIONAL
    Ψ :: Float64 # angle of roll w/r/t local x axis, for 3d frames
    R :: Matrix{Float64} # rotation matrix
    LCS :: Vector{Vec3{Float64}} #local coordinate system

    # Basic generator without material properties
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, type::Symbol)
        element = new(nodeIndex)
        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.type = type
        return element
    end

    # truss generator, 2 and 3 dimensions
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, E::Float64, A::Float64)
        element = new(nodeIndex)
        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.E = E
        element.A = A
        element.type = :truss
        return element
    end

    # 2D frame generator
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, E::Float64, A::Float64, I::Float64)
        element = new(nodeIndex)
        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.E = E
        element.A = A
        element.Iz = I
        element.type = :frame
        return element
    end

    # 3D frame generator
    function Element(nodes::Vector{Node}, nodeIndex::Vector{Int64}, E::Float64, A::Float64, G::Float64, Iz::Float64, Iy::Float64, J::Float64)
        element = new(nodeIndex)
        element.posStart = nodes[nodeIndex[1]].position
        element.posEnd = nodes[nodeIndex[2]].position
        element.length = norm(element.posEnd .- element.posStart)
        element.E = E
        element.A = A
        element.G = G
        element.Iz = Iz
        element.Iy = Iy
        element.J = J
        element.Ψ = pi/2 #default
        element.type = :frame
        element.LCS = lcs(element, element.Ψ)
        return element
    end

end

# creates a Load type with proper DOF indexing
mutable struct Load
    position :: Vector{Float64} #approximate position of load
    load :: Vector{Float64} #value of load [p1, p2, p3, ..., pdof]
    index :: Int64 #index of node assigned to load

    function Load(nodes::Vector{Node}, position::Vector{Float64}, load::Vector{Float64}; tol = 0.5)
        load = new(position, load)
        load.index = findmin([norm(n.position .- load.position) for n in nodes])[2]
        return load
    end
end

mutable struct Structure
    nodes :: Vector{Node} # array of Nodes
    elements :: Vector{Element} # array of Elements
    loads :: Vector{Load} # array of Loads
    dims :: Int64 # 2, 3
    displacedNodes :: Vector{Node} # array of displaced nodes
    displacedElements :: Vector{Element} # array of displaced elements
    DOFS :: Vector{Bool} # degrees of freedom 
    K :: SparseMatrixCSC{Float64,Int64} # global stiffness matrix
    F :: Vector{Float64} # external loads
    U :: Vector{Float64} # Displacements after analysis
    reactions :: Vector{Float64} # external reactions
    compliance :: Float64 # compliance after analysis
    nNodes :: Int64 # number of nodes
    nElements :: Int64 #number of elements
    nDOFS :: Int64 #number of DOFs (total)
    freeDOFS :: Vector{Int64} #array of free DOFs

    #general structure definition with all necesary components
    function Structure(nodes::Vector{Node}, elements::Vector{Element}, loads::Vector{Load})

        # initialize
        structure = new(nodes, elements, loads)

        # deduce dimensions
        structure.dims = length(nodes[1].position) 

        # global DOF vector
        structure.DOFS = vcat([node.DOFS for node in nodes]...)

        # general information
        structure.nNodes = length(nodes)
        structure.nElements = length(elements)
        structure.nDOFS = length(structure.DOFS)
        structure.freeDOFS = findall(structure.DOFS)

        # create global load vector
        structure.F = zeros(structure.nDOFS)
        nodalDOFLength = length(nodes[1].DOFS)

        for load in loads
            idx = load.index * nodalDOFLength - (nodalDOFLength - 1) .+ collect(0:nodalDOFLength-1)
            structure.F[idx] .= load.load
        end

        return structure
    end

    # Basic structure definition without loads
    function Structure(nodes::Vector{Node}, elements::Vector{Element})
        structure = new(nodes, elements)
        structure.dims = length(nodes[1].position)
        structure.DOFS = vcat([node.DOFS for node in nodes]...)

        structure.nNodes = length(nodes)
        structure.nElements = length(elements)
        structure.nDOFS = length(structure.DOFS)
        structure.freeDOFS = findall(structure.DOFS)

        return structure
    end

end

const forceDic = Dict([-1 => 1:3, 1 => 7:9])
const forceDicTruss = Dict([-1 => 1:3, 1 => 4:6])

mutable struct nodeForce
    nodeID :: Int64
    elementID :: Vector{Int64}
    position :: Point3{Float64}
    startPoints :: Vector{Point3{Float64}}
    centerLines :: Vector{Point3{Float64}}
    forces :: Vector{Vec3{Float64}}
    theta::Vector{Float64}
    phi::Vector{Float64}
    axialForces :: Vector{Float64}
    TC :: Vector{Float64}
    normalizingFactor :: Float64

    function nodeForce(ID::Int64, structure::Structure; factor = 2)
        # initialize
        nodeforce = new(ID)

        # extract node
        node = structure.nodes[ID]

        nodeforce.position = node.position
        # indices of connected elements
        nodeforce.elementID = [id[1] for id in node.elements]

        # extract elements
        elements = structure.elements[nodeforce.elementID]

        # does element start or end at this node?
        startEnd = [id[2] for id in node.elements]

        # tension or compression?
        TC = sign.([element.axialForce for element in elements])

        # axial force vectors acting on node
        if elements[1].type == :frame
            forces = - [Vec3(elements[i].globalForce[forceDic[startEnd[i]]]) for i = 1:length(elements)]
        else
            forces = - [Vec3(elements[i].globalForce[forceDicTruss[startEnd[i]]]) for i = 1:length(elements)]
        end
        # add external load if applicable
        if isdefined(node, :load)
            push!(forces, Vec3(node.load[1:3]))
            push!(TC, 1.0)
        end

        # add external reactions if applicable
        if isdefined(node, :reaction)
            push!(forces, Vec3(node.reaction[1:3]))
            push!(TC, -1.0)
        end

        # update fields
        nodeforce.forces = forces
        nodeforce.TC = TC
        nodeforce.axialForces = norm.(forces) .* TC

        # normalization for visualization
        nodeforce.normalizingFactor = maximum(norm.(nodeforce.forces))
        
        # remove zero force members (messes with analysis)
        zeroIndices = findall(norm.(forces) .< 1e-3)
        deleteat!(nodeforce.forces, zeroIndices)
        deleteat!(nodeforce.TC, zeroIndices)
        deleteat!(nodeforce.axialForces, zeroIndices)

        # start points as function of tension/compression and load orientation
        nodeforce.startPoints = [Point3(TC[i] .* factor .* nodeforce.forces[i].data ./ norm(nodeforce.forces[i].data)) for i = 1:length(nodeforce.forces)]

        # dashed lines for clearer visualization
        nodeforce.centerLines = vcat([[nodeforce.startPoints[i], Point3([0,0,0])] for i = 1:length(nodeforce.startPoints)]...)

        return nodeforce
    end
end

mutable struct Geometry
    structure :: Structure
    nodes :: Vector #{Point3{Float64}}
    nodalForces :: Vector{nodeForce}
    elements :: Vector{Vector} #Vector{Vector{Point3{Float64}}}
    axialForce :: Vector{Float64}
    areas :: Vector{Float64}
    displacedNodes :: Vector #Vector{Point3{Float64}}
    displacedElements :: Vector{Vector} #Vector{Vector{Point3{Float64}}}
    midpoints :: Vector #Vector{Point3{Float64}}
    loads :: Vector{Vector{Float64}}
    reactions :: Vector #Vector{Vec3{Float64}}
    nodeLabels :: Vector{String}
    elementLabels :: Vector{String}
    pinDOFS :: NTuple #{6, Vector{Float64}}
    fixDOFS :: NTuple #{6, Vector{Float64}}

    function Geometry(structure::Structure; nodeAnalysis = false, SF = :auto)
        geometry = new(structure)

        if structure.dims == 3
            geometry.nodes = Point3.([node.position for node in structure.nodes])
            geometry.elements = [[Point3(element.posStart), Point3(element.posEnd)] for element in structure.elements]
        else
            geometry.nodes = Point3.([[node.position; 0.0] for node in structure.nodes])
            geometry.elements = [[Point3([element.posStart; 0.0]), Point3([element.posEnd; 0.0])] for element in structure.elements]
        end

        areas = [element.A for element in structure.elements]
        geometry.areas = areas ./ maximum(areas)
        # geometry.midpoints = [(e.posStart + e.posEnd) / 2 for e in structure.elements]
        
        geometry.midpoints = sum.(geometry.elements) ./ 2

        geometry.nodeLabels = ["N" * string(i) for i = 1:structure.nNodes]
        geometry.elementLabels = ["E" * string(i) for i = 1:structure.nElements]

        geometry.pinDOFS = pinDofConverter(structure)
        geometry.fixDOFS = fixDofConverter(structure)

        if isdefined(structure, :displacedNodes)
            if structure.dims == 3
                geometry.displacedNodes = Point3.([node.position for node in structure.displacedNodes])
                geometry.displacedElements = [[Point3(element.posStart), Point3(element.posEnd)] for element in structure.displacedElements]
            else
                geometry.displacedNodes = Point3.([[node.position; 0.0] for node in structure.displacedNodes])
                geometry.displacedElements = [[Point3([element.posStart; 0.0]), Point3([element.posEnd; 0.0])] for element in structure.displacedElements]
            end
        end

        if isdefined(structure, :loads)
            geometry.loads = loadConverter(structure; scaleFactor = SF)
        end

        if isdefined(structure.elements[1], :localForce)
            geometry.axialForce = [element.axialForce for element in structure.elements]
            geometry.reactions = reactionConverter(structure)
        end

        if nodeAnalysis
            geometry.nodalForces = [nodeForce(i, structure) for i = 1:structure.nNodes]
        end

        return geometry
    end
        
end