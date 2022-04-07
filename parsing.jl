using Base: StackFrame
using CSV, JSON

function renaud2structure(structureID::String; A = 0.001,
    E = 200e6, #kN/m^2
    density = 80, #kN/m^3
    directory = "generated_truss_jsons/space_truss_")

    filename = directory * structureID * ".json"
    truss = JSON.parsefile(filename)

    # nodes and elements
    ns = truss["node_list"]
    es = truss["element_list"]

    # Create nodes
    nodes = Vector{Node}()
    for i = 1:length(ns)
        n = ns[i]
        position = [n["point"]["X"], n["point"]["Y"], n["point"]["Z"]]
        if n["is_grounded"] == 1
            fixtype = :fixed
        else
            fixtype = :free
        end
        push!(nodes, Node(position, :truss, fixtype))
    end

    # Create elements
    elements = Vector{Element}()
    for e in es
        indices = Int.(e["end_node_ids"] .+ 1)
        push!(elements, Element(nodes, indices, E, A))
    end

    # Create self-weight loads
    lengths = [e.length for e in elements]
    P = mean(lengths) * A * density / 2
    Pvec = [0, 0, -1.0] # 1kN downwards force 
    loads = Vector{Load}()
    for node in nodes
        if any(node.DOFS .== false)
            continue
        else
            push!(loads, Load(nodes, node.position, Pvec))
        end
    end
    structure = Structure(nodes, elements, loads)
    analyze(structure)
    return structure
end


function compas2structure(name::String; directory = "exampleStrutures/")

    filename = directory * name * ".json"
    test1 = JSON.parsefile(filename)

    # nodes (in mm)
    ns = test1["node_list"]

    # elements
    es = test1["element_list"]

    # material properties
    mp = test1["material_properties"]
    E = Float64(mp["youngs_modulus"]) # kn/cm2
    A = mp["cross_sec_area"] #cm^2
    G = Float64(mp["shear_modulus"]) #kN/cm^2
    Iz = mp["Iz"] #cm^4
    Iy = mp["Iy"] # cm^4
    J = mp["Jx"] #cm^4

    # Create nodes
    NODES = Vector{Node}()

    # for DJMM Bridge, flip structure

    for i = 1:length(ns)
        n = ns[i]
        position = [n["point"]["X"], n["point"]["Y"], n["point"]["Z"]] ./ 10
        if any(fixidx .== i)
            fixtype = :fixed
        else
            fixtype = :free
        end
        push!(NODES, Node(position, :frame, fixtype))
    end

    # Create elements
    ELEMENTS = Vector{Element}()
    for e in es
        indices = e["end_node_ids"] .+ 1
        push!(ELEMENTS, Element(NODES, indices, E, A, G, Iz, Iy, J))
    end

    # Create self-weight loads
    lengths = [e.length for e in ELEMENTS]
    P = mean(lengths) * A * mp["density"]/100^3 / 2
    Pvec = [0, 0, -P, 0, 0, 0]

    LOADS = Vector{Load}()
    for node in NODES
        push!(LOADS, Load(NODES, node.position, Pvec))
    end

    STRUCTURE = Structure(NODES, ELEMENTS, LOADS)

    analyze(STRUCTURE)

    return STRUCTURE
end

allids = Vector{String}()
for i = 0:99 
    if i <10
        push!(allids, "0000" * string(i))
    else
        push!(allids, "000" * string(i))
    end
end