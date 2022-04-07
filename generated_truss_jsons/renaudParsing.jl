using JSON

function renaud2structure(structureID::String; A = 0.001,
    E = 200e6, #kN/m^2
    density = 80, #kN/m^3
    directory = "exampleStructures/generated_truss_jsons/space_truss_")

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