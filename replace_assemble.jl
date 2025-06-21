#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

# Function to process files
function replace_assemble_in_file(filepath)
    content = read(filepath, String)
    
    # Replace the function definition
    content = replace(content, r"assemble\(systems::Vector\{(?:NonlinearSystem|ModelingToolkit\.System)\}; equations=\[\], name\)\s*=\s*(?:ModelingToolkit\.)?compose\(NonlinearSystem\(equations, \[\], \[\], name=name\), systems\)" => 
        "# Removed custom assemble function - using ModelingToolkit.compose directly")
    
    # Replace all usages of assemble with compose
    content = replace(content, r"@named\s+(\w+)\s*=\s*assemble\(\[(.*?)\],\s*equations\s*=\s*\[(.*?)\]\)" => 
        s"@named \1 = ModelingToolkit.compose(NonlinearSystem([\3], [], [], name=:\1), [\2])")
    
    write(filepath, content)
    println("Processed: $filepath")
end

# Process all files that might contain assemble
function process_all_files()
    base_dir = "/home/adrien/Dokument/GitHub/OptimalGshpFlow.jl"
    
    for (root, dirs, files) in walkdir(base_dir)
        for file in files
            if endswith(file, ".jl")
                filepath = joinpath(root, file)
                # Check if file contains "assemble"
                if !isnothing(match(r"assemble\(", read(filepath, String)))
                    replace_assemble_in_file(filepath)
                end
            end
        end
    end
end

process_all_files()
