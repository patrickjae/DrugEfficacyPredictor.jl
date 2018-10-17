function to_json(data)
    response = Dict{String, Any}()
    response["message"] = "JSON serialization not yet implemented for $(typeof(data))"
    response
end

function to_json(data::Experiment)
    response = Dict{String, Any}()
    response["results"] = Vector{Dict{String, Any}}()
    for (_, outcome) in data.results
        push!(response["results"], to_json(outcome))
    end
    response["test_results"] = Vector{String, Any}()
    for (_, test_outcome) in data.rest_results
        push!(response["test_results"], to_json(test_outcome))
    end
    response["cell_lines"] = Vector{Dict{String, Any}}()
    for cl in collect(values(data.cell_lines))
        push!(response["cell_lines"], to_json(cl))
    end
    response["drugs"] = Vector{Dict{String, Any}}()
    for drug in collect(values(data.drugs))
        push!(response["drugs"], to_json(drug))
    end
    response["genes"] = Vector{Dict{String, Any}}()
    all_genes = union(collect(values(data.genes)), collect(values(data.genes_by_hgnc)), collect(values(data.genes_by_ensembl)))
    for gene in all_genes
        push!(response["genes"], to_json(gene))
    end
    response["proteins"] = Vector{Dict{String, Any}}()
    for protein in collect(values(data.proteins))
        push!(response["proteins"], to_json(protein))
    end
    response["views"] = map(v -> string(v), data.views)
    response["pathways"] = Vector{Dict{String, Any}}()
    for pw in data.pathway_information
        push!(response["pathways"], to_json(pw))
    end
    response
end

function to_json(data::CellLine)
    response = Dict{String, Any}()
    response["id"] = data.id
    response["cancer_type"] = data.cancer_type
    response["views"] = Dict{String, Any}()
    for (view_type, view_data) in cl.views
        response[view_type] = to_json(view_data)
    end
    response
end

function to_json(data::Gene)
    response = Dict{String, Any}()
    if data.entrez_id != -1 response["entrez_id"] = data.entrez_id end
    if data.hgnc_id != "" response["hgnc_id"] = data.hgnc_id end
    if data.ensembl_id != "" response["ensembl_id"] = data.ensembl_id end
    response["cancer_gene"] = data.cancer_gene
    response
end

function to_json(data::Protein)
    response = Dict{String, Any}()
    response["hgnc_id"] = data.hgnc_id
    response["antibody_validated"] = data.antibody_validated
end

function to_json(data::Drug)
    response = Dict{String, Any}()
    response["id"] = data.id
    if length(data.affected_genes) > 0
        response["affected_genes"] = Vector{Dict{String, Any}}()
        for gene in data.affected_genes
            push!(response["affected_genes"], to_json(gene))
        end
    end
    if data.chemical_structure != ""
        response["chemical_structure"] = data.chemical_structure
    end
    response
end

function to_json(data::Outcome)
    response = Dict{String, Any}()
    response["drug_id"] = data.drug_id
    response["outcome_values"] = Dict{String, Any}()
    for (cl, val) in data.outcome_values
        response["outcome_values"][cl.id] = val
    end
    response["outcome_type"] = data.outcome_type
    response
end

function to_json(data::Pathway)
    response = Dict{String, Any}()
    response["id"] = data.id
    response["name"] = data.name
    response["genes"] = Vector{Dict{String, Any}}()
    for gene in data.genes
        push!(response["genes"], to_json(gene))
    end
    response
end