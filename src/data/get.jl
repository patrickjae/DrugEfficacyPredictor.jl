#aliases
get_cell_lines(experiment_id::String) = to_json(collect(values(get_experiment(experiment_id).cell_lines)))
get_drugs(experiment_id::String) = to_json(collect(values(get_experiment(experiment_id).drugs)))
get_genes(experiment_id::String) = to_json(union(collect(values(get_experiment(experiment_id).genes)), collect(values(get_experiment(experiment_id).genes_by_hgnc)), collect(values(get_experiment(experiment_id).genes_by_ensembl))))
get_proteins(experiment_id::String) = to_json(get_experiment(experiment_id).proteins)
get_pathways(experiment_id::String) = to_json(collect(values(get_experiment(experiment_id).pathway_information)))
get_outcomes(experiment_id::String) = to_json(collect(values(get_experiment(experiment_id).results)))

function to_json(data)
    response = Dict{String, Any}()
    response["message"] = "JSON serialization not yet implemented for $(typeof(data))"
    response
end

function to_json(experiment_id::String)
    data = get_experiment(experiment_id)
    response = Dict{String, Any}()
    response["results"] = to_json(collect(values(data.results)))
    response["cell_lines"] = to_json(collect(values(data.cell_lines)))
    response["drugs"] = to_json(collect(values(data.drugs)))
    response["genes"] = to_json(union(collect(values(data.genes)), collect(values(data.genes_by_hgnc)), collect(values(data.genes_by_ensembl))))
    response["proteins"] = to_json(data.proteins)
    response["views"] = map(v -> string(v), data.views)
    response["pathways"] = to_json(collect(values(data.pathway_information)))
    response
end

function to_json(data::Vector)
    ret = Vector{Dict{String, Any}}()
    for obj in data
        push!(ret, to_json(obj))
    end
    ret
end

function to_json(data::CellLine)
    response = Dict{String, Any}()
    response["id"] = data.id
    response["cancer_type"] = data.cancer_type
    response["views"] = Dict{String, Any}()
    for (view_type, view_data) in data.views
        response[string(view_type)] = to_json(view_data)
    end
    response["in_test_set"] = data.in_test_set
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
