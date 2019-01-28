"""
Creates a Cell line object from JSON. Creates Gene (or Protein) objects on the fly as needed.
"""
function create_cell_line(experiment_id::String, data::Dict{String, Any}; for_prediction::Bool=false)
	experiment = get_experiment(experiment_id)
	log_message("fetching/creating cell line...")
	cl = get_cell_line(experiment_id, data, for_prediction = for_prediction)
	if haskey(data, "views")
		# construct the views
		views = data["views"]
		# should be a dict of views
		for (k,v) in views
			view_type = eval(Meta.parse(k))
			if view_type == RPPA
				key_type = Protein
			else
				key_type = Gene
			end
			data_view = get_dataview!(cl, view_type, DataView{key_type, view_type}(cl.id))
			populate_data_view!(experiment_id, data_view, v)
			add_view!(experiment, view_type)
		end

	end
end

function create_cell_lines(experiment_id::String, data::Dict{String, Any}; for_prediction::Bool=false)
    log_message("importing cell lines")
    if haskey(data, "cell_lines")
        cell_lines = data["cell_lines"]
        for cl in cell_lines
            create_cell_line(experiment_id, cl, for_prediction = for_prediction)
        end
    end
end

function create_genes(experiment_id::String, data::Dict{String, Any})
    log_message("importing genes")
    if haskey(data, "genes")
        for gene in data["genes"]
            add_gene(experiment_id, gene)
        end
    end
end

function create_proteins(experiment_id::String, data::Dict{String, Any})
    log_message("importing proteins")
    if haskey(data, "proteins")
        for protein in data["proteins"]
            add_protein(experiment_id, protein)
        end
    end
end

function create_drugs(experiment_id::String, data::Dict{String, Any})
    log_message("importing drugs")
    if haskey(data, "drugs")
        for drug in data["drugs"]
            add_drug(experiment_id, drug)
        end
    end
end

function create_outcomes(experiment_id::String, data::Dict{String, Any})
    log_message("importing outcomes")
    if !haskey(data, "outcome_type")
        throw(ArgumentError("No outcome type specified."))
    end
    outcome_type = data["outcome_type"]
    # TODO: change to "outcomes"
    if !haskey(data, "outcomes")
        throw(ArgumentError("No outcomes specified, provide a dictionary of drug => outcome under the identifier 'outcomes'"))
    end
    outcomes = data["outcomes"]
	experiment = get_experiment(experiment_id)
    for (drug_name, outcome_data) in outcomes
        log_message("importing outcome for drug $drug_name")
        drug = get_drug!(experiment, drug_name)
        add_outcome(experiment_id, drug, outcome_data, outcome_type = outcome_type)
    end
end

function create_pathways(experiment_id::String, data::Dict{String, Any})
    log_message("importing pathway information")
    if haskey(data, "pathways")
        pathways = data["pathways"]
        for pathway_data in pathways
            add_pathway(experiment_id, pathway_data)
            log_message("done for $(pathway_data["name"])")
        end
    end
end

"""
Populates a data view with ExomeSeq data from JSON data.
If no protein change is provided, we assume the empty string.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene,ExomeSeq}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)

		# mandatory values
        if !haskey(entry, "protein_change")
            throw(ArgumentError("You have to specify at least a protein change for ExomeSeq data."))
        end

		# if !(haskey(entry,"reference_mismatch_avg") && haskey(entry,"variant_mismatch_avg")
		# 		&& haskey(entry,"reference_mismatch_sum") && haskey(entry, "variant_mismatch_sum")
		# 		&& haskey(entry, "reference_dist3effective_avg") && haskey(entry, "variant_dist3effective_avg"))
		# 	throw(ArgumentError("For ExomeSeq data, you have to specify at least values for reference_mismatch_avg, variant_mismatch_avg,
		# 							reference_mismatch_sum, variant_mismatch_sum, reference_dist3effective_avg and variant_dist3effective_avg!"))
		# end
        protein_change = entry["protein_change"]

		# optional values
        reference_mismatch_avg = get(entry, "reference_mismatch_avg", 0.)
        reference_mismatch_sum = get(entry, "reference_mismatch_sum", 0.)
        reference_dist3effective_avg = get(entry, "reference_dist3effective_avg", 0.)
        variant_mismatch_avg = get(entry, "variant_mismatch_avg", 0.)
        variant_mismatch_sum = get(entry, "variant_mismatch_sum", 0.)
        variant_dist3effective_avg = get(entry, "variant_dist3effective_avg", 0.)

		num_cosmic = get(entry, "num_cosmic", 0)
		is_cancer_gene = get(entry, "is_cancer_gene", false)
		variant_effect = get(entry, "variant_effect", "")
		nucleotid_change = get(entry, "nucleotid_change", "")
		variant_confidence = get(entry, "variant_confidence", 1.)
		norm_zygosity = get(entry, "norm_zygosity", "")
		norm_reference_count = get(entry, "norm_reference_count", 0.)
		norm_variant_count = get(entry, "norm_variant_count", 0.)
		tumor_zygosity = get(entry, "tumor_zygosity", "")
		tumor_reference_count = get(entry, "tumor_reference_count", 0.)
		tumor_variant_count = get(entry, "tumor_variant_count", 0.)
		details = get(entry, "details", "")

		exome_data = ExomeSeq(protein_change, reference_mismatch_sum=reference_mismatch_sum, reference_mismatch_avg=reference_mismatch_avg,
                                reference_dist3effective_avg=reference_dist3effective_avg, variant_mismatch_sum=variant_mismatch_sum,
                                variant_mismatch_avg=variant_mismatch_avg, variant_dist3effective_avg=variant_dist3effective_avg,
								num_cosmic=num_cosmic, variant_effect=variant_effect,
								nucleotid_change=nucleotid_change, variant_confidence=variant_confidence,
								norm_zygosity=norm_zygosity, norm_reference_count=norm_reference_count, norm_variant_count=norm_variant_count,
								tumor_zygosity=tumor_zygosity, tumor_reference_count=tumor_reference_count, tumor_variant_count=tumor_variant_count,
								details=details)
		add_measurement!(d, gene, exome_data)
	end
end

"""
Populates a data view with Methylation data from JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene, Methylation}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)

        # TODO: decide for key name
		# beta_value = entry["illumina_beta_value"]
        beta_value = entry["beta_value"]
		#optionals

        cgct1_value = get(entry, "cgct1", 1)
        cct1_value = get(entry, "cct1", 1)
		methylated_threshold = get(entry, "methylated_threshold", .2)
		illumina_id = get(entry, "illumina_id", "")

		methylation_data = Methylation(beta_value, cgct1_value, cct1_value, methylated_threshold, illumina_id=illumina_id)
		add_measurement!(d, gene, methylation_data)
	end
end

"""
Populate a data view with gene expression data from JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene, GeneExpression}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)

		expression_value = Float64(entry["expression_value"])

		add_measurement!(d, gene, GeneExpression(expression_value))
	end
end

"""
Populate a data view with RNASeq data from JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene, RNASeq}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)

		expression_value = entry["expression_value"]

		add_measurement!(d, gene, RNASeq(expression_value, expression_status))
	end
end

"""
Populate a data view with RNASeqCall data from JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene, RNASeqCall}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)

		expression_status = entry["expression_status"]

		add_measurement!(d, gene, RNASeqCall(expression_value, expression_status))
	end
end

"""
Populate a data view with RPPA data from JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Protein, RPPA}, data::Vector{Any})
	for entry in data
		protein = add_protein(experiment_id, entry)

		abundance = entry["protein_abundance"]

		add_measurement!(d, protein, RPPA(abundance))
	end
end

"""
Populate a data view with CNV data form JSON.
"""
function populate_data_view!(experiment_id::String, d::DataView{Gene, CNV}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment_id, entry)
		cnv_value = Float64(entry["gene_level_cnv"])
		add_measurement!(d, gene, CNV(cnv_value))
	end
end
