""" 
Creates a Cell line object from JSON. Creates Gene (or Protein) objects on the fly as needed.
"""
function import_cell_line(experiment::Experiment, data::Dict{String, Any})
	info("fetching/creating cell line...")
	cl = get_cell_line(experiment, data)
	if haskey(data, "views")
		# construct the views
		views = data["views"]
		# should be a dict of views
		for (k,v) in views
			view_type = eval(parse(k))
			if view_type == ExomeSeq
				key_type = Tuple{Gene, String}
			elseif view_type == RPPA
				key_type = Protein
			else
				key_type = Gene
			end
			data_view = get_dataview!(cl_obj, view_type, DataView{key_type, view_type}(cl.id))
			info("populating data view with key type '$key_type' and view type '$view_type'")
			populate_data_view!(experiment, data_view, v)
			info("adding view type to experiment")
			add_view!(experiment, view_type)
		end

	end

end


"""
Populates a data view with ExomeSeq data from JSON data.
If no protein change is provided, we assume the empty string.
"""
function populate_data_view!(experiment::Experiment, d::DataView{Gene,ExomeSeq}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment, entry)

		# mandatory values
		if !(haskey(entry,"reference_mismatch_avg") && haskey(entry,"variant_mismatch_avg") 
				&& haskey(entry,"reference_mismatch_sum") && haskey(entry, "variant_mismatch_sum")
				&& haskey(entry, "reference_dist3effective_avg") && haskey(entry, "variant_dist3effective_avg"))
			throw(ArgumentError("For ExomeSeq data, you have to specify at least values for reference_mismatch_avg, variant_mismatch_avg,
									reference_mismatch_sum, variant_mismatch_sum, reference_dist3effective_avg and variant_dist3effective_avg!"))
		end
		reference_mismatch_avg = entry["reference_mismatch_avg"]
		reference_mismatch_sum = entry["reference_mismatch_sum"]
		reference_dist3effective_avg = entry["reference_dist3effective_avg"]
		variant_mismatch_avg = entry["variant_mismatch_avg"]
		variant_mismatch_sum = entry["variant_mismatch_sum"]
		variant_dist3effective_avg = entry["variant_dist3effective_avg"]

		# optional values
		protein_change = get(entry, "protein_change", "")

		num_cosmic = get(entry, "num_cosmic", 0)
		is_cancer_gene = get(entry, "is_cancer_gene", false)
		variant_effect = get(entry, "variant_effect", "")
		nucleotid_change = get(entry, "nucleotid_change", "")
		variant_confidence = get(entry, "variant_confidence", 1.)
		norm_zygosity = get(entry, "norm_zygosity", "")
		norm_reference_count = get(entry, "norm_reference_count", 0)
		norm_variant_count = get(entry, "norm_variant_count", 0)
		tumor_zygosity = get(entry, "tumor_zygosity", "")
		tumor_reference_count = get(entry, "tumor_reference_count", 0)
		tumor_variant_count = get(entry, "tumor_variant_count", 0)
		details = get(entry, "details", "")

		exome_data = ExomeSeq(reference_mismatch_sum, reference_mismatch_avg, reference_dist3effective_avg,
								variant_mismatch_sum, variant_mismatch_avg, variant_dist3effective_avg,
								num_cosmic=num_cosmic, variant_effect=variant_effect, protein_change=protein_change,
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
function populate_data_view!(experiment::Experiment, d::DataView{Gene, Methylation}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment, entry)

		cgct1_value = entry["cgct1"]
		cct1_value = entry["cct1"]
		beta_value = entry["illumina_beta_value"]

		#optionals

		methylated_threshold = get(entry, "methylated_threshold", .2)
		illumina_id = get(entry, "illumina_id", "")

		methylation_data = Methylation(beta_value, cgct1_value, cct1_value, methylated_threshold, illumina_id=illumina_id)
		add_measurement!(d, gene, methylation_data)
	end
end

"""
Populate a data view with gene expression data from JSON.
"""
function populate_data_view!(experiment::Experiment, d::DataView{Gene, GeneExpression}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment, entry)

		expression_value = entry["expression_value"]

		add_measurement!(d, gene, GeneExpression(expression_value))
	end
end

"""
Populate a data view with RNASeq data from JSON.
"""
function populate_data_view!(experiment::Experiment, d::DataView{Gene, RNASeq}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment, entry)

		expression_value = entry["expression_value"]
		expression_status = entry["expression_status"]

		add_measurement!(d, gene, RNASeq(expression_value, expression_status))
	end
end

"""
Populate a data view with RPPA data from JSON.
"""
function populate_data_view!(experiment::Experiment, d::DataView{Protein, RPPA}, data::Vector{Any})
	for entry in data
		protein = add_protein(experiment, entry)

		abundance = entry["protein_abundance"]

		add_measurement!(d, protein, RPPA(abundance))
	end
end

"""
Populate a data view with CNV data form JSON.
"""
function populate_data_view!(experiment::Experiment, d::DataView{Gene, CNV}, data::Vector{Any})
	for entry in data
		gene = add_gene(experiment, data)
		cnv_value = entry["gene_level_cnv"]
		add_measurement!(d, gene, CNV(cnv_value))
	end
end