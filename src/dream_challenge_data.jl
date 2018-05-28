function import_dream_challenge_data(directory::String)
	info("importing dream challenge data from $directory")
	experiment = create_experiment()
	for f in readdir(directory)
		if !endswith(f,".txt") || endswith(f, "README.txt")
			continue
		end
		data_type_string = replace(f, "DREAM7_DrugSensitivity1_", "")
		data_type_string = replace(data_type_string, ".txt","")
		key_type = Gene
		if data_type_string == "Exomeseq"
			data_type = ExomeSeq
		elseif data_type_string == "GeneExpression"
			data_type = GeneExpression
		elseif data_type_string == "Methylation"
			data_type = Methylation
		elseif data_type_string == "RNAseq_expressed_calls" || data_type_string == "RNAseq_quantification"
			data_type = RNASeq
		elseif data_type_string == "RPPA"
			key_type = Protein
			data_type = RPPA
		elseif data_type_string == "SNP6_gene_level"
			data_type = CNV
		elseif data_type_string == "Drug_Response_Training"
			# responses
		else
			warn("couldn't determine data type: $data_type_string")
			continue
		end
		info("importing file $f")
		# make nullable
		df = CSV.read(joinpath(directory, f), DataFrame, delim='\t',rows_for_type_detect=1200, null="NA", weakrefstrings=false)
		if data_type_string == "Drug_Response_Training"
			df = CSV.read(joinpath(directory, f), DataFrame, delim='\t',rows_for_type_detect=1200, null="NA", weakrefstrings=false)
		end
		info("creating genes or proteins...")
		info("number of genes so far: Entrez: $(length(experiment.genes)), HGNC: $(length(experiment.genes_by_hgnc)), Ensembl: $(length(experiment.genes_by_ensembl))")
		if data_type_string == "Drug_Response_Training"
			drugs = names(df)[2:end]
			cell_lines = df[:CellLine].values
			cl_objects = map(cl_id -> get_cell_line(experiment, cl_id, "Breast cancer"), cell_lines)
			for d in drugs
				vals = df[d].values
				o = Outcome("IC50")
				add_results!(o, cl_objects[.!df[d].isnull], vals[.!df[d].isnull])
				add_outcome!(experiment, Drug(string(d)), o)
			end
			continue
		elseif data_type_string == "RPPA"
			i = 1
			all_proteins = df[[:Antibody_ID, :FullyValidated]]
			n = size(all_proteins, 1)
			for p_row in eachrow(all_proteins)
				is_fully_validated = p_row[:FullyValidated].value == "Yes" ? true : false
				add_protein(experiment, p_row[:Antibody_ID].value, fully_validated = is_fully_validated)
			end
			i += 1
			if i%5000 == 0 info("processed $i of $n proteins...") end
		else
			i = 1
			n = size(df, 1)
			df_names = names(df)
			for data_row in eachrow(df)
				gene_object = nothing
				has_entrez_id = false
				has_hgnc_id = false
				has_ensembl_id = false
				entrez_id = -1
				ensembl_id = ""
				hgnc_id = ""

				# check for different id types and try to retrieve gene object
				if in(:EntrezID, df_names)
					has_entrez_id = true
					entrez_id = data_row[:EntrezID].value
					try 
						gene_object = get_gene(experiment, entrez_id)
					end
				end

				if in(:Ensembl_ID, df_names)
					has_ensembl_id = true
					ensembl_id = data_row[:Ensembl_ID].value
					if gene_object == nothing
						try
							gene_object = get_gene(experiment, ensembl_id, "ensembl_id")
						end
					else
						add_gene_id!(experiment, gene_object, ensembl_id, "ensembl_id")
					end
				end

				if in(:HGNC_ID, df_names)
					has_hgnc_id = true
					hgnc_id = data_row[:HGNC_ID].value
					if gene_object == nothing
						try
							gene_object = get_gene(experiment, hgnc_id, "hgnc_id")
						end
					else
						add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id")
					end
				end

				if has_entrez_id
					if gene_object == nothing
						#no gene found for any id type
						gene_object = add_gene(experiment, entrez_id)
						if has_ensembl_id add_gene_id!(experiment, gene_object, ensembl_id, "ensembl_id") end
						if has_hgnc_id add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id") end
					else
						# check if object had no entrez id before
						if gene_object.entrez_id == -1 
							add_gene_id!(experiment, gene_object, entrez_id, "entrez_id")
						else
							# found gene object has another entrez id (which is unique) so create a new gene
							gene_object = add_gene(experiment, entrez_id)
							if has_ensembl_id add_gene_id!(experiment, gene_object, ensembl_id, "ensembl_id") end
							if has_hgnc_id add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id") end
						end
					end
				end

				if has_ensembl_id
					if gene_object == nothing
						# if still no gene object, no entrez id is present
						gene_object = add_gene(experiment, ensembl_id, "ensembl_id")
						if has_hgnc_id add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id") end
					else
						# ensembl ids also are uniqe
						if gene_object.ensembl_id == ""
							add_gene_id!(experiment, gene_object, ensembl_id, "ensembl_id")
						else
							gene_object = add_gene(experiment, ensembl_id, "ensembl_id")
							if has_entrez_id add_gene_id!(experiment, gene_object, entrez_id, "entrez_id") end
							if has_hgnc_id add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id") end
						end
					end
				end

				if has_hgnc_id
					if gene_object == nothing
						# if still no gene object, hgnc id is the only id present
						gene_object = add_gene(experiment, hgnc_id, "hgnc_id")
					else
						if gene_object.hgnc_id != hgnc_id add_gene_id!(experiment, gene_object, hgnc_id, "hgnc_id") end
					end
				end

				if in(Symbol("CancerGene?"), df_names)
					is_cancer_gene = data_row[Symbol("CancerGene?")].value == "No" ? false : true
					gene_object.cancer_gene = is_cancer_gene
				end
				i+=1
				if i%5000 == 0 info("processed $i of $n genes...") end
			end
		end
		info("number of genes now: Entrez: $(length(experiment.genes)), HGNC: $(length(experiment.genes_by_hgnc)), Ensembl: $(length(experiment.genes_by_ensembl))")
		info("importing cell line data")
		for cl in get_cell_line_names_from_data_frame(data_type, df)
			cl_obj = get_cell_line(experiment, cl, "Breast cancer")
			data_view = get_dataview!(cl_obj, data_type, DataView{key_type, data_type}(cl))
			populate_data_view!(data_view, df, experiment)
		end
		add_view!(experiment, data_type)
	end
	experiment
end

get_cell_line_names_from_data_frame(::Type{ExomeSeq}, df::DataFrame) = map(string, map(v->v.value, unique(df[:CellLine])))
get_cell_line_names_from_data_frame(::Type{GeneExpression}, df::DataFrame) = extract_column_names(df, 2)
get_cell_line_names_from_data_frame(::Type{Methylation}, df::DataFrame) = extract_column_names(df, 5)
get_cell_line_names_from_data_frame(::Type{<:Union{RNASeq, RPPA, CNV}}, df::DataFrame) = extract_column_names(df, 3)

extract_column_names(df::DataFrame, first_used_column::Int64) = map(string, names(df)[first_used_column:end])

function populate_data_view!(data_view::DataView{Gene,ExomeSeq}, df::DataFrame, experiment::Experiment)
	cl_df = view(df, df[:CellLine].values .== data_view.cell_line_id)
	for data_row in eachrow(cl_df)
		gene = get_gene(experiment, data_row[:HGNC_ID].value)
		num_cosmic = data_row[Symbol("#Cosmic")].value
		is_cancer_gene = data_row[Symbol("CancerGene?")].value == "No" ? false : true
		variant_effect = data_row[:Type].value
		protein_change = data_row[:Summary].hasvalue ? data_row[:Summary].value : "None"
		nucleotid_change = data_row[:AltBase].value
		variant_confidence = Float64(data_row[:Confidence].value)
		norm_zygosity = data_row[Symbol("Zygosity(norm)")].value
		norm_reference_count = data_row[Symbol("RefCount(norm)")].value == "." ? 0 : parse(Int64, data_row[Symbol("RefCount(norm)")].value)
		norm_variant_count = data_row[Symbol("altCount(norm)")].value == "." ? 0 : parse(Int64, data_row[Symbol("altCount(norm)")].value)
		tumor_zygosity = data_row[Symbol("Zygosity(tumor)")].value
		tumor_reference_count = data_row[Symbol("RefCount(tumor)")].value
		tumor_variant_count = data_row[Symbol("AltCount(tumor)")].value
		reference_mismatch_avg = data_row[Symbol("Avg#Mismatch(ref)")].value
		variant_mismatch_avg = data_row[Symbol("Avg#Mismatch(alt)")].value
		reference_mismatch_sum = data_row[Symbol("MismatchQualitySum(ref)")].value
		variant_mismatch_sum = data_row[Symbol("MismatchQualitySum(alt)")].value
		reference_dist3effective_avg = data_row[Symbol("DistanceEffective3'end(ref)")].value
		variant_dist3effective_avg = data_row[Symbol("DistanceEffective3'end(alt)")].value
		details = data_row[:Details].value
		# data_summary = reference_mismatch_avg + variant_mismatch_avg + reference_mismatch_sum + variant_mismatch_sum + reference_dist3effective_avg + variant_dist3effective_avg
		
		#create some exome data view, create a key type for exome sequencing as a combination
		#of gene and protein change
		# exome_data = ExomeSeq(num_cosmic, variant_effect, protein_change,
		# 		nucleotid_change, variant_confidence, 
		# 		norm_zygosity, norm_reference_count, norm_variant_count,
		# 		tumor_zygosity, tumor_reference_count, tumor_variant_count,
		# 		reference_mismatch_avg, variant_mismatch_avg,
		# 		reference_mismatch_sum, variant_mismatch_sum,
		# 		reference_dist3effective_avg, variant_dist3effective_avg,
		# 		details, data_summary)
		exome_data = ExomeSeq(reference_mismatch_sum, reference_mismatch_avg, reference_dist3effective_avg,
								variant_mismatch_sum, variant_mismatch_avg, variant_dist3effective_avg,
								num_cosmic=num_cosmic, variant_effect=variant_effect, protein_change=protein_change,
								nucleotid_change=nucleotid_change, variant_confidence=variant_confidence,
								norm_zygosity=norm_zygosity, norm_reference_count=norm_reference_count, norm_variant_count=norm_variant_count,
								tumor_zygosity=tumor_zygosity, tumor_reference_count=tumor_reference_count, tumor_variant_count=tumor_variant_count,
								details=details)
		add_measurement!(data_view, gene, exome_data)
	end
end

function populate_data_view!(data_view::DataView{Gene,GeneExpression}, df::DataFrame, experiment::Experiment)
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		if !data_row[2].hasvalue continue end
		gene = get_gene(experiment, data_row[1].value)
		gene_expression = GeneExpression(data_row[2].value)
		add_measurement!(data_view, gene, gene_expression)
	end
end

function populate_data_view!(data_view::DataView{Gene,Methylation}, df::DataFrame, experiment::Experiment)
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id), :CGct1, :Cct1, :Illumina_ID]]
	for data_row in eachrow(cl_df)
		gene = get_gene(experiment, data_row[1].value)
		methylation = Methylation(data_row[2].value, data_row[3].value, data_row[4].value, illumina_id = data_row[5].value)
		add_measurement!(data_view, gene, methylation)
	end
end

function populate_data_view!(data_view::DataView{Gene,RNASeq}, df::DataFrame, experiment::Experiment)
	cl_df = df[[:Ensembl_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		gene = get_gene(experiment, data_row[1].value, "ensembl_id")
		is_quant = eltype(data_row[2].value) == Float64 ? false : true
		empty_measurement = RNASeq(0.0, false)
		# this will return the existing measurment for the gene and spit out a warn message
		measurement = get_measurement!(data_view, gene, empty_measurement)
		if is_quant
			measurement.expression_status = Bool(data_row[2].value)
		else
			measurement.expression_value = data_row[2].value
		end
	end
end

function populate_data_view!(data_view::DataView{Protein,RPPA}, df::DataFrame, experiment::Experiment)
	cl_df = df[[:Antibody_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		rppa = RPPA(data_row[2].value)
		protein = get_protein(experiment, data_row[1].value)
		add_measurement!(data_view, protein, rppa)
	end
end

function populate_data_view!(data_view::DataView{Gene,CNV}, df::DataFrame, experiment::Experiment)
	cl_df = df[[:EntrezID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		if !data_row[2].hasvalue continue end
		gene_level_cnv = CNV(data_row[2].value)
		gene = get_gene(experiment, data_row[1].value)
		add_measurement!(data_view, gene, gene_level_cnv)
	end
end
