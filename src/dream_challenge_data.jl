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
			# problem: if not Type is different from Missense, the summary does not exist
			key_type = Tuple{Gene, String}
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
		elseif data_type_string != "RPPA"
			# all_genes = df[:HGNC_ID]
			is_cancer_gene = false
			i = 1
			n = size(df, 1)
			for data_row in eachrow(df)
				gene_object = add_gene(experiment, data_row[:HGNC_ID].value, is_cancer_gene = is_cancer_gene)
				if in(Symbol("CancerGene?"), names(df))
					is_cancer_gene = data_row[Symbol("CancerGene?")].value == "No" ? false : true
					gene_object.cancer_gene = is_cancer_gene
				end
				if in(:Ensembl_ID, names(df))
					e_id = data_row[:Ensembl_ID].value
					add_gene_id!(experiment, gene_object, e_id, "ensembl_id")
				end
				if in(:EntrezID, names(df))
					e_id = data_row[:EntrezID].value
					add_gene_id!(experiment, gene_object, e_id)
				end
				i+=1
				if i%5000 == 0 info("processed $i of $n genes...") end
			end
		else
			i = 1
			all_proteins = df[[:Antibody_ID, :FullyValidated]]
			n = size(all_proteins, 1)
			for p_row in eachrow(all_proteins)
				is_fully_validated = p_row[:FullyValidated].value == "Yes" ? true : false
				add_protein(experiment, p_row[:Antibody_ID].value, fully_validated = is_fully_validated)
			end
			i += 1
			if i%5000 == 0 info("processed $i of $n proteins...") end
		end
		info("importing cell line data")
		for cl in get_cell_line_names_from_data_frame(data_type, df)
			cl_obj = get_cell_line(experiment, cl, "Breast cancer")
			data_view = nothing
			try
				data_view = get_dataview(cl_obj, data_type)
			catch KeyError
				data_view = DataView{key_type, data_type}(cl)
			end
			populate_data_view!(data_view, df, experiment)

			add_dataview!(cl_obj, data_view)
		end
		add_view!(experiment, data_type)
	end
	experiment
end

get_cell_line_names_from_data_frame(::Type{T}, df::DataFrame) where T<:ExomeSeq = map(string, map(v->v.value, unique(df[:CellLine])))
get_cell_line_names_from_data_frame(::Type{T}, df::DataFrame) where T<:GeneExpression = extract_column_names(df, 2)
get_cell_line_names_from_data_frame(::Type{T}, df::DataFrame) where T<:Methylation = extract_column_names(df, 5)
get_cell_line_names_from_data_frame(::Type{T}, df::DataFrame) where T<:Union{RNASeq, RPPA, CNV} = extract_column_names(df, 3)

extract_column_names(df::DataFrame, first_used_column::Int64) = map(string, names(df)[first_used_column:end])

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Tuple{Gene, String}, T<:ExomeSeq}
	cl_df = view(df, df[:CellLine] .== data_view.cell_line_id)
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
		key = (gene, protein_change)
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

		add_measurement!(data_view, key, exome_data)
	end
end

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Gene, T <: GeneExpression}
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		if !data_row[2].hasvalue continue end
		gene = get_gene(experiment, data_row[1].value)
		gene_expression = GeneExpression(data_row[2].value)
		add_measurement!(data_view, gene, gene_expression)
	end
end

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Gene, T <: Methylation}
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id), :CGct1, :Cct1, :Illumina_ID]]
	for data_row in eachrow(cl_df)
		gene = get_gene(experiment, data_row[1].value)
		methylation = Methylation(data_row[2].value, data_row[3].value, data_row[4].value, illumina_id = data_row[5].value)
		add_measurement!(data_view, gene, methylation)
	end
end

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Gene, T <: RNASeq}
	tic()
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		gene = get_gene(experiment, data_row[1].value)
		is_quant = eltype(data_row[2].value) == Float64 ? false : true
		try 
			measurement = get_measurement(data_view, gene)
			if is_quant
				measurement.expression_status = data_row[2].value == 0 ? false : true
			else
				measurement.expression_value = data_row[2].value
			end
		catch KeyError
			if is_quant
				add_measurement!(data_view, gene, RNASeq(0.0, true))
			else
				add_measurement!(data_view, gene, RNASeq(data_row[2].value, false))
			end
		end
	end
	info("cell line $(data_view.cell_line_id) took $(toq()) seconds")
end

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Protein, T <: RPPA}
	cl_df = df[[:Antibody_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		rppa = RPPA(data_row[2].value)
		protein = get_protein(experiment, data_row[1].value)
		add_measurement!(data_view, protein, rppa)
	end
end

function populate_data_view!(data_view::DataView{K,T}, df::DataFrame, experiment::Experiment) where {K <: Gene, T <: CNV}
	cl_df = df[[:HGNC_ID, Symbol(data_view.cell_line_id)]]
	for data_row in eachrow(cl_df)
		if !data_row[2].hasvalue continue end
		gene_level_cnv = CNV(data_row[2].value)
		gene = get_gene(experiment, data_row[1].value)
		add_measurement!(data_view, gene, gene_level_cnv)
	end
end