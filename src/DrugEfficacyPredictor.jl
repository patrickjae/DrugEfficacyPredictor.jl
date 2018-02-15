module DrugEfficacyPredictor
include("types.jl")

########################
######### TODO #########
########################

# steps

# populate model

# expose API

# cross-validation during training

# split application and API

########################
######### TODO #########
########################


# create the overall data set
create_experiment() = Experiment()

# get a cell line object for an experiment
get_cell_line(experiment::Experiment, cell_line_id::String, cancer_type::String) = get!(experiment.cell_lines, cell_line_id, CellLine(cell_line_id, cancer_type))

# add a data view to the experiment
add_dataview!{T,D}(e::Experiment{T,D}, d::DataView{T,D}) = union!(e.measurements, d)
add_dataview{T,D}(e::Experiment{T,D}, d::DataView{T,D}) = union(e.measurements, d)

# add an empty data view to the experiment, associated with a cell line
add_dataview!(e::Experiment, cl::CellLine) = union!(e.measurements, DataView{T,D}(cl))
add_dataview(e::Experiment, cl::CellLine) = union(e.measurements, DataView{T,D}(cl))

# add a measurement, i.e. a gene/protein data pair to a data view
function add_measurement!{T,D}(d::DataView{T,D}, subject::T, measurement::D)
	d.data[subject] = measurement
	d
end

# add an outcome of a drug to the experiment
function add_outcome!(e::Experiment, d::Drug, o::Outcome)
	e.results[d] = o
	e
end

# add results for a drug on a cell line to an outcome
add_result!(o::Outcome, c::CellLine, result::Float64) = o.outcome_values[cl] = result

# add results for a given array of results
add_results!(o::Outcome, results::Dict{CellLine, Float64}) = o.outcome_values = results
function add_results!(o::Outcome, cell_lines::Vector{CellLine}, results::Vector{Float64})
	if length(cell_lines) != length(results)
		throw(DimensionMismatch("You must provide a matching number of results. Got $(length(cell_lines)) cell lines but $(length(results)) results"))
	end
	for (cl,res) in zip(cell_lines, results)
		add_result!(o, cl, res)
	end
end

##################### I/O #####################

function import_dream_challenge_data(directory::String)
	experiment = create_experiment()
	for f in readdir(directory)
		if !endswith(f,".txt") || endswith(f, "README.txt")
			continue
		end
		data_type_string = replace(f, "DREAM7_DrugSensitivity1_", "")
		key_type = Gene
		if data_type_string == "Exomeseq"
			data_type = ExomeSeq
		elseif data_type_string == "GeneExpression"
			data_type = GeneExpression
		elseif data_type_string == "Methylation"
			data_type = Methylation
		elseif data_type_string == "RNAseq_expressed_calls" || data_type_string == "RNAseq_quantification"
			data_type = RNAseq
		elseif data_type_string == "RPPA"
			key_type = Protein
			data_type = RPPA
		elseif data_type_string == "SNP6"
			data_type = CNV
		end
		process_data_frame(df, key_type, data_type, experiment)
	end
end

function process_data_frame(f::String, gene_type::Gene, exome_data::ExomeSeq, experiment::Experiment)
	type_vector = Vector{DataType}(27)
	fill!(type_vector, String)
	df = CSV.read(f, DataFrame, delim='\t', types=type_vector, null=".")
	all_cell_lines = unique(df[:CellLine])
	for cl in all_cell_lines
		# get the cell line object somehow, possibly query an experiment-wide dictionary
		cl_obj = CellLine(cl, "Breast cancer")

		#create the data view, unknown key type
		key_type = #combination of gene and protein change
		data_view = DataView{key_type, ExomeSeq}(cl_obj)

		cl_df = view(df, df[:CellLine] .== cl)
		for hgnc_id in unique(cl_df[:HGNC_ID])
			gene_view = view(cl_df, cl_df[:HGNC_ID] .== hgnc_id)
			for data_row in eachrow(gene_view)
				num_cosmic = data_row[Symbol("#COSMIC")]
				variant_effect = data_row[:Type]
				protein_change = data_row[:Summary]
				nucleotid_change = data_row[:AltBase]
				variant_confidence = data_row[:Confidence]
				norm_zygosity = data_row[Symbol("Zygosity(norm)")]
				norm_reference_count = data_row[Symbol("RefCount(norm)")]
				norm_variant_count = data_row[Symbol("AltCount(norm)")]
				tumor_zygosity = data_row[Symbol("Zygosity(tumor)")]
				tumor_reference_count = data_row[Symbol("RefCount(tumor)")]
				tumor_variant_count = data_row[Symbol("AltCount(tumor)")]
				reference_mismatch_avg = data_row[Symbol("Avg#Mismatch(ref)")]
				variant_mismatch_avg = data_row[Symbol("Avg#Mismatch(alt)")]
				reference_mismatch_sum = data_row[Symbol("MismatchQualitySum(ref)")]
				variant_mismatch_sum = data_row[Symbol("MismatchQualitySum(alt)")]
				reference_dist3effective_avg = data_row[Symbol("DistanceEffective3'End(ref)")]
				variant_dist3effective_avg = data_row[Symbol("DistanceEffective3'End(alt)")]
				details = data_row[:Details]

				#create some exome data view, create a key type for exome sequencing as a combination
				#of gene and protein change
				exome_data = (...)

				add_measurement(data_view, key, exome_data)
			end
		end
	end
	experiment
end

end # module
