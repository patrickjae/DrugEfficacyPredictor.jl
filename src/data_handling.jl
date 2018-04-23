######### EXPERIMENT ###########

# create the overall data set
create_experiment() = Experiment()


######### CELL LINES ###########
# get a cell line object for an experiment or creates a new one
get_cell_line(experiment::Experiment, cell_line_id::String, cancer_type::String) = get!(experiment.cell_lines, cell_line_id, CellLine(cell_line_id, cancer_type))

function get_cell_line(experiment::Experiment, data::Dict{String, Any})
	if !haskey(data, "id") || !haskey(data, "cancer_type")
		throw(ArgumentError("You need to provide a cell line id and the cancer type."))
	end
	get_cell_line(experiment, data["id"], data["cancer_type"])
end
######### VIEWS ###########
add_view!(experiment::Experiment, view::Type{<:ViewType}) = union!(experiment.views, [view])


function normalize_data_views(experiment::Experiment)
	for v in experiment.views
		all_values = Float64[]
		#collect all values
		for cl in values(experiment.cell_lines)
			if haskey(cl.views, v)
				all_values.push!(all_values, map(get_measurement_value, collect(values(cl.views[v]))))
			end
		end
		view_mean = mean(all_values)
		view_std = stdm(all_values, view_mean)
		for cl in values(experiment.cell_lines)
			if haskey(cl.views, v)
				f = (x) -> (x-view_mean)/view_std
				measurements = collect(values(cl.views[v]))
				norm_vals = map(f, map(get_measurement_value, measurements))
				map((m, nv) -> set_normalized_value(m, nv), zip(measurements, norm_vals))
			end
		end
	end
end

######### GENES ###########
function add_gene(experiment::Experiment, gene::Gene) 
	experiment.genes[gene.entrez_id] = gene
	experiment.genes_by_hgnc[gene.hgnc_id] = gene
	experiment.genes_by_ensembl[gene.ensembl_id] = gene
	gene
end

function add_gene(experiment::Experiment, gene_id::Union{Int64, String}, id_type::String="entrez_id"; is_cancer_gene::Bool=false)
	try
		return get_gene(experiment, gene_id, id_type)
	catch KeyError
		gene_object = Gene()
		gene_object.cancer_gene = is_cancer_gene
		add_gene_id!(experiment, gene_object, gene_id, id_type)
		return gene_object
	end
end

function add_gene(experiment::Experiment, data::Dict{String, Any})
	display(data)
	if !haskey(data, "gene_id")
		throw(ArgumentError("You need to provide a gene id."))
	end
	gene_id = data["gene_id"]
	id_type = get(data, "type_id", "entrez_id")
	is_cancer_gene = get(data, "is_cancer_gene", false)
	add_gene(experiment, gene_id, id_type, is_cancer_gene = is_cancer_gene)
end
function add_gene(experiment::Experiment, gene_id::Int64; is_cancer_gene::Bool=false)
	get!(experiment.genes, gene_id, Gene(gene_id, is_cancer_gene))
end

function add_gene_id!(experiment::Experiment, gene::Gene, id::String, id_type::String)
	if id_type == "hgnc_id"
		if (gene.hgnc_id != "") && (gene.hgnc_id != id) warn("overwriting HGNC_ID, new value is $(id), was $(gene.hgnc_id) before") end
		gene.hgnc_id = id
		experiment.genes_by_hgnc[id] = gene
	elseif id_type == "ensembl_id"
		if (gene.ensembl_id != "") && (gene.ensembl_id != id) warn("overwriting Ensembl ID, new value is $(id), was $(gene.ensembl_id) before") end
		gene.ensembl_id = id
		experiment.genes_by_ensembl[id] = gene
	else
		warn("unknown genen id type $id_type, ignoring")
	end
end

function add_gene_id!(experiment::Experiment, gene::Gene, entrez_id::Int64, id_type::String)
	gene.entrez_id = entrez_id
	experiment.genes[entrez_id] = gene
end

# get a gene object
function get_gene(experiment::Experiment, gene_id::String, id_type::String="hgnc_id") 
	if id_type == "hgnc_id"
		return experiment.genes_by_hgnc[gene_id]
	elseif id_type == "ensembl_id"
		return experiment.genes_by_ensembl[gene_id]
	end
end

get_gene(experiment::Experiment, gene_id::Int64, id_type::String="entrez_id") = experiment.genes[gene_id]


######### PROTEINS ###########
function add_protein(experiment::Experiment, hgnc_id::String; fully_validated::Bool=false)
	p_obj = Protein(hgnc_id, fully_validated)
	experiment.proteins[hgnc_id] = p_obj
	p_obj
end

get_protein(experiment::Experiment, hgnc_id::String; fully_validated::Bool=false) = get!(experiment.proteins, hgnc_id, Protein(hgnc_id, fully_validated))



######### DATA VIEWS ###########
# add a data view to the cell line
add_dataview!(cl::CellLine, d::DataView{T,D}) where {T,D} = cl.views[D] = d
add_dataview(cl::CellLine, d::DataView{T,D}) where {T,D}  = cl.views[D] = d

# get a data view for a certain data type from a cell line or create if not present
get_dataview(cl::CellLine, data_type::Type) = cl.views[data_type]


######### MEASUREMENTS ###########
# add a measurement, i.e. a gene/protein data pair to a data view
function add_measurement(d::DataView{T,D}, subject::T, measurement::D) where {T,D}
	try
		measurement = get_measurement(d, subject)
	catch KeyError
		d.measurements[subject] = measurement
		return d
	end
	warn("measurement already exists for $subject, ignoring...")
	warn("if you want to force the new value, use add_measurement! instead")
	d
end

function add_measurement!(d::DataView{T,D}, subject::T, measurement::D) where {T,D}
	d.measurements[subject] = measurement
	d
end

get_measurement(d::DataView{T,<:ViewType}, subject::T) where {T} = d.measurements[subject]
get_measurement_value(d::DataView{T, <:ViewType}, subject::T) where {T} = get_measurement_value(d.measurements[subject])

function get_measurement_value(d::ExomeSeq)
	(d.reference_mismatch_avg + d.variant_mismatch_avg + d.reference_mismatch_sum 
		+ d.variant_mismatch_sum + d.reference_dist3effective_avg + d.variant_dist3effective_avg)
end

get_measurement_value(d::GeneExpression) = d.expression_value
get_measurement_value(d::Methylation) = d.illumina_beta_value
get_measurement_value(d::RNASeq) = (d.expression_value, d.expression_status)
get_measurement_value(d::RPPA) = d.protein_abundance
get_measurement_value(d::CNV) = d.gene_level_cnv


######### NORMALIZATION ###########
function set_normalized_value(d::T, value::Float64) where {T <: ViewType}
	d.normalized_value = value
end

get_normalized_value(d::DataView{T, D}, subject::T) where {T, D <: ViewType} = d.measurements[subject].normalized_value
get_normalized_value(d::T) where {T <: ViewType} = d.normalized_value



######### RESULTS ###########
# add results for a drug on a cell line to an outcome
add_result!(o::Outcome, cl::CellLine, result::Float64) = o.outcome_values[cl] = result

# add results for a given array of results
add_results!(o::Outcome, results::Dict{CellLine, Float64}) = o.outcome_values = results

function add_results!(o::Outcome, cell_lines::Vector{CellLine}, results::Vector{Float64})
	if length(cell_lines) != length(results)
		throw(DimensionMismatch("You must provide a matching number of results. Got $(length(cell_lines)) cell lines but $(length(results)) results."))
	end
	for (cl,res) in zip(cell_lines, results)
		add_result!(o, cl, res)
	end
end


######### OUTCOME ###########
# add an outcome of a drug to the experiment
function add_outcome!(e::Experiment, d::Drug, o::Outcome)
	e.results[d] = o
	e
end
