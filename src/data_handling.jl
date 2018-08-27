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
	if !haskey(data, "gene_id")
		throw(ArgumentError("You need to provide a gene id."))
	end
	if data["gene_id"] == nothing
		throw(ArgumentError("The provided gene id is null."))
	end
	gene_id = data["gene_id"]
	id_type = get(data, "type_id", "entrez_id")
	is_cancer_gene = get(data, "is_cancer_gene", false)
	# info("gene_id: $gene_id, id_type: $id_type, cancer: $is_cancer_gene")
	add_gene(experiment, gene_id, id_type, is_cancer_gene = is_cancer_gene)
end

function add_gene(experiment::Experiment, gene_id::Int64; is_cancer_gene::Bool=false)
	get!(experiment.genes, gene_id, Gene(gene_id, is_cancer_gene))
end

function add_gene_id!(experiment::Experiment, gene::Gene, id::String, id_type::String)
	if id_type == "hgnc_id"
		if (gene.hgnc_id != "") && (gene.hgnc_id != id) 
			warn("overwriting HGNC_ID, new value is $(id), was $(gene.hgnc_id) before") 
		end
		gene.hgnc_id = id
		experiment.genes_by_hgnc[id] = gene
	elseif id_type == "ensembl_id"
		if (gene.ensembl_id != "") && (gene.ensembl_id != id) 
			warn("overwriting Ensembl ID, new value is $(id), was $(gene.ensembl_id) before (other IDs: $(gene.hgnc_id) (HGNC), $(gene.entrez_id) (Entrez))") 
		end
		gene.ensembl_id = id
		experiment.genes_by_ensembl[id] = gene
	else
		warn("unknown gene id type $id_type, ignoring")
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
	throw(KeyError("key \"$gene_id\" not found for id_type \"$id_type\""))
end

get_gene(experiment::Experiment, gene_id::Int64, id_type::String="entrez_id") = experiment.genes[gene_id]

######### PATHWAYS ###########
function add_pathway(experiment::Experiment, genes::Vector{K}=Vector{Gene}(); id::String="", name::String="") where K<:KeyType
	pw = Pathway{K}(genes, name = name, id = id)
	push!(experiment.pathway_information, pw)
	pw
end

function add_pathway(experiment::Experiment, data::Dict{String, Any})
	if !haskey(data, "id") throw(ArgumentError("No ID provided for current pathway")) end
	if !haskey(data, "name") throw(ArgumentError("No name provided for current pathway")) end
	if !haskey(data, "genes") || !haskey(data, "type_id") throw(ArgumentError("You need to specify gene type id (entrez_id, hgnc_id, ensembl_id) and a list of genes for the pathway.")) end

	id_type = data["type_id"]
	genes = Vector{Gene}()
	for g in data["genes"]
		try
			gene = get_gene(experiment, g, id_type)
			push!(genes, gene)
		catch KeyError
			warn("no data for gene_id $g")
		end
	end
	info(typeof(genes))
	add_pathway(experiment, genes, name = string(data["name"]), id = string(data["id"]))
end

function add_gene_to_pathway(pw::Pathway{K}, gene::K) where {K<:KeyType}
	push!(pw.genes, gene)
	pw
end
######### PROTEINS ###########
function add_protein(experiment::Experiment, hgnc_id::String; fully_validated::Bool=false)
	p_obj = Protein(hgnc_id, fully_validated)
	experiment.proteins[hgnc_id] = p_obj
	p_obj
end

get_protein(experiment::Experiment, hgnc_id::String; fully_validated::Bool=false) = get!(experiment.proteins, hgnc_id, Protein(hgnc_id, fully_validated))

function add_protein(experiment::Experiment, data::Dict{String, Any})
	if !haskey(data, "hgnc_id")
		throw(ArgumentError("You need to provide a protein HGNC_ID."))
	end
	protein_id = data["hgnc_id"]
	ab_validated = get(entry, "antibody_validated", false)
	get_protein(experiment, hgnc_id, fully_validated=ab_validated)
end

######### VIEWS ###########
add_view!(experiment::Experiment, view::Type{<:ViewType}) = union!(experiment.views, [view])

######### DATA VIEWS ###########

# get a data view for a certain data type from a cell line or create if not present
get_dataview(cl::CellLine, data_type::Type{<:ViewType}) = cl.views[data_type]
get_dataview!(cl::CellLine, data_type::Type{<:ViewType}, dv::DataView{<:KeyType, <:ViewType}) = get!(cl.views, data_type, dv)

######### MEASUREMENTS ###########
# add a measurement, i.e. a gene/protein data pair to a data view
function add_measurement(d::DataView{<:KeyType,<:ViewType}, subject::KeyType, measurement::ViewType)
	if subject ∉ d.used_keys
		d.measurements[subject] = measurement
		union!(d.used_keys, [subject])
		return d
	end
	warn("measurement already exists for $subject, ignoring...")
	warn("if you want to force the new value, use add_measurement! instead")
	d
end
	
function add_measurement!(d::DataView{<:KeyType, <:VectorViewType}, subject::KeyType, measurement::VectorViewType)
	union!(d.used_keys, [subject])
	insert!(d.measurements, subject, measurement)
end

function add_measurement!(d::DataView{<:KeyType,<:SingleViewType}, subject::KeyType, measurement::SingleViewType)
	union!(d.used_keys, [subject])
	d.measurements[subject] = measurement
	d
end

function get_measurement!(d::DataView{<:KeyType, <:SingleViewType}, subject::KeyType, measurement::SingleViewType)
	union!(d.used_keys, [subject])
	get!(d.measurements, subject, measurement)
end

get_measurement(d::DataView{T,<:SingleViewType}, subject::T) where {T<:KeyType} = d.measurements[subject]

get_measurement_value(d::DataView{T, <:ViewType}, subject::T) where {T<:KeyType} = get_measurement_value.(d.measurements[subject])

get_measurement_value(ds::Vector{<:ViewType}) = get_measurement_value.(ds)

function get_measurement_value(d::ExomeSeq)
	(d.reference_mismatch_sum + d.variant_mismatch_sum)
end

get_measurement_value(d::GeneExpression) = d.expression_value
get_measurement_value(d::Methylation) = d.illumina_beta_value
get_measurement_value(d::RNASeq) = d.expression_value
get_measurement_value(d::RPPA) = d.protein_abundance
get_measurement_value(d::CNV) = d.gene_level_cnv


######### NORMALIZATION ###########

function normalize_data_views(experiment::Experiment)
	#for each of the views
	for v in experiment.views
		# Boolean vals need not be normalized
		if v == RNASeqCall continue end
		# reset the normalization statistics for this view
		experiment.statistics[v] = OrderedDict{KeyType, Tuple{Float64, Float64}}()
		# collect the keys in this view over all cell lines that it
		all_keys_in_view = Vector{KeyType}()
		for cl in values(experiment.cell_lines)
			# if the cell line has data for this view, collect the keys
			if haskey(cl.views, v)
				all_keys_in_view = union(all_keys_in_view, cl.views[v].used_keys)
			end
		end
		#iterate over all keys (i.e. genes/proteins/gene-protein change) and normalize across cell lines
		for key_id in 1:length(all_keys_in_view)
			key = all_keys_in_view[key_id]
			#collect values
			all_values = Float64[]
			for cl in values(experiment.cell_lines)
				if haskey(cl.views,v)
					if key ∈ cl.views[v].used_keys
						push!(all_values, get_measurement_value(cl.views[v].measurements[key])...)
					end
				end
			end
			mean_val = mean(all_values)
			std_val = stdm(all_values, mean_val)
			if std_val == 0.
				# warn("standard deviation is zero in view $v for key $key")
				std_val += 1e-7
			end
			#iterate over values again and normalize 
			for cl in values(experiment.cell_lines)
				if haskey(cl.views,v)
					if key ∈ cl.views[v].used_keys
						set_normalized_value(cl.views[v].measurements[key], (get_measurement_value(cl.views[v].measurements[key]) .- mean_val)./std_val)
					end
				end
			end
			# store statistics
			experiment.statistics[v][key] = (mean_val, std_val)
		end
	end
	experiment.is_normalized = true
end


function set_normalized_value(d::ViewType, value::Float64)
	d.normalized_value = value
end

set_normalized_value(ds::Vector{<:VectorViewType}, values::Vector{Float64}) = foreach(pair -> set_normalized_value(pair[1], pair[2]), zip(ds, values))

get_normalized_value(d::DataView{T, <:SingleViewType}, subject::T) where {T <: KeyType} = d.measurements[subject].normalized_value
get_normalized_value(d::ViewType) = d.normalized_value
get_normalized_value(ds::Vector{<:ViewType}) = get_normalized_value.(ds)
get_normalized_value(d::DataView{<:KeyType, <:VectorViewType}, subject::KeyType) = map(e -> e.normalized_value, d.measurements[subject])


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

function add_test_outcome!(e::Experiment, d::Drug, o::Outcome)
	e.test_results[d] = o
	e
end
