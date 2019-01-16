# DrugEfficacyPredictor

## Uploading Data
We provide two initial data sets that can be used:

    - "dream_challenge": the data set used for the drug efficacy prediction DREAM challenge on which the algorithm is based and which is reported in [1]
    - "iorio": data on a much larger panel of cell lines and drugs as reported in [2]

If you want to use your own data, you can do so by uploading to the server:
1) create an experiment: ```curl -X POST http://{host}:{port}/experiments -d '{"experiment_id" : "your_experiment"}'```
    - then data part is optional; if omitted, a random experiment_id is created, same is tue when the id is already in use
2) upload your data: ```curl -X POST http://{host}:{port}/experiments/{experiment_id}/data/{data_type} -d @{path to dat in json}```
    - data_type can be one of the following: genes, proteins, cell_lines, drugs, pathways, outcomes
    - format for the JSON files is described below
3) proceed with steps 3-6 as in the example below

## Data formats

### creating experiments
```
{
    "experiment_id" : "your_experiment_id"
}
```

### Genes
```
{
    "genes" :
    [
        {
            "gene_id" : "some_gene_id",
            "type_id" : "type_of_id (can be hgnc_id, entrez_id or ensembl_id)"
            "is_cancer_gene" : true|false
        }, ...
    ]
}
```

### Proteins
```
{
    "proteins" :
    [
        {
            "hgnc_id" : "the_hgnc_id (currently, only hgnc ids are supported for proteins)",
            "antibody_validated" : true|false
        }, ...
    ]
}
```

### Cell lines

```
{
    "cell_lines" :
    [
        {
            "id" : "cell line id".
            "cancer_type" : "type of cancer for this cell line",
            "in_test_set" : true|false (whether this cell line is used for training or testing),
            "views" :
            {
                "view_type (can be GeneExpression, Methylation, RNASeq, RNASeqCall, ExomeSeq, RPPA or CNV)" :
                [
                    # measurement as described below
                ], ... (further views)
            }
        }
    ]
}
```

#### GeneExpression data
```
{
    "gene_id" : "the_gene_id",
    "type_id" : "gene_id_type",
    "is_cancer_gene" : true|false (optional here, uploading cell lines creates data on the fly, so can be useful if gene does not exist in the system yet),

    "expression_value" : float number (the actual expression measurement)
}
```

#### Methylation data
```
{
    # gene stuff

    # mandatory:
    "beta_value" : float number (the illumination beta value)
    # optional
    "cgct1" : integer,
    "cct1" : integer,
    "methylated_threshold" : float number (threshold above which beta value indicates methylation),
    "illumina_id" : string
}
```

#### RNASeq data
```
{
    # gene stuff
    "expression_value" : float number
}
```

#### RNASeqCall data
```
{
    # gene stuff
    "expression_status" : float number
}
```

#### ExomeSeq data
```
{
    # gene stuff

    # mandatory
    "protein_change" : string
    #optional
    "reference_mismatch_avg" : float,
    "reference_mismatch_sum" : float,
    "reference_dist3effective_avg" : float,
    "variant_mismatch_avg" : float,
    "variant_mismatch_sum" : float,
    "variant_dist3effective_avg" : float,
    "num_cosmic" : integer,
    "is_cancer_gene" : bool,
    "variant_effect" : string,
    "nucleotid_change" : string,
    "variant_confidence" : float,
    "norm_zygosity" : string,
    "norm_reference_count" : float,
    "norm_variant_count" : float,
    "tumor_zygosity" : string,
    "tumor_reference_count" : float,
    "tumor_variant_count" : float,
    "details" : string
}
```

#### RPPA data
```
{
    "hgnc_id" : string, # the protein's hgnc id
    "antibody_validated" : bool, #optional here, same as for gene stuff above
    "protein_abundance" : float
}
```

#### CNV data
```
{
    # gene stuff

    "gene_level_cnv" : float (the cnv level aggregated across a gene as identified by #gene stuff)
}
```

### Drugs
```
{
    "drugs" :
    [
        {
            # mandatory
            "id" : string,
            # optional
            "affected_genes" : [#gene stuff, ...],
            "chemical_structure" : string
        }, ...
    ]
}
```

### Pathways
```
{
    "pathways" :
    [
        {
            "id" : string,
            "name" : string,
            "genes" : [# gene stuff, ...]
        }, ...
    ]

}
```

### Outcomes
```
{
    "outcome_type" : string,
    "outcomes" :
    {
         "{drug_name}" :
        {
            "{cell_line_id}" : float,
            ...
        },
        ...
    }
}
```

## Example
1) start the server: ```DrugEfficacyPredictor.start_server(port)``` will start a server on localhost at the specified port
2) load example data or upload your own data
    - example data can be loaded via command line with ```load_iorio_data()``` (experiment id "iorio") or ```load_dream_challenge_data()``` (experiment id "dream_challenge")
3) start the training process: ```curl -X POST http://{host}:{port}/experiments/{experiment id}/train -d @{path to config file}```, an example config file is provided in the ```data/configs``` directory
4) query the server for training progress (optionally)
5) download the training result directory as zip file (once the training is finished): ```curl -X GET http://{host}:{port}/experiments/{experiment id}/results```
6) predict response data for arbitrary additional cell lines (if applicable):

## References
[1] Costello, J. C., Heiser, L. M., Georgii, E., & Gönen, M. (2014). A community effort to assess and improve drug sensitivity prediction algorithms. Nature.
[2] Iorio, F., Knijnenburg, T. A., Vis, D. J., Bignell, G. R., Menden, M. P., Schubert, M., et al. (2016). A Landscape of Pharmacogenomic Interactions in Cancer. Cell, 166(3), 740–754. http://doi.org/10.1016/j.cell.2016.06.017

## TODOs
- full grid search, no we are setting all gamma distributions with equal hyper params, same for normals
- cross validation
    + check wpc index on test set alone (in Costello et al. always checked on combined data set, training and test) (ok, no difference)
    + compute correlation of accuracy in cv runs to actual outcomes
    + compute wpc internally - call perl scripts for each fold run and on test sets (ok, dream challenge only)
- pipeline
    + create cross validation methods (ok)
    + combine result data and create different splits of the data (ok)
    + write splits to filesystem (ok)
    + run external wpc scripts on newly created data splits (ok, dream challenge only)
    + best results so far, when \[.1 \le E[\nu] \le 10\] for \[\nu \sim \mathcal G(\cdot,\cdot)\] a Gamma-distributed RV
