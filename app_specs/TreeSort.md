# Application specification: TreeSort

This is the application specification for the service with identifier TreeSort.

The backend script implementing the application is [App-TreeSort.pl](../service-scripts/App-TreeSort.pl).

The raw JSON file for this specification is [TreeSort.json](TreeSort.json).

This service performs the following task: TreeSort infers both recent and ancestral reassortment events along the branches of a phylogenetic tree of a fixed genomic segment.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------  |
| clades_path | Clades path | string | ??? | |
| deviation | Deviation | float? | | 2 |
| equal_rates | Equal rates | boolean | | |
| input_fasta_data | Input FASTA data | string | | |
| input_fasta_existing_dataset | Existing dataset | string | | |
| input_fasta_file_id | Input FASTA file ID | wsid | | |
| input_fasta_group_id | Input genome group | wsid | | |
| input_source | Input source | fasta_data, fasta_existing_dataset, fasta_file_id, fasta_group_id | :heavy_check_mark: | fasta_file_id |
| match_regex | A regular expression for matching | string (regex) |  |  |
| match_type  | Match on the "EPI_ISL_XXX" field, a regular expression, the strain, or the default | default, epi, regex, strain |  | default |
| method | Method | local, mincut | :heavy_check_mark: | local |
| no_collapse | No collapse | boolean | | |
| output_file | File basename | string | :heavy_check_mark: | |
| output_path | Output folder | string | :heavy_check_mark: | |
| p_value | P-value | float | | 0.001 |
| ref_segment | Reference Segment | string | :heavy_check_mark: | HA |
| ref_tree_inference | Reference tree inference method | FastTree, IQTree | | IQTree |
| segments | Segments | string | :heavy_check_mark: | All segments |






