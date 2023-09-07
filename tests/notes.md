# Coverage

## Align WITH fragment building as starting point

## test_controller_coverage_paired_end_fragments_align_fragments.py
coverage parameters: default, count fragments, normalize by fragments 
(--no_fragments=False, --no_norm_by_fragments=False)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_coverage_paired_end_input_files
#### Test:
coverage: reademption coverage --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_coverage_paired_expected

### test_controller_coverage_paired_end_no_fragments_align_fragments.py
coverage parameters: count reads, normalize by reads
(--no_fragments=True, --no_norm_by_fragments=False)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_coverage_paired_end_input_files
#### Test:
coverage_no_fragments: reademption coverage --no_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_coverage_paired_no_fragments_expected

### test_controller_coverage_paired_end_no_norm_fragments_align_fragments.py
coverage parameters: count fragments, normalize by reads
(--no_fragments=False, --no_norm_by_fragments=True)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_coverage_paired_end_input_files
#### Test
coverage_no_norm_fragments: reademption coverage --no_norm_by_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected_folder: reademption_analysis_dual_with_coverage_paired_fragments_but_no_norm_fragments_expected


## Align WITHOUT fragment building as starting point

### test_controller_coverage_paired_end_fragments_align_no_fragments.py
coverage parameters: default, count fragments, normalize by fragments
(--no_fragments=False, --no_norm_by_fragments=False)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_coverage_paired_end_no_fragments_input_files
#### Test:
coverage: reademption coverage --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_coverage_paired_expected

### test_controller_coverage_paired_end_no_fragments_align_no_fragments.py
coverage parameters: count reads, normalize by reads
(--no_fragments=True, --no_norm_by_fragments=False)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_coverage_paired_end_no_fragments_input_files
#### Test:
coverage_no_fragments: reademption coverage --no_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_coverage_paired_no_fragments_expected

### test_controller_coverage_paired_end_no_norm_fragments_align_no_fragments.py
coverage parameters: count fragments, normalize by reads
(--no_fragments=False, --no_norm_by_fragments=True)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_coverage_paired_end_no_fragments_input_files
#### Test
coverage_no_norm_fragments: reademption coverage --no_norm_by_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected_folder: reademption_analysis_dual_with_coverage_paired_fragments_but_no_norm_fragments_expected




# Gene Quantification

## Align WITH fragment building as starting point

## test_controller_gene_quanti_paired_end_fragments_align_fragments.py
gene_quanti parameters: default, count fragments, normalize by fragments
(--no_fragments=False, --no_norm_by_fragments=False)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_gene_quanti_paired_end_input_files
#### Test:
gene_quanti: reademption gene_quanti --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_gene_quanti_paired_expected


### test_controller_gene_quanti_paired_end_no_fragments_align_fragments.py
gene_quanti parameters: count reads, normalize by reads
(--no_fragments=True, --no_norm_by_fragments=False)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_gene_quanti_paired_end_input_files
#### Test:
gene_quanti_no_fragments: reademption gene_quanti --no_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_gene_quanti_paired_no_fragments_expected


### test_controller_gene_quanti_paired_end_no_norm_fragments_align_fragments.py
gene_quanti parameters: count fragments, normalize by reads
(--no_fragments=False, --no_norm_by_fragments=True)
#### Setup:
align_split: reademption align --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_with_gene_quanti_paired_end_input_files
#### Test
gene_quanti_no_norm_fragments: reademption gene_quanti --no_norm_by_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected_folder: reademption_analysis_dual_with_gene_quanti_paired_fragments_but_no_norm_fragments_expected


## Align WITHOUT fragment building as starting point

### test_controller_gene_quanti_paired_end_fragments_align_no_fragments.py
gene_quanti parameters: default, count fragments, normalize by fragments
(--no_fragments=False, --no_norm_by_fragments=False)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_gene_quanti_paired_end_no_fragments_input_files
#### Test:
gene_quanti: reademption gene_quanti --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_gene_quanti_paired_expected

### test_controller_gene_quanti_paired_end_no_fragments_align_no_fragments.py
gene_quanti parameters: count reads, normalize by reads
(--no_fragments=True, --no_norm_by_fragments=False)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_gene_quanti_paired_end_no_fragments_input_files
#### Test:
gene_quanti_no_fragments: reademption gene_quanti --no_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected folder: reademption_analysis_dual_with_gene_quanti_paired_no_fragments_expected

### test_controller_gene_quanti_paired_end_no_norm_fragments_align_no_fragments.py
gene_quanti parameters: count fragments, normalize by reads
(--no_fragments=False, --no_norm_by_fragments=True)
#### Setup:
align_split_no_fragments: reademption align --no_fragment_building --split -p 1 -q -a 100 --paired_end --project_path ${PROJECT_NAME}
input folder: reademption_analysis_dual_paired_end_with_gene_quanti_paired_end_no_fragments_input_files
#### Test
gene_quanti_no_norm_fragments: reademption gene_quanti --no_norm_by_fragments --paired_end -p 4  --project_path ${PROJECT_NAME}
expected_folder: reademption_analysis_dual_with_gene_quanti_paired_fragments_but_no_norm_fragments_expected
