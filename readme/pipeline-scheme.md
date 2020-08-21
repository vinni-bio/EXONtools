## EXONtools pipeline

```text
EXONtools.py
└── LAB EXPERIMENT 1: De novo transcriptome sequencing
        └── Stage A. Processing raw reads
               ├── Step A1. Demultiplexing (demultiplex_reads)
               ├── Step A2. Formatting (format_reads)
               ├── Step A3. Error correction (correct_reads)
               ├── Step A4. Duplication removal (deduplicate_reads)
               ├── Step A5. Trimming and filtering (clean_reads)
               ├── Step A6. Filtering low complexity (filter_reads)
               ├── Step A7. Merging paired reads (merge_reads)
               ├── Step A8. Contamination removal (decontaminate_reads)
               └── Step A9. Formatting (format_reads)

        └── Stage B. Pseudoreference construction
               ├── Step B1. Assembling reads (assemble_reads)
               ├── Step B2. Generating consensus assembly (consensus_assembly)
               ├── Step B3. Mapping reads to assembly (map_reads)
               ├── Step B4. Base calling (call_bases)
               ├── Step B5. Assembly annotation (annotate_contigs)
               ├── Step B6. Evaluate mapping quality (evaluate_mapping)
               └── Step B7. Evaluate assembly quality (evaluate_assembly)

         └── Stage C. Hybridization bait design
               ├── Step C1. Predicting exon boundaries (search_exons)
               ├── Step C2. Mapping exons to pseudoreference (map_exons)
               └── Step C3. Designing hybridization baits (design_baits)
           
└── LAB EXPERIMENT 2: Hybridization bait synthesis and exon capture

└── LAB EXPERIMENT 3: Exon capture sequencing
        └── Stage D. Processing raw reads
               ├── Step D1. Demultiplexing (demultiplex_reads)
               ├── Step D2. Formatting (format_reads)
               ├── Step D3. Error correction (correct_reads)
               ├── Step D4. Duplication removal (deduplicate_reads)
               ├── Step D5. Trimming and filtering (clean_reads)
               ├── Step D6. Filtering low complexity (filter_reads)
               ├── Step D7. Merging paired reads (merge_reads)
               ├── Step D8. Contamination removal (decontaminate_reads)
               └── Step D9. Formatting (format_reads)
        └── Stage E. SNP dataset construction
               ├── Step E1. Assembling reads (assemble_reads)
               ├── Step E2. Generating consensus assembly (consensus_assembly)
               ├── Step E3. Mapping reads to assembly (map_reads)
               ├── Step E4. Base calling (call_bases)
               ├── Step E5. Assembly annotation (annotate_contigs)
               ├── Step E6. Clustering exons (stack_exons)
               ├── Step E7. Aligning exons (align_stacks)
               ├── Step E8. Trimming and filtering exons (clean_stacks)
               ├── Step E9. Predicting exon boundaries (split_stacks)
               ├── Step E10. SNP calling from exon regions (call_snps)  
               ├── Step E11. SNP calling from intron regions (call_snps)
               └── Step E12. SNP filtering (call_snps)               
```
[GO BACK](https://github.com/vinni-bio/EXONtools#table-of-contents)