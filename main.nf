// enabling nextflow DSL v2
nextflow.enable.dsl=2

process GET_METADATA {

    publishDir "${params.resultsDir}/metadata/", mode: 'copy', overwrite: true

    output:
        path("metadata.csv")

    script:
    """
        cubseq-get-ENA-metadata.R \\
            --taxon-id ${params.taxonId} \\
            --limit-search ${params.limitSearch} \\
            --remove-run '${params.removeRun}' \\
            --max-sra-bytes ${params.max_sra_bytes} \\
            --date-min ${params.dateMin} \\
            --date-max ${params.dateMax} \\
            metadata.csv
    """

    stub:
    """
        echo "tax_id,study_accession,experiment_accession,sample_accession,run_accession,library_source,instrument_platform,library_strategy,library_layout,read_count,base_count,fastq1,fastq2" >> metadata.csv
        echo "562,PRJDB13439,DRX362496,SAMD00466964,DRR376595,TRANSCRIPTOMIC,ILLUMINA,RNA-Seq,PAIRED,74553652,22366095600,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376595/DRR376595_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376595/DRR376595_2.fastq.gz" >> metadata.csv
        echo "562,PRJDB13439,DRX362497,SAMD00466965,DRR376596,TRANSCRIPTOMIC,ILLUMINA,RNA-Seq,PAIRED,102884857,30865457100,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376596/DRR376596_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376596/DRR376596_2.fastq.gz" >> metadata.csv
        echo "562,PRJDB13439,DRX362498,SAMD00466966,DRR376597,TRANSCRIPTOMIC,ILLUMINA,RNA-Seq,PAIRED,24784574,7435372200,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376597/DRR376597_1.fastq.gz,ftp.sra.ebi.ac.uk/vol1/fastq/DRR376/DRR376597/DRR376597_2.fastq.gz" >> metadata.csv
    """

}

process GET_READS {

    tag "${run_acc}"
    
    input:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), val(fastq1), val(fastq2)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}-1.fastq.gz"), path("${run_acc}-2.fastq.gz")

    script:
    """
        curl -L -Z ${fastq1} -o ${run_acc}-1.fastq.gz && \\
        curl -L -Z ${fastq2} -o ${run_acc}-2.fastq.gz
    """

    stub:
    """
        touch ${run_acc}-1.fastq.gz
        touch ${run_acc}-2.fastq.gz
    """

}

process PREPROCESS_READS {

    tag "${run_acc}"

    input:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(read1), path(read2)
    
    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${read1.simpleName}.trimmed.fastq.gz"), path("${read2.simpleName}.trimmed.fastq.gz"), emit: trimmed_fastq
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.fastp.json"), path("${run_acc}.fastp.html"), emit: qc

    script:
    """
        fastp -w ${task.cpus} \\
            --in1 ${read1} \\
            --in2 ${read2} \\
            --out1 ${read1.simpleName}.trimmed.fastq.gz \\
            --out2 ${read2.simpleName}.trimmed.fastq.gz \\
            --json ${run_acc}.fastp.json \\
            --html ${run_acc}.fastp.html
    """

    stub:
    """
        touch ${read1.simpleName}.trimmed.fastq.gz
        touch ${read2.simpleName}.trimmed.fastq.gz
        touch ${run_acc}.fastp.json
        touch ${run_acc}.fastp.html
    """

}

process GENERATE_GENOME_INDEX {

    tag "genome_index"

    input:
        path(genome_fasta)
        path(annotation_gtf)

    output:
        path("genome-index")

    script:
    """
        STAR \\
            --runThreadN ${task.cpus} \\
            --runMode genomeGenerate \\
            --genomeDir genome-index \\
            --genomeFastaFiles ${genome_fasta} \\
            --sjdbGTFfile ${annotation_gtf} \\
            --sjdbOverhang ${params.star.sjdbOverhang} \\
            --genomeSAindexNbases ${params.star.genomeSAindexNbases}
    """

    stub:
    """
        mkdir genome-index
    """

}

process ALIGN_READS {

    publishDir "${params.resultsDir}/bams/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(genome_index)
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(trimmed_read1), path(trimmed_read2)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.Aligned.sortedByCoord.out.bam")

    script:
    """
        STAR \\
            --runThreadN ${task.cpus} \\
            --genomeDir ${genome_index} \\
            --readFilesIn ${trimmed_read1} ${trimmed_read2} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${run_acc}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --alignIntronMax ${params.star.alignIntronMax} \\
            --limitBAMsortRAM ${params.star.limitBAMsortRAM} \\
            --outBAMsortingBinsN ${params.star.outBAMsortingBinsN}
    """

    stub: 
    """
        touch ${run_acc}.Aligned.sortedByCoord.out.bam
    """
    
}

process QUANTIFY_READS {

    publishDir "${params.resultsDir}/featureCounts/${run_acc}/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(annotation_gtf)
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(bam_file)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.featureCounts.txt"), path("${run_acc}.featureCounts.txt.summary")

    script:
    """
        featureCounts \\
            -T ${task.cpus} \\
            -p \\
            -t ${params.featureCounts.type.feature} \\
            -g ${params.featureCounts.type.attribute} \\
            -a ${annotation_gtf} \\
            -o ${run_acc}.featureCounts.txt \\
            ${bam_file}
    """

    stub:
    """
        touch ${run_acc}.featureCounts.txt
        touch ${run_acc}.featureCounts.txt.summary
    """

}

process FREEBAYES_CALL_VARIANTS {

    publishDir "${params.resultsDir}/freebayes-vcf/", mode: 'copy', overwrite: true
    
    tag "${run_acc}"

    input:
        path(genome_fasta)
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(bam_file)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.var.vcf")

    script:
    """
        freebayes \\
            --fasta-reference ${genome_fasta} \\
            --ploidy ${params.freebayes.ploidy} \\
            --bam ${bam_file} \\
            --vcf ${run_acc}.var.vcf
    """

    stub:
    """
        touch ${run_acc}.var.vcf
    """

}

process BCFTOOLS_FILTER_VCF {

    publishDir "${params.resultsDir}/vcf/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(vcf)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.filter.norm.var.vcf.gz")

    script:
    """
        # index vcf with tabix
        cat ${vcf} | bgzip > ${run_acc}.var.vcf.gz
        tabix --preset vcf ${run_acc}.var.vcf.gz

        # normalise to decompose MNVs into consecutive SNVs
        bcftools norm \\
            --threads ${task.cpus} \\
            --atomize ${run_acc}.var.vcf.gz \\
            --output-type z \\
            --output ${run_acc}.norm.var.vcf.gz

        # filter vcf by specified params (e.g. quality score, include only SNPs)
        bcftools filter \\
            --threads ${task.cpus} \\
            --include '${params.bcftools.filter_vcf.args}' \\
            --output-type z \\
            --output ${run_acc}.filter.norm.var.vcf.gz \\
            ${run_acc}.norm.var.vcf.gz
    """

    stub:
    """
        touch ${run_acc}.filter.norm.var.vcf.gz
    """

}

process BCFTOOLS_CREATE_CONSENSUS {

    publishDir "${params.resultsDir}/transcriptome-consensus/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(genome_fasta)
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(filtered_vcf)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.consensus.fa")

    script:
    """
        # index vcf with tabix
        tabix --preset vcf ${filtered_vcf}

        # create consensus sequence
        bcftools consensus \\
            --fasta-ref ${genome_fasta} \\
            --output ${run_acc}.consensus.fa \\
            ${filtered_vcf}
    """

    stub:
    """
        touch ${run_acc}.consensus.fa
    """

}

process GFFREAD_GET_WT_TRANSCRIPTOME {

    publishDir "${params.resultsDir}/wt-transcriptome/", mode: 'copy', overwrite: true

    input:
        path(genome_fasta)
        path(annotation_gtf)

    output:
        path("wt-transcriptome.fa"), emit: transcriptome
        path("wt-transcriptome.gtf"), emit: annotation
    // TODO: -w param output file not returning spliced exons?
    script:
    """
        gffread \\
            -g ${genome_fasta} \\
            -o wt-transcriptome.gtf \\
            -w wt-transcriptome.fa \\
            --gtf ${annotation_gtf} \\
            -L
    """

    stub:
    """
        touch wt-transcriptome.fa
        touch wt-transcriptome.gtf
    """

}

process GFFREAD_GET_MUT_TRANSCRIPTOME {

    publishDir "${params.resultsDir}/mut-transcriptome/${run_acc}/", mode: 'copy', overwrite: true

    input:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(consensus_fasta)
        path(annotation_gtf)

    output:
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.mut-transcriptome.fa"), emit: transcriptome
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path("${run_acc}.mut-transcriptome.gtf"), emit: annotation
    // TODO: -w param output file not returning spliced exons?
    script:
    """
        gffread \\
            -g ${consensus_fasta} \\
            -o ${run_acc}.mut-transcriptome.gtf \\
            -w ${run_acc}.mut-transcriptome.fa \\
            --gtf ${annotation_gtf} \\
            -L
    """

    stub:
    """
        touch ${run_acc}.mut-transcriptome.fa
        touch ${run_acc}.mut-transcriptome.gtf
    """

}

process SALMON_INDEX {

    input: 
        path(transcriptome)
        path(genome)

    output: 
        path("transcriptome-index")

    script:
    """
        # extract names of genome targets
        grep '^>' < ${genome} | cut -d " " -f 1 > decoys.txt
        sed -i.bak -e 's/>//g' decoys.txt

        # concatenate transcriptome and genome reference file for index
        cat ${transcriptome} ${genome} > gentrome.fa.gz

        salmon index \\
            --threads ${task.cpus} \\
            --transcripts gentrome.fa.gz \\
            --decoys decoys.txt \\
            --index transcriptome-index \\
            --keepDuplicates \\
            ${params.salmon.index.args}
    """

    stub:
    """
        mkdir transcriptome-index
    """

}

process SALMON_QUANT {

    publishDir "${params.resultsDir}/salmon-quant/", mode: 'copy', overwrite: true

    tag "${run_acc != null ? run_acc : 'no_tag' }"

    input: 
        path(transcriptome_index)
        tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(read1), path(read2)

    output:
        path("${run_acc}")

    script:
    """
        salmon quant \\
            --threads ${task.cpus} \\
            --libType ${params.salmon.quant.libtype} \\
            --index ${transcriptome_index} \\
            --mates1 ${read1} \\
            --mates2 ${read2} \\
            --output ${run_acc} \\
            ${params.salmon.quant.args}
    """

    stub:
    """
        mkdir ${run_acc}
    """
    
}

// CUB PRE-PROCESSING
process SUMMARIZE_TO_GENE {

    publishDir "${params.resultsDir}/dataset/", mode: 'copy', overwrite: true

    input:
        path(metadata_csv)
        path(annotation_gtf)
        path(quant_dir, stageAs: 'quant_dir/quant*.sf')

    output:
        path("${params.summarize_to_gene.output}")

    script:
    """
        cubseq-summarize-to-gene.R \\
            ${metadata_csv} \\
            --gtf ${annotation_gtf} \\
            --quant-dir ${quant_dir} \\
            --counts-from-abundance ${params.summarize_to_gene.counts_from_abundance} \\
            --output ${params.summarize_to_gene.output}
    """

    stub:
    """
        mkdir txi-summarized-experiment.rds 
    """

}

process RANK_GENES {

    publishDir "${params.resultsDir}/gene-rank-analysis/", mode: 'copy', overwrite: true

    input:
        path(txi)
        path(metadata)
        path(gtf)
        path(gff)

    output:
        path("rank-tpm.csv"), emit: rank_tpm
        path("chisq-results.csv"), emit: chisq_results
        path("heg-scores.csv"), emit: heg_scores
        path("leg-scores.csv"), emit: leg_scores
        path("proteins_tpm.rds"), emit: proteins_tpm

    script:
    """
        cubseq-fisher-transform.R \\
            ${txi} \\
            --metadata ${metadata} \\
            --gtf ${gtf} \\
            --gff ${gff} \\
            --fdr-threshold ${params.rank_genes.fdr}
    """

    stub:
    """
        touch rank-tpm.csv
        touch chisq-results.csv
        touch heg-scores.csv
        touch leg-scores.csv
        touch proteins_tpm.rds
    """

}

process GET_HEG_LEG_FASTA {

    publishDir "${params.resultsDir}/heg-mut-transcriptome/", pattern: '*.heg-mut-transcriptome.fa', mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/leg-mut-transcriptome/", pattern: '*.leg-mut-transcriptome.fa', mode: 'copy', overwrite: true

    input:
        path(metadata)
        path(gtf)
        path(mut_fasta)
        path(heg_geneID)
        path(leg_geneID)

    output:
        path("*.heg-mut-transcriptome.fa"), emit: heg
        path("*.leg-mut-transcriptome.fa"), emit: leg
    
    script:
    """
        cubseq-get-heg-leg-fasta.R \\
            ${metadata} \\
            --gtf ${gtf} \\
            --mut-transcriptome-dir ${mut_fasta} \\
            --heg-geneID ${heg_geneID} \\
            --leg-geneID ${leg_geneID}
    """

    stub:
    """
        mkdir heg-mut-transcriptome
        mkdir leg-mut-transcriptome
    """

}

process GET_PROTEIN_FASTA {

    publishDir "${params.resultsDir}/protein-mut-transcriptome/", pattern: '*.protein-mut-transcriptome.fa', mode: 'copy', overwrite: true

    input:
        path(metadata)
        path(gtf)
        path(mut_fasta)

    output:
        path("*.protein-mut-transcriptome.fa"), emit: protein_fasta
    
    script:
    """
        cubseq-get-protein-fasta.R \\
            ${metadata} \\
            --gtf ${gtf} \\
            --mut-transcriptome-dir ${mut_fasta}
    """

    stub:
    """
        mkdir protein-mut-transcriptome
    """

}

process GET_CU_FREQUENCIES {

    publishDir "${params.resultsDir}/cu-data/", mode: 'copy', overwrite: true

    input:
        path(heg_fasta_dir)
        path(leg_fasta_dir)
        path(protein_fasta_dir)

    output:
        path("heg-codon-counts.csv"), emit: heg
        path("leg-codon-counts.csv"), emit: leg
        path("protein-codon-counts.csv"), emit: protein
        path("kazusa-*-codon-counts.csv"), emit: kazusa
        path("cocoputs-*-codon-counts.csv"), emit: cocoputs
        path("*-rf.csv"), emit: rf
    
    script:
    """
        cubseq-get-cu.R \\
            --heg-fasta-dir ${heg_fasta_dir} \\
            --leg-fasta-dir ${leg_fasta_dir} \\
            --protein-fasta-dir ${protein_fasta_dir} \\
            --tax-id ${params.get_cu_frequencies.tax_id}
    """

    stub:
    """
        touch heg-codon-counts.csv
        touch leg-codon-counts.csv
        touch protein-codon-counts.csv
        touch kazusa-*-codon-counts.csv
        touch cocoputs-*-codon-counts.csv
        touch *-rf.csv
    """

}

// process GET_KO_ID {

//     publishDir "${params.resultsDir}/ko-mut-transcriptome/${run_acc}/", mode: 'copy', overwrite: true

//     tag "${run_acc}"

//     input:
//         tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(mut_fasta)
    
//     output:
//         path("${run_acc}.ko-mut-transcriptome.fa")

//     script:
//     """
//         cubseq-get-KO.R \\
//             ${mut_fasta} \\
//             --org-id ${params.get_ko_id.org} \\
//             --output ${run_acc}.ko-mut-transcriptome.fa
//     """

//     stub:
//     """
//         touch ${run_acc}.ko-mut-transcriptome.fa
//     """
    
// }

// process CORDON_COUNT_CODONS {

//     publishDir "${params.resultsDir}/coRdon-analysis/${run_acc}/", mode: 'copy', overwrite: true

//     tag "${run_acc}"

//     input:
//         tuple val(study_acc), val(sample_acc), val(experiment_acc), val(run_acc), path(mut_transcripts_fasta)

//     output:
//         path("${run_acc}.codon_table.txt"), emit: codon_table
//         path("${run_acc}.cub_measure.txt"), emit: cub_measure
//         path("${run_acc}.cub_expressivity.txt"), emit: cub_expressivity

//     """
//         cubseq-coRdon.R \\
//             ${mut_transcripts_fasta} \\
//             --cub-statistic ${params.coRdon.cub_statistic} \\
//             --len-threshold ${params.coRdon.len_threshold} \\
//             ${run_acc}.codon_table.txt \\
//             ${run_acc}.cub_measure.txt \\
//             ${run_acc}.cub_expressivity.txt
//     """

//     stub:
//     """
//         touch ${run_acc}.codon_table.txt
//         touch ${run_acc}.cubTable.txt
//         touch ${run_acc}.cub_expressivity.txt
//     """

// }


process GET_CU_FIGURES {

    publishDir "${params.resultsDir}/figures/", mode: 'copy', overwrite: true

    input:
        path(cu_table)
        path(aa_properties_data)

    output:
        path("figure1.pdf")
        path("figure2.pdf")
        path("figure3.pdf")
        path("figure4.pdf")
    
    script:
    """
        cubseq-get-cu-figures.R \\
            --cu-table ${cu_table} \\
            --aa-property ${aa_properties_data}
    """

    stub:
    """
        touch figure1.pdf
        touch figure2.pdf
        touch figure3.pdf
        touch figure4.pdf
    """

}

process GET_PCA {

    publishDir "${params.resultsDir}/pca-figures/", pattern: "*.pdf", mode: 'copy', overwrite: true
    publishDir "${params.resultsDir}/pca-data/", pattern: "*.{csv,rds}", mode: 'copy', overwrite: true

    input:
        path(txi_file)
        path(metadata)

    output:
        path("pca-panel.pdf")
        path("pca-study.pdf")
        path("pca-pc-importance.pdf")
        path("pca-metadata.csv")
        path("pca-variance.csv")
        path("pca-summary.rds")
    
    script:
    """
        cubseq-get-pca.R \\
            ${txi_file} \\
            --metadata ${metadata}
    """

    stub:
    """
        touch pca-panel.pdf
        touch pca-study.pdf
        touch pca-pc-importance.pdf
        touch pca-metadata.csv
        touch pca-variance.csv
        touch pca-summary.rds
    """

}

process STRAIN_CU_HEATMAP {

    publishDir "${params.resultsDir}/strain-cu/", mode: 'copy', overwrite: true

    input:
        path(metadata)
        path(fasta_dir)

    output:
        path("strain-tax.rds")
        path("strain-rf.csv")

    script:
    """
        cubseq-strain-cu-heatmap.R \\
            ${metadata} \\
            --fasta-dir ${fasta_dir} \\
            --heatmap-width ${params.strain_cu_heatmap.width} \\
            --heatmap-height ${params.strain_cu_heatmap.height}
    """

    stub:
    """
        touch strain-tax.rds
        touch strain-rf.csv
    """

}

process STRAIN_CU_CAT_PLOT {

    publishDir "${params.resultsDir}/strain-cat/", mode: 'copy', overwrite: true

    tag "${strain}"

    input:
        path(metadata_csv)
        path(annotation_gtf)
        tuple val(strain), path(quant_files, stageAs: 'quant*.sf')

    output:
        tuple val(strain), path("${strain}-txi-summarized-experiment.rds")

    script:
    """
        cubseq-summarize-to-gene.R \\
            ${metadata_csv} \\
            --gtf ${annotation_gtf} \\
            --counts-from-abundance ${params.summarize_to_gene.counts_from_abundance} \\
            --output ${strain}-txi-summarized-experiment.rds
    """

    stub:
    """
        touch ${strain}-txi-summarized-experiment.rds
    """

}

process STRAIN_HEG_CAT_PLOT {

    publishDir "${params.resultsDir}/strain-cat/", mode: 'copy', overwrite: true

    input:
        path(txi)
        path(metadata_csv)
        path(annotation_gtf)

    output:
        path("cat-plot.pdf")

    script:
    """
        cubseq-chi-meta-analysis.R \\
            ${txi} \\
            --metadata ${metadata_csv} \\
            --gtf ${annotation_gtf} \\
            --fdr-threshold ${params.strain_heg_cat_plot.fdr} \\
            --num-strains ${params.strain_heg_cat_plot.num_strains}
    """

    stub:
    """
        touch cat-plot.pdf
    """

}

process GET_STRAIN_CAT_PLOT {

    publishDir "${params.resultsDir}/strain-cat-data/", mode: 'copy', overwrite: true

    input:
        path(txi)
        path(metadata_csv)

    output:
        path("heg_tax_concordance.rds")
        path("leg_tax_concordance.rds")
        path("cat-plot-HEG.pdf")
        path("cat-plot-LEG.pdf")

    script:
    """
        cubseq-strain-cat-plot.R \\
            ${txi} \\
            --metadata ${metadata_csv} \\
            --tax-names ${params.get_strain_cat_plot.tax_names} \\
            --heg-scores ${params.get_strain_cat_plot.heg_scores} \\
            --leg-scores ${params.get_strain_cat_plot.leg_scores} \\
            --num-strains ${params.get_strain_cat_plot.num_strains} \\
            --fdr-threshold ${params.get_strain_cat_plot.fdr}
    """

    stub:
    """
        touch heg_tax_concordance.rds
        touch leg_tax_concordance.rds
        touch cat-plot-HEG.pdf
        touch cat-plot-LEG.pdf
    """

}

// variant analysis
process ANNOTATE_VCF {

    publishDir "${params.resultsDir}/annotate-vcf/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(gtf)
        tuple val(run_acc), path(vcf)

    output:
        tuple val(run_acc), path("${run_acc}.annotated-vcf.gtf")

    script:
    """
        bedtools annotate \\
            -both \\
            -i ${gtf} \\
            -files ${vcf} > ${run_acc}.annotated-vcf.gtf
    """

    stub:
    """
        touch ${run_acc}.annotated-vcf.gtf
    """

}

process PLOT_TOP_MUTATED_GENES {

    publishDir "${params.resultsDir}/top-mutated-genes/", mode: 'copy', overwrite: true

    input:
        path(metadata)
        path(mut_gtf_dir)

    output:
        path("vcf-bedtools-annotate-counts.rds")
        path("gene-mutation-fraction.csv")
        path("sample-mutation-counts.csv")
        path("top_mutated_gene_count_barplot_20.pdf")
        path("top-mutated-genes-barplot_100.pdf")
        path("top-mutated-genes-barplot_50.pdf")

    script:
    """
        cubseq-get-top-mutated-genes.R \\
            ${metadata} \\
            --annotated-mut-dir ${mut_gtf_dir}
    """

    stub:
    """
        touch vcf-bedtools-annotate-counts.csv
        touch gene-mutation-fraction.csv
        touch sample-mutation-counts.csv
        touch top_mutated_gene_count_barplot_20.pdf
        touch top-mutated-genes-barplot_100.pdf
        touch top-mutated-genes-barplot_50.pdf
    """

}

process CALL_CONSEQUENCE {

    publishDir "${params.resultsDir}/vcf-csq-calls/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(fasta)
        path(gff)
        tuple val(run_acc), path(vcf)

    output:
        tuple val(run_acc), path("${run_acc}.csq.vcf.gz")

    script:
    """
        # perform consequence calling
        bcftools csq \\
            --threads ${task.cpus} \\
            -f ${fasta} \\
            -g ${gff} \\
            ${vcf} \\
            -Oz \\
            -o ${run_acc}.csq.vcf.gz
    """

    stub:
    """
        touch ${run_acc}.csq.vcf.gz
    """

}

process CALL_CONSEQUENCE_ON_MERGED {

    publishDir "${params.resultsDir}/merged-vcf-csq-calls/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        path(fasta)
        path(gff)
        path(vcf)

    output:
        path("csq.vcf.gz")

    script:
    """
        # perform consequence calling
        bcftools csq \\
            --threads ${task.cpus} \\
            -f ${fasta} \\
            -g ${gff} \\
            ${vcf} \\
            -Oz \\
            -o csq.vcf.gz
    """

    stub:
    """
        touch csq.vcf.gz
    """

}

process MERGE_VCF {

    publishDir "${params.resultsDir}/merged-vcf/", mode: 'copy', overwrite: true

    input:
        path(vcf)

    output:
        path("merged.vcf.gz")
    
    script:
    """
        # index vcf files
        parallel tabix -p vcf ::: ${vcf}

        # merge vcf files
        bcftools merge --threads ${task.cpus} --force-samples -Oz ${vcf} > merged.vcf.gz
    """

    stub:
    """
        touch merged.vcf.gz
    """

}

process GET_VCF_STATS {

    publishDir "${params.resultsDir}/vcf-stats/", mode: 'copy', overwrite: true

    input:
        path(merged_vcf)
        path(fasta)

    output:
        path("merged-vcf-stats.vchk")
    
    script:
    """
        # parse merged VCF to produce stats
        bcftools stats --threads ${task.cpus} --fasta-ref ${fasta} ${merged_vcf} > merged-vcf-stats.vchk
    """

    stub:
    """
        touch merged-vcf-stats.vchk
    """

}

process DROP_VCF_GENOTYPES {

    publishDir "${params.resultsDir}/vcf-drop-genotypes/", mode: 'copy', overwrite: true

    input:
        path(merged_vcf)

    output:
        path("merged-drop-genotype.vcf.gz")
    
    script:
    """
        # drop genotypes from merged VCF
        bcftools view \\
            --threads ${task.cpus} \\
            --drop-genotypes \\
            --output-type z \\
            --output merged-drop-genotype.vcf.gz \\
            ${merged_vcf}
    """

    stub:
    """
        touch merged-drop-genotype.vcf.gz
    """

}

process PLOT_VCF_STATS {

    publishDir "${params.resultsDir}/vcf-stats-plots/", mode: 'copy', overwrite: true

    input:
        path(vcf_stats)

    output:
        path("vcf-stats-plots")
    
    script:
    """
        # plot VCF stats 
        plot-vcfstats --prefix vcf-stats-plots --vectors --main-title ${params.plot_vcf_stats.main_title} ${vcf_stats}
    """

    stub:
    """
        mkdir vcf-stats-plots
    """

}

process CONCAT_VCF {

    publishDir "${params.resultsDir}/concat-vcf/", mode: 'copy', overwrite: true

    input:
        path(vcf)

    output:
        path("concat.vcf.gz")
    
    script:
    """
        # index vcf files
        parallel tabix -p vcf ::: ${vcf}

        # merge vcf files
        bcftools concat --threads ${task.cpus} --allow-overlaps -Oz ${vcf} > concat.vcf.gz
    """

    stub:
    """
        touch concat.vcf.gz
    """

}

process VCFR_QC_MERGE {

    publishDir "${params.resultsDir}/vcfr-qc-merge/", mode: 'copy', overwrite: true

    input:
        path(merged_vcf)
        path(gff)
        path(fasta)

    output:
        path("vcfr-chrom-plot.pdf")
        path("vcfr-chromoqc-plot.pdf")

    script:
    """
        cubseq-vcfr-qc.R \\
            --vcf ${merged_vcf} \\
            --gff ${gff} \\
            --fasta ${fasta}
    """

    stub:
    """
        touch vcfr-chrom-plot.pdf
        touch vcfr-chromoqc-plot.pdf
    """

}

process VCFR_QC_CONCAT {

    publishDir "${params.resultsDir}/vcfr-qc-concat/", mode: 'copy', overwrite: true

    input:
        path(concat_vcf)
        path(gff)
        path(fasta)

    output:
        path("vcfr-chrom-plot.pdf")
        path("vcfr-chromoqc-plot.pdf")

    script:
    """
        cubseq-vcfr-qc.R \\
            --vcf ${concat_vcf} \\
            --gff ${gff} \\
            --fasta ${fasta}
    """

    stub:
    """
        touch vcfr-chrom-plot.pdf
        touch vcfr-chromoqc-plot.pdf
    """

}

process PLOT_VCF_QC {

    publishDir "${params.resultsDir}/vcf-qc-plot/", mode: 'copy', overwrite: true

    input:
        path(vcf)

    output:
        path("vcf-qc-raw.pdf")
        path("vcf-qc-99_log10.pdf")
        path("vcf-qc-999_log10.pdf")
        path("vcf-qc-log10-95_log10.pdf")

    script:
    """
        cubseq-vcf-qc.R ${vcf}
    """

    stub:
    """
        touch vcf-qc-raw.pdf
        touch vcf-qc-99_log10.pdf
        touch vcf-qc-999_log10.pdf
        touch vcf-qc-log10-95_log10.pdf
    """

}

process PLOT_MUTATION_CIRCOS {

    publishDir "${params.resultsDir}/mutation-circos-plot/", mode: 'copy', overwrite: true

    input:
        path(mutation_freq_data)
        path(gtf)
        path(gene_tpm_scores)

    output:
        path("circlize_mutations.pdf")

    script:
    """
        cubseq-mutation-panel.R \\
            ${mutation_freq_data} \\
            --gtf ${gtf} \\
            --gene_tpm_scores ${gene_tpm_scores} \\
            --fdr-threshold ${params.rank_genes.fdr}
    """

    stub:
    """
        touch circlize_mutations.pdf
    """

}

process ANNOTATE_VCF_NEW {

    publishDir "${params.resultsDir}/merged-vcf-annotated/", mode: 'copy', overwrite: true

    input:
        path(gtf)
        path(vcf)

    output:
        path("annotated.vcf.gz")

    script:
    """
        # annotate VCF file with gene names
        // bcftools annotate \\
        //     --threads ${task.cpus} \\
        //     -a genes.bed.gz \\
        //     -c CHROM,FROM,TO,GENE \\
        //     -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \\
        //     variants.vcf.gz

        bcftools annotate \\
            --threads ${task.cpus} \\
            --annotations ${gtf} \\
            --columns gene_name \\
            --output output.vcf \\
            ${vcf}

    """

    stub:
    """
        touch annotated.vcf.gz
    """

}

// CU count data
process SUMMARISE_CODON_COUNTS {

    publishDir "${params.resultsDir}/summarise-codon-counts/", mode: 'copy', overwrite: true

    input:
        path(codon_counts_csv)

    output:
        path("summary-counts.csv")
        path("aa-counts.csv")

    script:
    """
        cubseq-summarise-codon-counts.R ${codon_counts_csv}
    """

    stub:
    """
        touch summary-counts.csv
        touch aa-counts.csv
    """

}

// tRNA analysis
process GET_FEATURECOUNTS_TPM {

    publishDir "${params.resultsDir}/featureCounts-tpm/", mode: 'copy', overwrite: true

    tag "${run_acc}"

    input:
        tuple val(run_acc), path(featurecounts_file)
        path(gtf)

    output:
        tuple val(run_acc), path("${run_acc}.featureCounts-tpm.txt")

    script:
    """
        cubseq-get-featurecounts-tpm.R ${featurecounts_file} \\
            --gtf ${gtf} \\
            --output ${run_acc}.featureCounts-tpm.txt
    """

    stub:
    """
        touch ${run_acc}.featureCounts-tpm.txt
    """

}

process SUMMARIZE_TO_GENE_FEATURECOUNTS {

    publishDir "${params.resultsDir}/featurecounts-txi/", mode: 'copy', overwrite: true

    input:
        path(metadata)
        path(gtf)
        path(featurecounts_tpm_dir)

    output:
        path("txi-featurecounts-summarized-experiment.rds")

    script:
    """
        cubseq-summarize-featurecounts-to-gene.R \\
            ${metadata} \\
            --gtf ${gtf} \\
            --featurecounts-tpm-dir ${featurecounts_tpm_dir} \\
            --counts-from-abundance ${params.summarize_to_gene_featurecounts.counts_from_abundance} \\
            --output txi-featurecounts-summarized-experiment.rds
    """

    stub:
    """
        touch txi-featurecounts-summarized-experiment.rds
    """

}


workflow {

    // get metadata
    metadata_ch = GET_METADATA()
                    .splitCsv(header: true, sep: '\t')
                    .multiMap { row -> 
                        metadata: tuple(row.study_accession, row.sample_accession, row.experiment_accession, row.run_accession, row.fastq1, row.fastq2)
                        reads: row.run_accession
                    }

    // get sample fastq
    GET_READS(metadata_ch.metadata)
    PREPROCESS_READS(GET_READS.out)
    
    // perform alignment
    GENERATE_GENOME_INDEX(file(params.genome.reference), file(params.genome.annotation))
    ALIGN_READS(GENERATE_GENOME_INDEX.out, PREPROCESS_READS.out.trimmed_fastq)
    // TODO: add ALIGNMENT_QC process

    // count number reads that map to gene in genome annotation file
    QUANTIFY_READS(file(params.genome.annotation), ALIGN_READS.out)

    // variant calling
    FREEBAYES_CALL_VARIANTS(file(params.genome.reference), ALIGN_READS.out)

    // filter vcf
    BCFTOOLS_FILTER_VCF(FREEBAYES_CALL_VARIANTS.out)

    // generate consensus fasta file
    BCFTOOLS_CREATE_CONSENSUS(file(params.genome.reference), BCFTOOLS_FILTER_VCF.out)

    // generate transcriptome fasta for wildtype and mutated sequences
    GFFREAD_GET_WT_TRANSCRIPTOME(file(params.genome.reference), file(params.genome.annotation))
    GFFREAD_GET_MUT_TRANSCRIPTOME(BCFTOOLS_CREATE_CONSENSUS.out, file(params.genome.annotation))

    // create transcriptome index
    SALMON_INDEX(GFFREAD_GET_WT_TRANSCRIPTOME.out.transcriptome, file(params.genome.reference))

    // estimate transcript level abundance
    SALMON_QUANT(SALMON_INDEX.out, PREPROCESS_READS.out.trimmed_fastq)

    // summarise transcript-level abundance estimates to gene level
    SUMMARIZE_TO_GENE(GET_METADATA.out, file(params.genome.annotation), SALMON_QUANT.out.collect())

}

workflow SALMON_REDO {

    // construct file paths using run_accession from metadata file
    metadata_ch = Channel
                    .fromPath(params.metadata, checkIfExists: true)
                    .splitCsv(header: true, sep: '\t')
                    .multiMap { row -> 
                        metadata: tuple(row.study_accession, row.sample_accession, row.experiment_accession, row.run_accession, row.fastq1, row.fastq2)
                        reads: row.run_accession
                    }

    // get sample fastq
    GET_READS(metadata_ch.metadata)
    PREPROCESS_READS(GET_READS.out)

    // generate transcriptome fasta for wildtype and mutated sequences
    GFFREAD_GET_WT_TRANSCRIPTOME(file(params.genome.reference), file(params.genome.annotation))

    // create transcriptome index
    SALMON_INDEX(GFFREAD_GET_WT_TRANSCRIPTOME.out.transcriptome, file(params.genome.reference))

    // estimate transcript level abundance
    SALMON_QUANT(SALMON_INDEX.out, PREPROCESS_READS.out.trimmed_fastq)

    // summarise transcript-level abundance estimates to gene level
    SUMMARIZE_TO_GENE(file(params.metadata), file(params.genome.annotation), SALMON_QUANT.out.collect())

}

workflow CUB_PREPROCESSING {

    // summarise transcript-level abundance estimates to gene level
    // SUMMARIZE_TO_GENE(file(params.metadata), file(params.genome.annotation), file(params.summarize_to_gene.quant.dir))

    // rank genes using Fisher transform
    // RANK_GENES(SUMMARIZE_TO_GENE.out, file(params.metadata), file(params.genome.annotation), file(params.gff))
    RANK_GENES(file(params.txi), file(params.metadata), file(params.genome.annotation), file(params.gff))

    // filter FASTA files for significantly highly and lowly expressed genes
    GET_HEG_LEG_FASTA(file(params.metadata), file(params.genome.annotation), file(params.mut_fasta), RANK_GENES.out.heg_scores, RANK_GENES.out.leg_scores)
    // GET_HEG_LEG_FASTA(file(params.get_heg_leg_fasta.metadata), file(params.genome.annotation), file(params.mut_fasta), file(params.get_heg_leg_fasta.heg_scores), file(params.get_heg_leg_fasta.leg_scores))

    // calculate relative codon frequencies
    // GET_CU_FREQUENCIES(GET_HEG_LEG_FASTA.out.heg.collect(), GET_HEG_LEG_FASTA.out.leg.collect())
    // GET_CU_FREQUENCIES(GET_HEG_LEG_FASTA.out.heg, GET_HEG_LEG_FASTA.out.leg)
    // GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir))

    // testing
    // heg_fasta_ch = GET_HEG_LEG_FASTA.out.heg.collect()
    // leg_fasta_ch = GET_HEG_LEG_FASTA.out.leg.collect()
    // fileNamesChannel = EmitFiles.out.fileNames.collectFile()
    // GET_CU_FREQUENCIES(GET_HEG_LEG_FASTA.out.heg.collect{it}, GET_HEG_LEG_FASTA.out.leg.collect{it})

    // GET_HEG_LEG_FASTA.out.heg.collectFile().view()
    // heg_fasta_ch.view()
    // leg_fasta_ch.view()

    // // generate CU figures
    // GET_CU_FIGURES(GET_CU_FREQUENCIES.out.rf, file(params.get_cu_figures.aa_property_table))
    // // GET_CU_FIGURES(file(params.get_cu_figures.rf), file(params.get_cu_figures.aa_property_table))

}

workflow GET_CU {

    // calculate relative codon frequencies
    GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir))

    // generate CU figures
    GET_CU_FIGURES(GET_CU_FREQUENCIES.out.rf, file(params.get_cu_figures.aa_property_table))

}

workflow GET_GLOBAL_CUB {

    // filter fasta files to keep protein-coding genes
    GET_PROTEIN_FASTA(file(params.metadata), file(params.genome.annotation), file(params.mut_fasta))

    // calculate relative codon frequencies
    GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir), file(params.get_cu_frequencies.protein_fasta_dir))
    // GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir), GET_PROTEIN_FASTA.out[0])

}

workflow CUB_PREPROCESSING_REDO {

    // construct file paths using run_accession from metadata file
    metadata_ch = Channel
                    .fromPath(params.metadata, checkIfExists: true)
                    .splitCsv(header: true, sep: '\t')
                    .multiMap { row -> 
                        metadata: tuple(row.study_accession, row.sample_accession, row.experiment_accession, row.run_accession, row.fastq1, row.fastq2)
                        reads: row.run_accession
                    }

    // get sample fastq
    GET_READS(metadata_ch.metadata)
    PREPROCESS_READS(GET_READS.out)

    // create transcriptome index
    SALMON_INDEX(GFFREAD_GET_WT_TRANSCRIPTOME.out.transcriptome, file(params.genome.reference))

    // salmon quant - include ALL reads (tRNAs)
    SALMON_QUANT(SALMON_INDEX.out, metadata_ch.reads)

    // // summarise to gene level - tximport
    // SUMMARIZE_TO_GENE(file(params.metadata), file(params.genome.annotation), file(params.summarize_to_gene.quant.dir))

    // // rank genes

    // // filter FASTA files for significant highly expressed genes
    // GET_HEG_LEG_FASTA(file(params.get_heg_leg_fasta.metadata), file(params.genome.annotation), file(params.mut_fasta), file(params.get_heg_leg_fasta.heg_scores), file(params.get_heg_leg_fasta.leg_scores))

    // // calculate relative codon frequencies
    // GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir))

}

// workflow CUB_PREPROCESSING {

//     // summarise transcript-level abundance estimates to gene level
//     SUMMARIZE_TO_GENE(file(params.metadata), file(params.genome.annotation), file(params.summarize_to_gene.quant.dir))

//     // get significant HEG and LEG transcript IDs
//     //GET_HEG_LEG(SUMMARIZE_TO_GENE.out, file(params.genome.annotation))

//     // get concatenated fasta file for HEG and LEG
//     // GET_HEG_LEG_FASTA(metadata_ch.metadata, file(params.top.genes.id), GFFREAD_GET_MUT_TRANSCRIPTOME.out)
//     //GET_HEG_LEG_FASTA(file(params.metadata), file(params.mut_fasta), GET_HEG_LEG.out.heg, GET_HEG_LEG.out.leg)

//     // // Annotate fasta files with KO ID
//     // GET_KO_ID(GFFREAD_GET_MUT_TRANSCRIPTOME.out.transcriptome)

//     // // count codons
//     // CORDON_COUNT_CODONS(GFFREAD_GET_MUT_TRANSCRIPTOME.out.collect())

//     // calculate relative codon frequencies
//     // GET_CU_FREQUENCIES(GET_HEG_LEG_FASTA.out.heg, GET_HEG_LEG_FASTA.out.leg)
//     GET_CU_FREQUENCIES(file(params.get_cu_frequencies.heg_fasta_dir), file(params.get_cu_frequencies.leg_fasta_dir))

//     // generate CU figures
//     GET_CU_FIGURES(GET_CU_FREQUENCIES.out.rf, file(params.get_cu_figures.aa_property_table))

//     // generate pca figures
//     //GET_PCA(SUMMARIZE_TO_GENE.out, file(params.metadata))

// }

workflow PCA {

    GET_PCA(file(params.get_pca.txi), file(params.get_pca.metadata))

}

// TODO: can use nextflow "splitFasta" and "filter" operators to filter out HE genes

workflow STRAIN_ANALYSIS {

    // we associate each run_accession with its strain (the key)
    // we construct file paths for each sample fasta file
    // then group file paths by strain
    // use resulting channel as input for process

    // TODO: update so that workflow takes in filtered metadata

    // construct path to HEG fasta files and group by strains using metadata
    // strain_ch = Channel
    //     .fromPath(params.strain.test.metadata, checkIfExists: true)
    //     .splitCsv(header: true, sep: '\t')
    //     .multiMap { row ->
    //         fasta_files: tuple(row.strain, file(params.get_cu_frequencies.heg_fasta_dir + "/" + row.run_accession + ".heg-mut-transcriptome.fa", checkIfExists: true))
    //         quant_files: tuple(row.strain, file(params.strain.test.quant_dir + "/" + row.run_accession + "/quant.sf", checkIfExists: true))
    //     }

    // strain_ch
    //     .groupTuple()
    //     .view()

    // STRAIN_CU_CAT_PLOT(file(params.metadata), file(params.genome.annotation), strain_ch.quant_files.groupTuple())

    // strain_ch.quant_files.groupTuple().view()

    // STRAIN_HEG_CAT_PLOT(file(params.real.txi), file(params.real.metadata), file(params.genome.annotation))

    STRAIN_CU_HEATMAP(file(params.strain_cu_heatmap.metadata), file(params.strain_cu_heatmap.fasta))

}

workflow VARIANT_ANALYSIS {

    // construct path to vcf files using run_accession from metadata file
    vcf_ch = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.run_accession, file(params.vcf_dir + "/" + row.run_accession + ".filter.norm.var.vcf.gz", checkIfExists: true))
        }

    // vcf_ch.view()

    CALL_CONSEQUENCE(file(params.genome.reference), file(params.gff), vcf_ch)

}

workflow MERGED_CSQ {

    CALL_CONSEQUENCE_ON_MERGED(file(params.genome.reference), file(params.gff), file(params.vcf_dir))

}

workflow VCF_PROCESSING {

    // construct path to vcf files using run_accession from metadata file
    // vcf_ch = Channel
    //     .fromPath(params.metadata, checkIfExists: true)
    //     .splitCsv(header: true, sep: '\t')
    //     .map { row ->
    //         tuple(row.run_accession, file(params.vcf + "/" + row.run_accession + ".filter.norm.var.vcf.gz", checkIfExists: true))
    //     }

    vcf_ch = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            file(params.vcf + "/" + row.run_accession + ".filter.norm.var.vcf.gz", checkIfExists: true)
        }

    // MERGE_VCF(vcf_ch.collect())
    // VCFR_QC_MERGE(MERGE_VCF.out, file(params.gff), file(params.genome.reference))

    // GET_VCF_STATS(MERGE_VCF.out, file(params.genome.reference))
    GET_VCF_STATS(file(params.get_vcf_stats.merged_vcf), file(params.genome.reference))
    PLOT_VCF_STATS(GET_VCF_STATS.out)

    // // CONCAT_VCF(vcf_ch.collect())
    // // VCFR_QC_CONCAT(CONCAT_VCF.out, file(params.gff), file(params.genome.reference))

    // VCFR_QC(MERGE_VCF.out, file(params.gff), file(params.genome.reference))
    // VCFR_QC(file(params.vcfr_qc.merged_vcf), file(params.gff), file(params.genome.reference))

    // ANNOTATE_VCF(file(params.genome.annotation), vcf_ch)

    // ANNOTATE_VCF.out.map{ it -> [it, it] }.view()

    // PLOT_TOP_MUTATED_GENES(file(params.metadata), file(params.annotate_vcf.dir))

}

workflow VCF {

    vcf_ch = Channel
        .fromPath(params.metadata, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            file(params.vcf + "/" + row.run_accession + ".filter.norm.var.vcf.gz", checkIfExists: true)
        }

    MERGE_VCF(vcf_ch.collect())

    PLOT_VCF_QC(MERGE_VCF.out)
    
}

workflow VCF_DROP_GENOTYPES {

    DROP_VCF_GENOTYPES(file(params.bcftools_view.merged_vcf))

}

workflow PLOT_VCF {

    PLOT_VCF_QC(file(params.plot_vcf.merged_vcf))

}

workflow MUTATION_PANEL {

    PLOT_TOP_MUTATED_GENES(file(params.metadata), file(params.annotate_vcf.dir))
    
}

workflow MUTATION_PLOT {

    // GET_MUTATION_PANEL_DATA(file(params.plot_vcf.merged_vcf))

    PLOT_MUTATION_CIRCOS(file(params.mutation_freq_data), file(params.genome.annotation), file(params.gene_tpm_scores))

}

workflow SUMMARISE_CC {

    SUMMARISE_CODON_COUNTS(file(params.summarise_cc.codon_counts_csv))

}

// workflow DEBUG_BCFTOOLS {

//     // construct path to vcf files using run_accession from metadata file
//     vcf_ch = Channel
//         .fromPath(params.debug_bcftools.metadata, checkIfExists: true)
//         .splitCsv(header: true, sep: '\t')
//         .map { row ->
//             tuple(row.study_accession, row.sample_accession, row.experiment_accession, row.run_accession, file(params.debug_bcftools.vcf + "/" + row.run_accession + ".var.vcf", checkIfExists: true))
//         }

//     // filter vcf
//     BCFTOOLS_FILTER_VCF(vcf_ch)

// }

// workflow CORRECT_SNV {

//     // This process is to correctly filter MNVs and only retain true SNVs (i.e. not to decompose MNVs to SNVs).

//     // construct path to BAM files using run_accession from metadata file
//     bam_ch = Channel
//         .fromPath(params.metadata, checkIfExists: true)
//         .splitCsv(header: true, sep: '\t')
//         .map { row ->
//             tuple(row.study_accession, row.sample_accession, row.experiment_accession, row.run_accession, file(params.correct_snv.bam + "/" + row.run_accession + ".Aligned.sortedByCoord.out.bam", checkIfExists: true))
//         }
    
//     // bam_ch.view()

//     // variant calling
//     FREEBAYES_CALL_VARIANTS(file(params.genome.reference), bam_ch)

//     // filter vcf
//     SNV_BCFTOOLS_FILTER_VCF(FREEBAYES_CALL_VARIANTS.out)

//     // generate consensus fasta file
//     BCFTOOLS_CREATE_CONSENSUS(file(params.genome.reference), SNV_BCFTOOLS_FILTER_VCF.out)

//     // generate transcriptome fasta for wildtype and mutated sequences
//     GFFREAD_GET_WT_TRANSCRIPTOME(file(params.genome.reference), file(params.genome.annotation))
//     GFFREAD_GET_MUT_TRANSCRIPTOME(BCFTOOLS_CREATE_CONSENSUS.out, file(params.genome.annotation))

// }

// tRNA analysis
workflow TRNA_ANALYSIS {

    // construct path to vcf files using run_accession from metadata file
    featurecounts_ch = Channel
        .fromPath(params.featurecounts_metadata, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.run_accession, file(params.featurecounts_dir + "/" + row.run_accession + "/" + row.run_accession + ".featureCounts.txt", checkIfExists: true))
        }
    
    // featurecounts_ch.view()

    GET_FEATURECOUNTS_TPM(featurecounts_ch, file(params.genome.annotation))

    // summarise transcript-level abundance estimates to gene level
    SUMMARIZE_TO_GENE_FEATURECOUNTS(file(params.featurecounts_metadata), file(params.genome.annotation), file(params.featurecounts_tpm_dir))

}
