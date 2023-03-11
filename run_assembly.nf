#!/usr/bin/env nextflow

/*
Assemble Illumina short reads using Unicycler or SPAdes.

[Use guide]
To run this pipeline in a screen session:
    nextflow -Djava.io.tmpdir=$PWD run run_assembly.nf --fastq "./reads/*_{1,2}.fastq.gz" --outdir assembly \
    --conda_unicycler "unicycler0.5.0" -c run_assembly.config -profile pbs

[Declaration]
Copyright (C) 2020-2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public License v3.0
Publication: 16 June 2022; latest update: 11 March 2023
*/

/*------------------------------------------------------------------------------
                       C O N F I G U R A T I O N
------------------------------------------------------------------------------*/

nextflow.enable.dsl = 2

def mkdir(dir_path) {  // Creates a directory and returns a File object
    def dir_obj = new File(dir_path)
    if ( !dir_obj.exists() ) {
        result = dir_obj.mkdir()
        println result ? "Successfully created directory ${dir_path}" : "Cannot create directory ${dir_path}"
    } else {
        println "Directory ${dir_path} exists."
    }
    return dir_obj
}

outdir = mkdir(params.outdir)
logdir = mkdir(params.outdir + "/log")

/*------------------------------------------------------------------------------
                           P R O C E S S E S 
------------------------------------------------------------------------------*/
process Unicycler {
    publishDir "${outdir}", pattern: "assembly.fasta", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.fna" }
    publishDir "${outdir}", pattern: "assembly.gfa", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.gfa" }
    publishDir "${outdir}/log", pattern: "unicycler.log", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.log" }
    
    input:
    tuple val(genome), file(fastqs)

    output:
    file("assembly.fasta")
    file("assembly.gfa")
    file("unicycler.log")
    
    script:    
    """
    module load anaconda3/personal
    source activate ${params.conda_unicycler}
    unicycler -1 ${fastqs[0]} -2 ${fastqs[1]} --mode normal --threads ${params.cpus} --keep 0 --out .
    """
}

process SPAdes {
    publishDir "${outdir}", pattern: "output/scaffolds.fasta", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}__scaffolds.fna" }
    publishDir "${outdir}", pattern: "output/assembly_graph_with_scaffolds.gfa", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}__scaffolds.gfa" }
    publishDir "${outdir}", pattern: "output/scaffolds.paths", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}__scaffolds.paths" }
    publishDir "${outdir}/log", pattern: "output/spades.log", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.log" }

    input:
    tuple val(genome), file(fastqs)

    output:
    file("output/scaffolds.fasta")
    file("output/assembly_graph_with_scaffolds.gfa")
    file("output/scaffolds.paths")
    file("output/spades.log")
    
    script:    
    """
    module load anaconda3/personal
    source activate ${params.conda_spades}
    spades.py -1 ${fastqs[0]} -2 ${fastqs[1]} -o output --phred-offset 33 --${params.spades_mode} --threads ${params.cpus} -k '${params.spades_kmers}'
    """
}

/*------------------------------------------------------------------------------
                           Main 
------------------------------------------------------------------------------*/
workflow {
    readsets = Channel.fromFilePairs(params.fastq)
    if (params.assembler == "unicycler") {
        println "Use Unicycler to assemble genomes"
        Unicycler(readsets)
    } else {
        println "Use SPAdes to assemble genomes"
        SPAdes(readsets)
    }
}

/* References
1. https://github.com/nf-core/denovohybrid
2. https://github.com/nextflow-io/patterns/blob/master/docs/publish-rename-outputs.adoc
3. https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir
*/