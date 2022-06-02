#!/usr/bin/env nextflow

/*
Assemble Illumina short reads using Unicycler.

[Use guide]
To run this pipeline in a screen session:
    nextflow -Djava.io.tmpdir=$PWD run unicycler.nf --fastq "./reads/*_{1,2}.fastq.gz" \
                              --assemblyDir assembly --condaUnicycler "$HOME/anaconda3/envs/unicycler" \
                              -c unicycler.config -profile pbs

[Declaration]
Copyright (C) 2020-2022 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public License v3.0
Publication: 28 March 2020; latest update: 20 May 2022
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

/*------------------------------------------------------------------------------
                           P R O C E S S E S 
------------------------------------------------------------------------------*/
process Unicycler {
    publishDir "${outdir}", pattern: "assembly.fasta", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.fasta" }
    publishDir "${outdir}", pattern: "assembly.gfa", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.gfa" }
    publishDir "${outdir}", pattern: "unicycler.log", mode: "copy", overwrite: true, saveAs: { filename -> "${genome}.log" }
    
    input:
    tuple val(genome), file(fastqs)

    output:
    file("assembly.fasta")
    file("assembly.gfa")
    file("unicycler.log")
    
    script:    
    """
    module load anaconda3/personal
    source activate unicycler0.5.0
    unicycler -1 ${fastqs[0]} -2 ${fastqs[1]} --mode normal --threads ${params.cpus} --keep 0 --out .
    """
}

/*------------------------------------------------------------------------------
                           Main 
------------------------------------------------------------------------------*/
workflow {
    readsets = Channel.fromFilePairs(params.fastq)
    Unicycler(readsets)
}

/* References
1. https://github.com/nf-core/denovohybrid
2. https://github.com/nextflow-io/patterns/blob/master/docs/publish-rename-outputs.adoc
3. https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir
*/