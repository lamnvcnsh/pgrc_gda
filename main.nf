params.input_dir = 'data/*'
params.output_dir =  'output'

process processFiles {
    tag "Process $vcf"
    input:
    path vcf

    output:
    val "${vcf.baseName}"

    script:
    """
    echo ${vcf.baseName}
    """
}

process addChr {
    // publishDir "${params.output_dir}/add_chr", mode: 'copy'
    input:
    path vcf
    
    output:
    path "${vcf.baseName}_addchr.vcf", emit: vcf_chr

    """
    awk '!/^#/ {print "chr"\$0; next} 1' ${vcf} > ${vcf.baseName}_addchr.vcf
    """
}


process runPharmcat {
    container 'pgkb/pharmcat:latest'
    publishDir "${params.output_dir}/pharmcat", mode: 'copy'

    input:
    // path vcf_file from file(params.vcf_file)
    path vcf_file

    output:
    path "${vcf_file.baseName}"

    """
    mkdir -p /pharmcat/data
    cp -r $vcf_file /pharmcat/data/
    /pharmcat/pharmcat_pipeline /pharmcat/data/${vcf_file.name} -o ${vcf_file.baseName}
    """
}


// process runStargzer {
//     publishDir "${params.output_dir}/pharmcat", mode: 'copy'
//     input:
//     path vcf_file

// }


workflow {
    ch = Channel.fromPath(params.input_dir)
    add_chr = addChr(ch)
    runPharmcat(add_chr)
}
