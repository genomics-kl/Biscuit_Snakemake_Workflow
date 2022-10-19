rule biscuiteer:
    input:
        mergecg_gz = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
        mergecg_tbi = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz.tbi',
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
        vcf_tbi = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz.tbi',
    output:
        f'{output_directory}/analysis/biscuiteer/{{sample}}.bsseq.rds'
    params:
    log:
        out = f'{output_directory}/logs/biscuiteer/{{sample}}.o',
        err = f'{output_directory}/logs/biscuiteer/{{sample}}.e',
    benchmark:
        f'{output_directory}/benchmarks/biscuiteer/{{sample}}.txt',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuiteer']
    shell:
        """
        Rscript --vanilla -e "library(biscuiteer); saveRDS(readBiscuit(BEDfile='{input.mergecg_gz}', VCFfile='{input.vcf}', merged=TRUE), '{output}'); sessionInfo();" 1> {log.out} 2> {log.err}
        """

rule unionize_bsseq:
    input:
        expand(f'{output_directory}/analysis/biscuiteer/{{samples.sample}}.bsseq.rds', samples = SAMPLES.itertuples())
    output:
        f'{output_directory}/analysis/biscuiteer/all.bsseq.rds'
    params:
        bsseqs = lambda wildcards, input: ",".join("'" + x + "'" for x in input)
    log:
        out = f'{output_directory}/logs/unionize_bsseq/out.o',
        err = f'{output_directory}/logs/unionize_bsseq/err.e',
    benchmark:
        f'{output_directory}/benchmarks/unionize_bsseq/bench.txt',
    threads: 8
    resources:
        mem_gb = config['hpcParameters']['intermediateMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuiteer']
    shell:
        """
        Rscript --vanilla -e "library(biscuiteer); saveRDS(do.call(unionize, lapply(c({params.bsseqs}), readRDS)), '{output}'); sessionInfo();" 1> {log.out} 2> {log.err}
        """
