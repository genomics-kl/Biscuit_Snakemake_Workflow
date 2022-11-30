rule biscuiteer:
    input:
        mergecg_gz = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz',
        mergecg_tbi = f'{output_directory}/analysis/pileup/{{sample}}_mergecg.bed.gz.tbi',
        vcf = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz',
        vcf_tbi = f'{output_directory}/analysis/pileup/{{sample}}.vcf.gz.tbi',
    output:
        #raw=f'{output_directory}/analysis/biscuiteer/{{sample}}.bsseq.raw.rds',
        f'{output_directory}/analysis/biscuiteer/{{sample}}.bsseq.rds'
    params:
        genome=config['biscuiteer']['genome']
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
        Rscript --vanilla -e "library(biscuiteer); library(bsseq); library(BiocParallel); bisc <- sort(readBiscuit(BEDfile='{input.mergecg_gz}', VCFfile='{input.vcf}', merged=TRUE, genome='{params.genome}')); saveRDS(bisc, '{output}'); sessionInfo();" 1> {log.out} 2> {log.err}
        """

rule unionize_bsseq:
    input:
        expand('{output_dir}/analysis/biscuiteer/{samples.sample}.bsseq.rds', samples = SAMPLES.itertuples(), output_dir=output_directory)
    output:
        smooth=f'{output_directory}/analysis/unionize_bsseq/all.bsseq.smooth.rds'
    params:
        bsseqs = lambda wildcards, input: ",".join("'" + x + "'" for x in input),
        #output_raw=f'{output_directory}/all.bsseq.raw.rds'
    log:
        out = f'{output_directory}/logs/unionize_bsseq/out.o',
        err = f'{output_directory}/logs/unionize_bsseq/out.e',
    benchmark:
        f'{output_directory}/benchmarks/unionize_bsseq/bench.txt',
    threads: 16
    resources:
        mem_gb = config['hpcParameters']['maxMemoryGb'],
        walltime = config['walltime']['medium'],
    conda:
        '../envs/biscuit.yaml'
    envmodules:
        config['envmodules']['biscuiteer']
    shell:
        """
        Rscript --vanilla -e "library(biscuiteer); library(bsseq); library(BiocParallel); all_bisc <- do.call(unionize, lapply(c({params.bsseqs}), readRDS)); all_bisc.smooth <- BSmooth(BSseq = all_bisc, BPPARAM = MulticoreParam(workers = '{threads}'), verbose = TRUE, keep.se = TRUE); saveRDS(all_bisc.smooth, '{output.smooth}', compress=FALSE); sessionInfo();" 1> {log.out} 2> {log.err}

        """
