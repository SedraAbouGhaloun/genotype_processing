# LOCAL CONFIGURRATION:

plink2 = "/sc-resources/ukb/data/shared/open/software/plink2"
plink = "/sc-resources/ukb/data/shared/open/software/plink"



## CHROM = range(1, 23)

rule all:
    input:
        expand("{out_path}/r2/{ancestry_name}/chr_{chr}_{r2}.ld",allow_missing=True,out_path=str(config["out"]),ancestry_name=(config["ancestry_name"]), r2=str(config["r2"]),chr=config["CHROM"])
    params:
        output_tmp_file=expand("{out_path}", allow_missing=True, out_path=str(config["out"]))
    shell:
        """
       ## rm -r  {params.output_tmp_file}/tmp_chr
       ## rm -r {params.output_tmp_file}/bgen11
        """
        
        
rule separate_chr:
    input:
        snps_file=config["snps_file"]
    output:
        tmp_chrs=expand("{out_path}/tmp_chr/chr{{chr}}.csv", allow_missing=True, out_path=str(config["out"]))
    log:
        expand("{out_path}/logs/separate_chr{{chr}}.log", allow_missing=True, out_path=str(config["out"]))

    params:
        output_tmp_file=expand("{out_path}", allow_missing=True, out_path=str(config["out"]))
    shell:
        """
        mkdir -p {params.output_tmp_file}/tmp_chr
        python umbrella/pre/gwas/scripts/split_chr.py --i_path {input.snps_file} --o_path {params.output_tmp_file}/tmp_chr/
        """

rule write_bgen_v11:
    threads: 1
    input:
        bgen=expand("{genetic_data}/ukb22828_c{{chr}}_b0_v3.bgen",genetic_data=config["genetic_data"]),
        sample_file=expand("{genetic_data}/ukb22828_c10_b0_v3.sample",genetic_data=config["genetic_data"]),
        var_subset=expand("{out_path}/tmp_chr/chr{{chr}}.csv",allow_missing=True,out_path=str(config["out"]))
    output:
        bgen=expand("{out_path}/bgen11/chr_{{chr}}.bgen",allow_missing=True,out_path=str(config["out"])),
        sample=expand("{out_path}/bgen11/chr_{{chr}}.sample", allow_missing = True, out_path = str(config["out"]))
    log:
        expand("{out_path}/logs/write_bgen_v11_{{chr}}.log",allow_missing=True,out_path=str(config["out"]))
    params:
        out_path = str(config["out"])
    shell:
        """
        mkdir -p {params.out_path}/bgen11
        {plink2} --bgen {input.bgen} ref-last --sample {input.sample_file} --export bgen-1.1 --extract {input.var_subset} \
        --out {params.out_path}/bgen11/chr_{wildcards.chr} --debug >> {log} # --memory 50000  --threads 1
        """



rule calculate_LD:
    threads: 1
    input:
        bgen=expand("{out_path}/bgen11/chr_{{chr}}.bgen",allow_missing=True,out_path=str(config["out"])),
        sample_file=expand("{genetic_data}/ukb22828_c10_b0_v3.sample",genetic_data=config["genetic_data"])
    log:
        expand("{out_path}/logs/calculate_LD_{{ancestry_name}}_{{chr}}.log",allow_missing=True,out_path=str(config["out"]))
    output:
        expand("{out_path}/r2/{{ancestry_name}}/chr_{{chr}}_{r2}.ld",allow_missing=True,out_path=str(config["out"]),r2=str(config["r2"]))
    params:
        ld_window = config["ld_window_kb"],
        r2 = config["r2"],
        ancestry_path = config["ancestry_path"],
        out_path = str(config["out"])
    shell: 
        """
        mkdir -p {params.out_path}/r2/{wildcards.ancestry_name}
        {plink} --bgen {input.bgen} --sample {input.sample_file} --keep {params.ancestry_path}/{wildcards.ancestry_name}.txt --r2 --ld-window-kb {params.ld_window} --ld-window-r2 {params.r2}  --out {params.out_path}/r2/{wildcards.ancestry_name}/chr_{wildcards.chr}_{params.r2} --debug >> {log} 
        """
        
