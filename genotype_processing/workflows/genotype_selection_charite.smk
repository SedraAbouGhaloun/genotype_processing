# LOCAL CONFIGURRATION:
SNPEFFJAR = '/sc-resources/ukb/data/shared/open/software/snpEff/snpEff.jar'
plink2 = '/sc-resources/ukb/data/shared/open/software/plink2'

GENOMEBUILD = 'GRCh37.75'

NFEATURES = ['10k', '100k', '250k', '500k', '2500k', '1M', '5M']
SPLIT = ["train", "test", "valid", "debug", "extended_tr", "outliers_te"]

rule all:
    threads: 1
    input:
        expand("{out_path}/plink/plink_temp_chr_{CHROM}.{nfeatures}.raw", allow_missing = True, out_path = str(config["out"]), CHROM=config["CHROM"], nfeatures=str(config["nfeatures"])),



#########################################################################################################################

rule symlink_basics:
    output:
        expand("{out_path}/pvar_ex.tsv", allow_missing=True,out_path=str(config["out"]))
    params:
        variant_list = str(config["variant_list"])
    shell:
        """
        ln -s {params.variant_list} {output} | true
        """

#########################################################################################################################

rule download_prs:
    output:
        expand("{out_path}/{pgsid}.txt", allow_missing=True,out_path=str(config["out"]), pgsid=str(config["pgsid"]))
    log:
        expand("{out_path}/logs/download_prs.log",allow_missing=True,out_path=str(config["out"]))
    params:
        pgsid=str(config["pgsid"]),
        out_dir=str(config["out"])
    shell:
        """
        mkdir -p {params.out_dir} | true
        wget -P {params.out_dir} http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{params.pgsid}/ScoringFiles/Harmonized/{params.pgsid}_hmPOS_GRCh37.txt.gz > {log}

        gzip -d {params.out_dir}/{params.pgsid}_hmPOS_GRCh37.txt.gz
        mv {params.out_dir}/{params.pgsid}_hmPOS_GRCh37.txt {params.out_dir}/{params.pgsid}.txt
        """

#########################################################################################################################


rule rsids_for_prs:
    input:
        pgs=expand("{out_path}/{pgsid}.txt", allow_missing=True,out_path=str(config["out"]), pgsid=str(config["pgsid"])),
        pvar_ex=expand("{out_path}/pvar_ex.tsv", allow_missing=True,out_path=str(config["out"]))
    output:
        snpefflike=expand("{out_path}/{pgsid}.txt.snpefflike", allow_missing=True,out_path=str(config["out"]), pgsid=str(config["pgsid"])),
        subset=expand("{out_path}/subsets_pgs/variant_subset.{nfeatures}.tsv", out_path = config["out"], nfeatures=config["nfeatures"])
    log:
        expand("{out_path}/logs/rsids_for_prs.log",allow_missing=True,out_path=str(config["out"]))
    params:
        out_path=str(config["out"]),
        nfeatures=str(config["nfeatures"])
    resources:
        mem_mb=100_000, # MB
        runtime=1440, # 48 hours
        slurm_extra= '-o $HOME/logs/smk_rfp_%j.out -e $HOME/logs/smk_rfp_%j.out --cpus-per-task=20'
    shell:
        """
        ln -s $HOME/logs/smk_rfp_$SLURM_JOB_ID.out {log}.$SLURM_JOB_ID.out | true
        ln -s $HOME/logs/smk_rfp_$SLURM_JOB_ID.err {log}.$SLURM_JOB_ID.err | true

        mkdir -p {params.out_path}/subsets_pgs/ | true
        python -m scripts.rsids_for_prs --in_pgs_file {input.pgs} --in_pvar_extended_ids {input.pvar_ex} --out {params.out_path} --nvar {params.nfeatures} > {log}
        """


#########################################################################################################################

rule annotate_pgs:
    threads: 1
    input:
        expand("{out_path}/subsets{subset_suffix}/variant_subset.{{suffix}}.tsv", out_path = config["out"], subset_suffix=config["subset_suffix"])
    params:
        jar = expand("{snpeff}", snpeff = SNPEFFJAR),
        out = expand("{out_path}/annotated", out_path = str(config["out"])),
        out_path = str(config["out"]),
        genome = expand("{genome}", genome = GENOMEBUILD),
    output:
        expand("{out_path}/annotated/snpeff_temp.{{suffix}}.txt", allow_missing = True, out_path = config["out"])
    log:
         expand("{out_path}/logs/annotated/snpeff_temp.{{suffix}}.txt.log", allow_missing = True, out_path = config["out"])
    resources:
        mem_mb=100_000, # MB
        runtime=1440, # 48 hours
        slurm_extra= '-o $HOME/logs/smk_snf_%j.out -e $HOME/logs/smk_snf_%j.out --cpus-per-task=20'
    shell:
        """
        ln -s $HOME/logs/smk_snf_$SLURM_JOB_ID.out {log}.$SLURM_JOB_ID.out | true
        ln -s $HOME/logs/smk_snf_$SLURM_JOB_ID.err {log}.$SLURM_JOB_ID.err | true

        mkdir -p {params.out}
        java -Xmx10g -jar {params.jar} -v {params.genome} {input} > {params.out}/snpeff_temp.{wildcards.suffix}.txt
        """


#########################################################################################################################


rule adjust_annotated_file:
    threads: 1
    input:
        expand("{out_path}/annotated/snpeff_temp.{{suffix}}.txt", out_path = str(config["out"]))
    output:
        expand("{out_path}/annotated_df/snpeff_df.{{suffix}}.txt", out_path = str(config["out"]))
    log:
        expand("{out_path}/logs/annotated_df/snpeff_df.{{suffix}}.txt.log", out_path = str(config["out"]))
    shell:
        """
        2>{log}
        python -m scripts.adjust_annotated_file --input {input} --output {output}
        """


#########################################################################################################################



rule extract_chr_specific_rsids:
    threads: 1
    input:
        expand("{out_path}/annotated_df/snpeff_df.{{suffix}}.txt", allow_missing = True, out_path = str(config["out"]))
    output:
        ID = expand("{out}/rsids/pgs_rsids_chr_{chrom}.{{suffix}}.csv", chrom = config["CHROM"], out = str(config["out"])),
        altID = expand("{out}/rsids/pgs_rsids_alt_chr_{chrom}.{{suffix}}.csv",chrom = config["CHROM"], out = str(config["out"]))
    log:
        expand("{out}/logs/rsids/pgs_rsids_chrs.{{suffix}}.csv.log", out = str(config["out"])),
    shell:
        """
        2>{log}
        python -m scripts.extract_rsids --df {input} --outID {output.ID} --outaltID {output.altID}
        """


#########################################################################################################################

rule plink_for_chr_1:
    threads: 1
    input:
        bgen=expand("{genetic_data}/ukb22828_c{{CHROM}}_b0_v3.bgen",genetic_data=config["genetic_data"]),
        sample_file=expand("{genetic_data}/ukb22828_c{{CHROM}}_b0_v3.sample",genetic_data=config["genetic_data"]),
        ID = expand("{out_path}/rsids/pgs_rsids_chr_{{CHROM}}.{{suffix}}.csv", out_path = str(config["out"]))
    output:
        expand("{out_path}/plink/plink_temp_chr_{{CHROM}}.{{suffix}}.pgen", allow_missing = True, out_path = str(config["out"]))
    log:
        expand("{out_path}/logs/plink/plink_temp_chr_{{CHROM}}.{{suffix}}.raw.log", allow_missing = True, out_path = str(config["out"]))
    params:
        pgen_out_path=expand("{out_path}/plink/plink_temp_chr_{{CHROM}}.{{suffix}}", allow_missing = True, out_path = str(config["out"])),
        out_path = str(config["out"])
    resources:
        mem_mb=100_000, # MB
        runtime=2880, # 48 hours
        cpus_per_task=20,
        slurm_extra= '-o $HOME/logs/smk_pfc_%j.out -e $HOME/logs/smk_pfc_%j.out --cpus-per-task=20'
    shell:
        """
        ln -s $HOME/logs/smk_pfc_$SLURM_JOB_ID.out {log}.$SLURM_JOB_ID.out | true
        ln -s $HOME/logs/smk_pfc_$SLURM_JOB_ID.err {log}.$SLURM_JOB_ID.err | true

        if [ ! -s {input.ID} ]
        then
            touch {output}
        else
            {plink2} --bgen {input.bgen} ref-last --sample {input.sample_file} --keep-allele-order --make-pgen --indiv-sort natural  --rm-dup force-first 'list'  --out {params.pgen_out_path} --extract {input.ID} --threads 20 --memory 99000 >> {log}
        fi
        """

rule plink_for_chr_2:
    threads: 1
    input:
        pgen=expand("{out_path}/plink/plink_temp_chr_{{CHROM}}.{{suffix}}.pgen", allow_missing = True, out_path = str(config["out"])),
        altID = expand("{out_path}/rsids/pgs_rsids_alt_chr_{{CHROM}}.{{suffix}}.csv", out_path = str(config["out"])),
        ID = expand("{out_path}/rsids/pgs_rsids_chr_{{CHROM}}.{{suffix}}.csv", out_path = str(config["out"]))
    output:
        expand("{out_path}/plink/plink_temp_chr_{{CHROM}}.{{suffix}}.raw", allow_missing = True, out_path = str(config["out"]))
    log:
        expand("{out_path}/logs/plink/plink_temp_chr_{{CHROM}}.{{suffix}}.raw.log", allow_missing = True, out_path = str(config["out"]))
    resources:
        mem_mb=100_000, # MB
        runtime=2880, # 48 hours
        slurm_extra= '-o $HOME/logs/smk_pfc2_%j.out -e $HOME/logs/smk_pfc2_%j.out --cpus-per-task=20'
    params:
        pgen_in_path=expand("{out_path}/plink/plink_temp_chr_{{CHROM}}.{{suffix}}", allow_missing = True, out_path = str(config["out"])),
        tmp_out_path=expand("{out_path}/plink/plink_temp_chr_{{CHROM}}_{{suffix}}_tmp", allow_missing = True, out_path = str(config["out"])),
        out_path = str(config["out"]),
    shell:
        """
        ln -s $HOME/logs/smk_pfc2_$SLURM_JOB_ID.out {log}.$SLURM_JOB_ID.out | true
        ln -s $HOME/logs/smk_pfc2_$SLURM_JOB_ID.err {log}.$SLURM_JOB_ID.err | true

        # TODO: --rm-dup ref-first flag necessary? yes
        if [ ! -s {input.pgen} ]
        then
            touch {output}
        else
            {plink2} --pfile {params.pgen_in_path} --export-allele {input.altID} --extract {input.ID} --debug --export A --out {params.tmp_out_path}  --memory 99000 --threads 20 >> {log}
            mv {params.tmp_out_path}.raw {output}
        fi
        """
#########################################################################################################################

rule write_split_npy:
    threads: 1
    input:
        var = expand("{out_path}/annotated_df/snpeff_df.{{nfeatures}}.txt", allow_missing = True, out_path = str(config["out"])),
        raw =  expand("{out_path}/plink/plink_temp_chr_9.{{nfeatures}}.raw", allow_missing = True, out_path = str(config["out"])), # 1 is placeholder for all 22
        chrs = expand("{out_path}/plink/plink_temp_chr_{CHROM}.{{nfeatures}}.raw", allow_missing = True, out_path = str(config["out"]), CHROM=config["CHROM"]), # 1 is placeholder for all 22
        pick_list = expand("/sc-projects/sc-proj-dh-ukb-intergenics/analysis/results/genopheno/eid_lists_subset_imputed/{split}.controlled.txt", allow_missing = True, split=['all'])
    output:
        expand("{out_path}/splits/all_chrs_{split}.{{dsname}}-{{nfeatures}}.X.npy", allow_missing = True, out_path = str(config["out"]), split=["all"]),
        expand("{out_path}/splits/all_chrs_{split}.{{dsname}}-{{nfeatures}}.var.csv", allow_missing = True, out_path = str(config["out"]), split=["all"])
    log:
        expand("{out_path}/logs/write_split_npy.{{dsname}}-{{nfeatures}}.log",allow_missing=True, out_path=str(config["out"]))
    params:
        nfeatures = str(config["nfeatures"]),
        out_path = str(config["out"])
    resources:
        mem_mb=100_000, # MB
        runtime=2880, # 48 hours
        cpus_per_task=20,
        slurm_extra= '-o $HOME/logs/smk_npy_%j.out -e $HOME/logs/smk_npy_%j.out --cpus-per-task=20'
    conda:
        'genotype_selection'
    shell:
        """
        ln -s $HOME/logs/smk_npy_$SLURM_JOB_ID.out {log}.$SLURM_JOB_ID.out | true
        ln -s $HOME/logs/smk_npy_$SLURM_JOB_ID.err {log}.$SLURM_JOB_ID.err | true

        python -m scripts.to_npy_split  \
        --i_raw {params.out_path}/plink/plink_temp_chr_{{}}.{wildcards.nfeatures}.raw \
        --i_var {input.var} \
        --pickpath {input.pick_list} \
        --suffix {wildcards.dsname}-{wildcards.nfeatures} \
        --output {params.out_path} \
        --eid_mapping_file /sc-resources/ukb/data/shared/open/fam_files/link_file.tsv \
        --compress \
        --eid_map_to EID.49966 \
        --eid_map_from EID.44448 > {log}
        """

#########################################################################################################################

rule require_nfeatures:
    threads: 1
    input:
        expand("{out_path}/splits/all_chrs_{split}.{dsname}-{nfeatures}.X.npy", allow_missing=True, out_path=str(config["out"]), split=['all'], nfeatures=str(config["nfeatures"]), dsname=str(config["dsname"])),




#########################################################################################################################
#                                                    Evaluation rules                                                   #
#########################################################################################################################




rule calculate_pca_for_train_split:
    input:
        expand("{out_path}/splits/all_chrs_{split}.{dsname}-{nfeatures}.obs.csv", allow_missing=True,  out_path=str(config["out"]), split=['all'], nfeatures=str(config["nfeatures"]), dsname=str(config["dsname"])),
        
    output:
        expand("{out_path}/explained_variance_10k_{split}.{dsname}-{nfeatures}.npy", allow_missing=True, out_path=str(config["out"]), nfeatures=str(config["nfeatures"]),split=['all'], dsname=str(config["dsname"])),

    params:
        out_path = str(config["out"]),
        ds_name = str(config["dsname"]),
        split=['all'],
        nfeatures = str(config["nfeatures"]),

    resources:
        mem_mb=300000,
        runtime=480 
    shell:
        """
        python -m scripts.pca {params.out_path}/splits/all_chrs_{params.split}.{params.ds_name} {params.nfeatures}
        """
