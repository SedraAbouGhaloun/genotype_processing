# Define config
configfile: "config.yaml"

# Input dataset path from the config
dataset_path = config["dataset_path"]

name = dataset_path.split("/")[-1]
nfeatures = "100k"

rule all:
    input:
        f"/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{name}/ld/r2/afr/chr_10_0.05.ld",
        f"{dataset_path}/explained_variance_10k.npy" 
########################################################################################################################

rule count_genes:
    output:
        unique_genes=f"{dataset_path}/{name}_unique_genes.txt"
    params:
        name=name
    log:
        f"{dataset_path}/logs/{name}_unique_genes.log"
    shell:
        """
        python -m umbrella.pre.gwas.eval_scripts.count_genes {params.name} 
        """


########################################################################################################################

rule ld_preprocessing:
    input:
        unique_genes=rules.count_genes.output.unique_genes
    output:
        processed_data=f"/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{name}/ld/{name}.csv"
    params:
        name=name
    log:
        f"/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{name}/logs/{name}_ld_preprocessing.log"
    shell:
        """
        mkdir -p /sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{params.name}/ld
        python  -m umbrella.pre.gwas.eval_scripts.ld_preprocessing {params.name}
      """

rule run_ld_workflow:
    input:
        processed_data=rules.ld_preprocessing.output.processed_data
    output:
        processed_data=f"/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{name}/ld/r2/afr/chr_10_0.05.ld"
    params:
        name=name

    shell:
        """
        snakemake --cluster "sbatch --mem=120G --time 48:00:00" \
                  --jobs 100 \
                  --latency-wait 60 \
                  -s /home/aboughas/code/umbrella/umbrella/pre/gwas/ld.smk \
                  --configfile /home/aboughas/code/umbrella/umbrella/pre/gwas/ld_config.yaml \
                  --config out=/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{params.name}/ld \
                           snps_file=/sc-scratch/sc-scratch-dh-ukb-intergenics/geno_datasets/{params.name}/ld/{params.name}.csv \
                           ancestry_path=/sc-projects/sc-proj-dh-ukb-intergenics/analysis/development/aboughas/code/ld \
                           ld_window_kb=500 \
                           r2=0.05 \
                  --reason \
                  --quiet \
                  -p 
        """


########################################################################################################################

rule calculate_pca_for_train_split:
    input:
        f"{dataset_path}/splits/all_chrs_all.{name}-{nfeatures}.X.npy"
        
    output:
        f"{dataset_path}/explained_variance_10k.npy"

    params:
        dataset_path = dataset_path,
        name = name,
        nfeatures = nfeatures

    resources:
        mem_mb=300000,
        runtime=480,
        slurm_extra= '-o $HOME/logs/smk_snf_%j.out -e $HOME/logs/smk_snf_%j.out --cpus-per-task=20'

    shell:
        """
        python -m umbrella.pre.gwas.scripts.pca {params.dataset_path}/splits/all_chrs_all.{params.name} {params.nfeatures}
        """

########################################################################################################################



