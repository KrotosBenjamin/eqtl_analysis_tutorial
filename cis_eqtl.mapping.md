
# Table of Contents

-   [cis-eQTL analysis with tensorQTL](#org2be8d66)
    -   [Nominal cis-eQTL analysis](#orgfa5fd99)
        -   [Genotypes](#org19a1028)
        -   [Covariates](#orgb114999)
        -   [Phenotype (i.e., expression)](#orgdf90ecc)
        -   [Main mapping](#orgf37992b)
    -   [Permutation analysis](#org7ecf3d3)
    -   [Post hoc](#orge4218d4)
-   [Conditional analysis and fine mapping](#org98f70fd)
    -   [Conditional analysis](#orga7a7d3c)



<a id="org2be8d66"></a>

# cis-eQTL analysis with tensorQTL

First thing that needs to be done is installing `tensorQTL`.

I am currently using version 1.0.7 with PyTorch installed
using cuda. Please make sure to install PyTorch correctly
otherwise you will not be able to leverage the GPUs on
your machine.

I provide a several helper scripts so that this can
be run on the command line. For now, let's review
the steps in more detail.

I'll setup the CUDA now.

    nvidia-smi

    Fri Sep 29 06:35:12 2023       
    +-----------------------------------------------------------------------------+
    | NVIDIA-SMI 470.57.02    Driver Version: 470.57.02    CUDA Version: 11.4     |
    |-------------------------------+----------------------+----------------------+
    | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
    | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
    |                               |                      |               MIG M. |
    |===============================+======================+======================|
    |   0  Tesla V100S-PCI...  Off  | 00000000:3B:00.0 Off |                    0 |
    | N/A   41C    P0    26W / 250W |      4MiB / 32510MiB |      0%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+
    |   1  Tesla V100S-PCI...  Off  | 00000000:5E:00.0 Off |                    0 |
    | N/A   32C    P0    25W / 250W |      4MiB / 32510MiB |      0%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+
    |   2  Tesla V100S-PCI...  Off  | 00000000:86:00.0 Off |                    0 |
    | N/A   33C    P0    24W / 250W |      4MiB / 32510MiB |      0%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+
    |   3  Tesla V100S-PCI...  Off  | 00000000:AF:00.0 Off |                    0 |
    | N/A   32C    P0    24W / 250W |      4MiB / 32510MiB |      0%      Default |
    |                               |                      |                  N/A |
    +-------------------------------+----------------------+----------------------+
                                                                                   
    +-----------------------------------------------------------------------------+
    | Processes:                                                                  |
    |  GPU   GI   CI        PID   Type   Process name                  GPU Memory |
    |        ID   ID                                                   Usage      |
    |=============================================================================|
    |  No running processes found                                                 |
    +-----------------------------------------------------------------------------+

    export CUDA_VISIBLE_DEVICES=0


<a id="orgfa5fd99"></a>

## Nominal cis-eQTL analysis

The main function for `tensorQTL` is from the `cis`
module. But we also need to read in the expression
data (i.e., `read_phenotype_bed`) and genotype
data (i.e., `genotypeio`).

We also need to make sure we have R loaded as the
permutation function uses R to calculate q-values.

I'll need to load the version of python and R that
work with `rpy2` for my configuration on JHPCE.

    import pandas as pd
    from functools import lru_cache
    from tensorqtl import cis, read_phenotype_bed, genotypeio

First lets load the covariates, expression, and genotype data.


<a id="org19a1028"></a>

### Genotypes

There are several genotype reading (IO) function available,
including the ablilty to load the data from a VCF file.
For PLINK (BED/BIM/FAM) format, the `PlinkReader` function
will assume the genotype names are in the `IID` column.

If this is not true, you'll need to redo your genotype
processing.

    @lru_cache()
    def get_genotypes():
        plink_prefix_path = "input/protected_data/genotypes"
        pr = genotypeio.PlinkReader(plink_prefix_path)
        variant_df = pr.bim.set_index("snp")[["chrom", "pos"]]
        variant_df.loc[:, "chrom"] = "chr" + variant_df.chrom
        return pr.load_genotypes(), variant_df
    
    
    genotype_df, variant_df = get_genotypes()

    

    genotype_df.shape

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<tbody>
<tr>
<td class="org-right">7678274</td>
<td class="org-right">435</td>
</tr>
</tbody>
</table>

    variant_df.head(2)

                                 chrom      pos
    snp                                        
    chr1_833068_A_G_rs12562034    chr1   833068
    chr1_1057648_A_G_rs144629300  chr1  1057648


<a id="orgb114999"></a>

### Covariates

    @lru_cache()
    def get_covars():
        covar_file = "data/genes.combined_covariates.txt"
        return pd.read_csv(covar_file, sep='\t', index_col=0).T
    
    get_covars().shape  

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<tbody>
<tr>
<td class="org-right">435</td>
<td class="org-right">33</td>
</tr>
</tbody>
</table>


<a id="orgdf90ecc"></a>

### Phenotype (i.e., expression)

    @lru_cache()
    def get_phenotype():
        expr_bed = "data/genes.expression.bed.gz"
        return read_phenotype_bed(expr_bed)
    
    phenotype_df, phenotype_pos_df = get_phenotype()
    phenotype_df.shape

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<tbody>
<tr>
<td class="org-right">22398</td>
<td class="org-right">435</td>
</tr>
</tbody>
</table>


<a id="orgf37992b"></a>

### Main mapping

The main function for nominal analysis is `cis.map_nominal`.

    help(cis.map_nominal)

    Help on function map_nominal in module cis:
    
    map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=None, maf_threshold=0, interaction_df=None, maf_threshold_interaction=0.05, group_s=None, window=1000000, run_eigenmt=False, output_dir='.', write_top=True, write_stats=True, logger=None, verbose=True)
        cis-QTL mapping: nominal associations for all variant-phenotype pairs
        
        Association results for each chromosome are written to parquet files
        in the format <output_dir>/<prefix>.cis_qtl_pairs.<chr>.parquet
        
        If interaction_df is provided, the top association per phenotype is
        written to <output_dir>/<prefix>.cis_qtl_top_assoc.txt.gz unless
        write_top is set to False, in which case it is returned as a DataFrame

    mkdir output

We'll run the analysis using MAF threshold of 0.01 and window of 0.5Mbp.

    prefix = "caudate"
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
    		prefix, covariates_df=get_covars(),
    		maf_threshold=0.01, window=500000,
    		output_dir="output/")

      * 7678274 variants
      * applying in-sample 0.01 MAF filter
    
        ** dropping 169 phenotypes without variants in cis-window
      * Computing associations
        Mapping chromosome chr1
    
        * writing output
        Mapping chromosome chr2
    
        * writing output
        Mapping chromosome chr3
    
        * writing output
        Mapping chromosome chr4
    
        * writing output
        Mapping chromosome chr5
    
        * writing output
        Mapping chromosome chr6
    
        * writing output
        Mapping chromosome chr7
    
        * writing output
        Mapping chromosome chr8
    
        * writing output
        Mapping chromosome chr9
    
        * writing output
        Mapping chromosome chr10
    
        * writing output
        Mapping chromosome chr11
    
        * writing output
        Mapping chromosome chr12
    
        * writing output
        Mapping chromosome chr13
    
        * writing output
        Mapping chromosome chr14
    
        * writing output
        Mapping chromosome chr15
    
        * writing output
        Mapping chromosome chr16
    
        * writing output
        Mapping chromosome chr17
    
        * writing output
        Mapping chromosome chr18
    
        * writing output
        Mapping chromosome chr19
    
        * writing output
        Mapping chromosome chr20
    
        * writing output
        Mapping chromosome chr21
    
        * writing output
        Mapping chromosome chr22
    
        * writing output
        Mapping chromosome chrX
    
        time elapsed: 3.23 min
        * writing output
    done.

The results are saved as compressed, parquet files. And the total
time was less than 5 minutes.

    ls output/

    caudate.cis_qtl_pairs.chr10.parquet
    caudate.cis_qtl_pairs.chr11.parquet
    caudate.cis_qtl_pairs.chr12.parquet
    caudate.cis_qtl_pairs.chr13.parquet
    caudate.cis_qtl_pairs.chr14.parquet
    caudate.cis_qtl_pairs.chr15.parquet
    caudate.cis_qtl_pairs.chr16.parquet
    caudate.cis_qtl_pairs.chr17.parquet
    caudate.cis_qtl_pairs.chr18.parquet
    caudate.cis_qtl_pairs.chr19.parquet
    caudate.cis_qtl_pairs.chr1.parquet
    caudate.cis_qtl_pairs.chr20.parquet
    caudate.cis_qtl_pairs.chr21.parquet
    caudate.cis_qtl_pairs.chr22.parquet
    caudate.cis_qtl_pairs.chr2.parquet
    caudate.cis_qtl_pairs.chr3.parquet
    caudate.cis_qtl_pairs.chr4.parquet
    caudate.cis_qtl_pairs.chr5.parquet
    caudate.cis_qtl_pairs.chr6.parquet
    caudate.cis_qtl_pairs.chr7.parquet
    caudate.cis_qtl_pairs.chr8.parquet
    caudate.cis_qtl_pairs.chr9.parquet
    caudate.cis_qtl_pairs.chrX.parquet


<a id="org7ecf3d3"></a>

## Permutation analysis

Next, we want to use permutation analysis to identify
significant associations. This will take more time, but
it is still much faster than other alternatives.

    help(cis.map_cis)

    Help on function map_cis in module cis:
    
    map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=None, group_s=None, maf_threshold=0, beta_approx=True, nperm=10000, window=1000000, random_tiebreak=False, logger=None, seed=None, verbose=True, warn_monomorphic=True)
        Run cis-QTL mapping

    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df,
    		     phenotype_pos_df, covariates_df=get_covars(),
    		     maf_threshold=0.01, window=500000, seed=13131313)

    cis-QTL mapping: empirical p-values for phenotypes
      * 435 samples
      * 22398 phenotypes
      * 33 covariates
      * 7678274 variants
      * applying in-sample 0.01 MAF filter
      * using seed 13131313
    
        ** dropping 169 phenotypes without variants in cis-window
      * computing permutations
    
      return 2*stats.t.cdf(-np.abs(np.sqrt(tstat2)), dof)
    WARNING: scipy.optimize.newton failed to converge (running scipy.optimize.minimize)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    