
# Table of Contents

-   [eQTL background](#org5e40ba6)
    -   [Why eQTL analysis](#org2caf46b)
-   [Prepare data](#orge5331b5)
    -   [BrainSEQ data](#orgbd74bc7)
    -   [Study design](#org4233d62)
    -   [Phenotype data](#orgc5a90db)
    -   [Normalized counts](#org5ad67ee)
        -   [Export gene annotation](#org86016fa)
    -   [Genotypes](#org03600d7)
        -   [Population structure](#orgf142af6)
    -   [R session information](#org0f31c6c)



<a id="org5e40ba6"></a>

# eQTL background

Quantitative trait loci (QTL) analysis relays on genetic variant
data (i.e., SNP) and another kind of molecular trait (i.e.,
expression, protein, or DNA methylation). The idea is to test
if a specific molecular trait changes based on allele and
allele dosage.

In general, this pipeline can be used for any kind of quantitative
data with some adaptation. However, I will be focusing on
expression QTL (eQTL) for the remainder of this tutorial.


<a id="org2caf46b"></a>

## Why eQTL analysis

![img](./img/eqtl_summary.v2.png)

eQTL are genetic variants that explain variation in expression levels.

-   Translate GWAS to biology
-   Useful for gene prioritization
-   Detection of novel gene-trait associations
-   Inferring direction of association

Therefore, to test for these associations, we need two main information
from the same samples: 1) expression and 2) genotypes.


<a id="orge5331b5"></a>

# Prepare data

For this analysis, we will be using the BrainSEQ data
([PMID: 26687217](https://www.ncbi.nlm.nih.gov/pubmed/26687217)), which saves processed data in R objects.
One important note is that publiclly available human
genotype data is always controlled access due to privacy
issues. This takes time to get access too, so please
address this issue first.


<a id="orgbd74bc7"></a>

## BrainSEQ data

Processed data are RangedSummarizedExperiment R Objects
(hg38; GENCODE v25).

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">Brain region</th>
<th scope="col" class="org-left">RNA type</th>
<th scope="col" class="org-left">Data</th>
<th scope="col" class="org-left">Reference</th>
</tr>
</thead>

<tbody>
<tr>
<td class="org-left">DLPFC</td>
<td class="org-left">Poly-A</td>
<td class="org-left"><a href="https://s3.us-east-2.amazonaws.com/jaffe-nat-neuro-2018/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda">Gene counts</a></td>
<td class="org-left"><a href="https://pubmed.ncbi.nlm.nih.gov/30050107/">PMID: 30050107</a></td>
</tr>
</tbody>

<tbody>
<tr>
<td class="org-left">DLPFC/Hippocampus</td>
<td class="org-left">total RNA</td>
<td class="org-left"><a href="https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_gene_unfiltered.Rdata">Gene counts</a></td>
<td class="org-left"><a href="https://pubmed.ncbi.nlm.nih.gov/31174959/">PMID: 31174959</a></td>
</tr>
</tbody>

<tbody>
<tr>
<td class="org-left">Caudate nucleus</td>
<td class="org-left">total RNA</td>
<td class="org-left"><a href="https://caudate-eqtl.s3.us-west-2.amazonaws.com/caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda">Gene counts</a></td>
<td class="org-left"><a href="https://pubmed.ncbi.nlm.nih.gov/36319771/">PMID: 36319771</a></td>
</tr>
</tbody>
</table>

These object have counts, phenotype information, and
sequencing information that can be used for study
design.

Genotype data is controlled access.


<a id="org4233d62"></a>

## Study design

Before starting eQTL, it is important to have a goal in
mind. The most common study design is cis-eQTL analysis.
This is assessing variants near (within 1M bp) the transcriptional
start site (TSS). In addition to this, we could also test
for interaction models (genetic variant X phenotype) or
trans-eQTL models. For this tutorial, we will do the more
common cis-eQTL approach. As such, we will use as many
samples as possible to increase our power of detection.

Note: The BrainSEQ data set has neurotypical controls and
individuals with neuropsychiatric disorders (i.e.,
schizophrenia and bipolar disorder). However, diagnosis
status or antipsychotic status does not effect the
eQTL analysis compared with neurotypical controls
([PMID: 36319771](https://pubmed.ncbi.nlm.nih.gov/36319771/)).


<a id="orgc5a90db"></a>

## Phenotype data

Now lets extract phenotype data for the caudate
nucleus.

Let's check to make sure we're in the right place.

    pwd

Now, we'll download the gene-level data from the eQTL data
browser website.

With the data downloaded, we can extract the phenotype information.
First, I'll load in a helpful function from [jaffelab](https://github.com/LieberInstitute/jaffelab).

Using SummarizedExperiment, I can select just the variables
we want to keep.
Now, let's take a look at it.

As we want to use as many samples as possible, we will only do some
basic filtering for our study design:

1.  Including only individual age > 13, and
2.  Limit to self-identified Black and White Americans

Now, we'll save this as a text file to work with python.


<a id="org5ad67ee"></a>

## Normalized counts

We next need normalized counts data. The authors of
tensorQTL and fastQTL transform their counts data using
`edgeR` TMM method.

    x <- edgeR::calcNormFactors(x, method="TMM")

They used a helper set of functions to convert this R
function into python. However, since we are already
importing data in R, we can skip the steps of converting
counts and a set of normalized expression and applying
the helper function to transform it to normalized expression
with python.

Now, we can save normalized counts data.


<a id="org86016fa"></a>

### Export gene annotation


<a id="org03600d7"></a>

## Genotypes

We have our genotype data as both VCF and PLINK format.
For this tutorial, I will assume the genotypes are already
quality controlled and in PLINK format (BED/FAM/BIM).


<a id="orgf142af6"></a>

### Population structure

In addition to having genotypes, we also need information on
population structure. To generate this data, we'll use PLINK
to generate MDS data from pruned data.

    echo "**** Make temporary directory ***"
    mkdir -p tmp

    module load plink/2.0
    
    echo "**** Prune genotypes ****"
    plink2 --bfile input/TOPMed_LIBD_AA_EA \
           --indep-pairwise 500kb 0.5 \
           --out tmp/genotypes

    echo "**** Filtered genotypes ****"
    plink2 --bfile input/TOPMed_LIBD_AA_EA \
           --extract tmp/genotypes.prune.in --make-bed \
           --maf 0.05 --out tmp/TOPMed_LIBD_AA_EA

Note: This step is memory intensive.

    module load plink/1.90b6.6
    
    plink --bfile tmp/TOPMed_LIBD_AA_EA --cluster \
          --mds-plot 10 --out input/TOPMed_LIBD_AA_EA

    rm -rf tmp/


<a id="org0f31c6c"></a>

## R session information

