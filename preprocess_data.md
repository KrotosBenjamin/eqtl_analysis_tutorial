
# Table of Contents

-   [Pre-process data and align samples](#org5aa9ab1)
    -   [Sample selection and GCT format](#org0c0cfa8)
        -   [Organize data](#orgf2da91c)
        -   [Python session information](#org41faec2)
    -   [Genotype formatting](#org0f31a4d)
    -   [Expression formatting](#org716e3ce)
    -   [Generate covariates](#org77a8e66)



<a id="org5aa9ab1"></a>

# Pre-process data and align samples

One of the biggest errors I have often run into with using
either fastQTL or tensorQTL is an incorrect order of samples
across expression, genotype, and covariates data. So, this
section focus is getting the input data into a format that
will work with tensorQTL.


<a id="org0c0cfa8"></a>

## Sample selection and GCT format

This set of functions are used to:

1.  select individuals with genotypes
2.  generate a list to map expression to genotypes IDs (RNum to BrNum)
3.  chromosomes to be assessed
4.  convert normalized counts to GCT format

The GCT format is used by the authors of fastQTL and tensorQTL.
It is not necessary as long as the final input for tensorQTL is in
the right format.

Example script is provided: <./scripts/01.prepare_gct.py>.


<a id="orgf2da91c"></a>

### Organize data

    import pandas as pd
    from functools import lru_cache
    def to_gct(filename, df):
        description_df = pd.DataFrame({'Description': df.index.values},index=df.index)
        dfo = pd.concat([description_df, df], axis=1)
        dfo.index.name = 'Names'
        with open(filename, "wt") as out:
    	print("#1.2", file=out)
    	print(df.shape[0], df.shape[1], sep="\t", file=out)
    	dfo.to_csv(out, sep="\t")

    @lru_cache()
    def get_pheno():
        return pd.read_csv("data/caudate_phenotypes.csv", index_col=0)
    
    get_pheno().iloc[0:2, 0:6]

             BrNum    RNum   Region  RIN    Age Sex
    R12864  Br1303  R12864  Caudate  9.6  42.98   F
    R12865  Br1320  R12865  Caudate  9.5  53.12   M

    @lru_cache()
    def get_fam():
        ## Edit for location of genotypes
        fam_file = "input/TOPMed_LIBD_AA_EA.fam"
        return pd.read_csv(fam_file, sep="\t", header=None,
    		       names=["ID","BrNum","V2","V3","V4","V5"])
    
    get_fam().head(2)

                      ID   BrNum  V2  V3  V4  V5
    0  3998646007_R01C01  Br2585   0   0   2  -9
    1  3998646007_R02C01  Br2565   0   0   2  -9

    @lru_cache()
    def load_data():
        pheno_df = get_pheno()
        pheno_df["ids"] = pheno_df.RNum
        pheno_df.set_index("ids", inplace=True)
        norm_df = pd.read_csv("data/caudate.normalized_expression.tsv",
    			  sep="\t", index_col=0)
        samples = list(set(norm_df.columns).intersection(set(pheno_df["RNum"])))
        return pheno_df.loc[samples,:], norm_df.loc[:,samples]
    
    pheno_df, norm_df = load_data()
    print(pheno_df.shape)
    print(norm_df.shape)

    (444, 11)
    (22465, 444)

Now, we'll extract the selected samples.

    def select_idv(pheno_df, norm_df):
        samples = list(set(pheno_df.loc[norm_df.columns,:].BrNum)\
    		   .intersection(set(get_fam().BrNum)))
        new_fam = get_fam()[(get_fam()["BrNum"].isin(samples))]\
    	.drop_duplicates(subset="BrNum")
        new_fam.to_csv("data/keepFam.txt", sep='\t', index=False, header=False)
        return pheno_df.loc[:, ["RNum", "BrNum"]]\
    		   .reset_index().set_index("BrNum")\
    		   .loc[new_fam.BrNum].reset_index().set_index("ids")
    
    
    new_pheno = select_idv(pheno_df, norm_df)
    new_pheno.head(2)

             BrNum    RNum
    ids                   
    R12995  Br2585  R12995
    R13019  Br5073  R13019

    to_gct("data/norm.gct", norm_df.loc[:,new_pheno.index])
    new_pheno.loc[:, ["RNum", "BrNum"]]\
    	 .to_csv("data/sample_id_to_brnum.tsv", sep="\t", index=False)
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('data/vcf_chr_list.txt', header=False, index=None)


<a id="org41faec2"></a>

### Python session information

    import session_info
    session_info.show()

    -----
    pandas              1.5.3
    session_info        1.0.0
    -----
    Python 3.10.10 | packaged by conda-forge | (main, Mar 24 2023, 20:08:06) [GCC 11.3.0]
    Linux-3.10.0-1160.el7.x86_64-x86_64-with-glibc2.17
    -----
    Session information updated at 2023-09-28 11:34


<a id="org0f31a4d"></a>

## Genotype formatting

Now that we have samples selected and mapping files, we can format our
genotype data. Note: this will be placed in a protected location.

I'll be working on JHPCE for this. This should also order the samples.

    module load plink/2.0
    plink2 --bfile input/TOPMed_LIBD_AA_EA \
           --keep data/keepFam.txt --make-bed \
           --out input/protected_data/genotypes

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-right" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">PLINK</td>
<td class="org-left">v2.00a3LM</td>
<td class="org-left">64-bit</td>
<td class="org-left">Intel</td>
<td class="org-left">(17</td>
<td class="org-left">Dec</td>
<td class="org-right">2021)</td>
<td class="org-left">www.cog-genomics.org/plink/2.0/</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">(C)</td>
<td class="org-left">2005-2021</td>
<td class="org-left">Shaun</td>
<td class="org-left">Purcell,</td>
<td class="org-left">Christopher</td>
<td class="org-left">Chang</td>
<td class="org-right">GNU</td>
<td class="org-left">General</td>
<td class="org-left">Public</td>
<td class="org-left">License</td>
<td class="org-left">v3</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Logging</td>
<td class="org-left">to</td>
<td class="org-left">input/protected_data/genotypes.log.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Options</td>
<td class="org-left">in</td>
<td class="org-left">effect:</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;bfile</td>
<td class="org-left">input/TOPMed_LIBD_AA_EA</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;keep</td>
<td class="org-left">data/keepFam.txt</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;make-bed</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;out</td>
<td class="org-left">input/protected_data/genotypes</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Start</td>
<td class="org-left">time:</td>
<td class="org-left">Thu</td>
<td class="org-left">Sep</td>
<td class="org-left">28</td>
<td class="org-left">11:29:56</td>
<td class="org-right">2023</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">499853</td>
<td class="org-left">MiB</td>
<td class="org-left">RAM</td>
<td class="org-left">detected;</td>
<td class="org-left">reserving</td>
<td class="org-left">249926</td>
<td class="org-right">MiB</td>
<td class="org-left">for</td>
<td class="org-left">main</td>
<td class="org-left">workspace.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Allocated</td>
<td class="org-left">7915</td>
<td class="org-left">MiB</td>
<td class="org-left">successfully,</td>
<td class="org-left">after</td>
<td class="org-left">larger</td>
<td class="org-right">attempt(s)</td>
<td class="org-left">failed.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Using</td>
<td class="org-left">up</td>
<td class="org-left">to</td>
<td class="org-left">64</td>
<td class="org-left">threads</td>
<td class="org-left">(change</td>
<td class="org-right">this</td>
<td class="org-left">with</td>
<td class="org-left">&#x2013;threads).</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">1938</td>
<td class="org-left">samples</td>
<td class="org-left">(725</td>
<td class="org-left">females,</td>
<td class="org-left">1209</td>
<td class="org-left">males,</td>
<td class="org-right">4</td>
<td class="org-left">ambiguous;</td>
<td class="org-left">1938</td>
<td class="org-left">founders)</td>
<td class="org-left">loaded</td>
<td class="org-left">from</td>
</tr>


<tr>
<td class="org-left">input/TOPMed_LIBD_AA_EA.fam.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">7678274</td>
<td class="org-left">variants</td>
<td class="org-left">loaded</td>
<td class="org-left">from</td>
<td class="org-left">input/TOPMed_LIBD_AA_EA.bim.</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Note:</td>
<td class="org-left">No</td>
<td class="org-left">phenotype</td>
<td class="org-left">data</td>
<td class="org-left">present.</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;keep:</td>
<td class="org-left">435</td>
<td class="org-left">samples</td>
<td class="org-left">remaining.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">435</td>
<td class="org-left">samples</td>
<td class="org-left">(141</td>
<td class="org-left">females,</td>
<td class="org-left">294</td>
<td class="org-left">males;</td>
<td class="org-right">435</td>
<td class="org-left">founders)</td>
<td class="org-left">remaining</td>
<td class="org-left">after</td>
<td class="org-left">main</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">filters.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Writing</td>
<td class="org-left">input/protected_data/genotypes.fam</td>
<td class="org-left">&#x2026;</td>
<td class="org-left">done.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Writing</td>
<td class="org-left">input/protected_data/genotypes.bim</td>
<td class="org-left">&#x2026;</td>
<td class="org-left">done.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Writing</td>
<td class="org-left">input/protected_data/genotypes.bed</td>
<td class="org-left">&#x2026;</td>
<td class="org-left">0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-right">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">End</td>
<td class="org-left">time:</td>
<td class="org-left">Thu</td>
<td class="org-left">Sep</td>
<td class="org-left">28</td>
<td class="org-left">11:30:31</td>
<td class="org-right">2023</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>
</tbody>
</table>


<a id="org716e3ce"></a>

## Expression formatting

For expression formatting, we need to:

1.  convert to BED format with gene information (i.e., chromosome, start, end)
2.  replace expression ids with genotype ids
3.  compress and index expression file

For this, we will used an adapted version of [eqtl\_prepare\_expression.py](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/eqtl_prepare_expression.py) from
the fastQTL/tensorQTL authors. Details on how they used this in the
GTEx QTL workflow can be found [here](https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md).

The modified helper script takes the following input:

1.  normalized data: GCT format
2.  BED file with gene annotation
3.  sample ID mapping file
4.  chromosomes to analyze

    python3 ./scripts/02.prepare_expression.py --help

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />

<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">usage:</td>
<td class="org-left">02.prepare_expression.py</td>
<td class="org-left">[-h]</td>
<td class="org-left">[-o</td>
<td class="org-left">OUTPUT_DIR]</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">[&#x2013;sample_id_list</td>
<td class="org-left">SAMPLE_ID_LIST]</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">[&#x2013;feature</td>
<td class="org-left">FEATURE]</td>
<td class="org-left">[&#x2013;bed_file</td>
<td class="org-left">BED_FILE]</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">norm_gct</td>
<td class="org-left">sample_participant_lookup</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">vcf_chr_list</td>
<td class="org-left">prefix</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Generate</td>
<td class="org-left">normalized</td>
<td class="org-left">expression</td>
<td class="org-left">BED</td>
<td class="org-left">files</td>
<td class="org-left">for</td>
<td class="org-left">eQTL</td>
<td class="org-left">analyses</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">positional</td>
<td class="org-left">arguments:</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">norm_gct</td>
<td class="org-left">GCT</td>
<td class="org-left">file</td>
<td class="org-left">with</td>
<td class="org-left">normalized</td>
<td class="org-left">expression</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">sample_participant_lookup</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Lookup</td>
<td class="org-left">table</td>
<td class="org-left">linking</td>
<td class="org-left">samples</td>
<td class="org-left">to</td>
<td class="org-left">participants</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">vcf_chr_list</td>
<td class="org-left">List</td>
<td class="org-left">of</td>
<td class="org-left">chromosomes</td>
<td class="org-left">in</td>
<td class="org-left">VCF</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">prefix</td>
<td class="org-left">Prefix</td>
<td class="org-left">for</td>
<td class="org-left">output</td>
<td class="org-left">file</td>
<td class="org-left">names</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">options:</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">-h,</td>
<td class="org-left">&#x2013;help</td>
<td class="org-left">show</td>
<td class="org-left">this</td>
<td class="org-left">help</td>
<td class="org-left">message</td>
<td class="org-left">and</td>
<td class="org-left">exit</td>
</tr>


<tr>
<td class="org-left">-o</td>
<td class="org-left">OUTPUT_DIR,</td>
<td class="org-left">&#x2013;output_dir</td>
<td class="org-left">OUTPUT_DIR</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">Output</td>
<td class="org-left">directory</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;sample_id_list</td>
<td class="org-left">SAMPLE_ID_LIST</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">File</td>
<td class="org-left">listing</td>
<td class="org-left">sample</td>
<td class="org-left">IDs</td>
<td class="org-left">to</td>
<td class="org-left">include</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;feature</td>
<td class="org-left">FEATURE</td>
<td class="org-left">gene,</td>
<td class="org-left">transcript</td>
<td class="org-left">or</td>
<td class="org-left">exon</td>
<td class="org-left">&#xa0;</td>
<td class="org-left">&#xa0;</td>
</tr>


<tr>
<td class="org-left">&#x2013;bed_file</td>
<td class="org-left">BED_FILE</td>
<td class="org-left">this</td>
<td class="org-left">is</td>
<td class="org-left">the</td>
<td class="org-left">bed</td>
<td class="org-left">file</td>
<td class="org-left">annotation</td>
</tr>
</tbody>
</table>

    module load htslib
    module load samtools
    
    BED="./data/gene.bed"
    python3 ./scripts/02.prepare_expression.py \
    	--feature gene --bed_file $BED -o data/ \
    	./data/norm.gct ./data/sample_id_to_brnum.tsv \
    	./data/vcf_chr_list.txt genes

    Loading expression data
    Map data
      * 22465 genes.
    bed_template_df.shape (22465, 4)
      * 22398 genes remain after removing contigs absent from VCF.
    Writing BED file

    ls data/genes*

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />
</colgroup>
<tbody>
<tr>
<td class="org-left">data/genes.expression.bed.gz</td>
</tr>


<tr>
<td class="org-left">data/genes.expression.bed.gz.tbi</td>
</tr>
</tbody>
</table>


<a id="org77a8e66"></a>

## Generate covariates

In concurrent with expression and genotype formatting, we
also need to generate covariates for our gene expression
data.

I've provided a helper [script](./scripts/03.generate_covs.R) that will load counts assuming
R object is a `RangedSummarizedExperiment`.

Below, I walk through the script with edits for just
gene level analysis and the caudate data we've already
loaded into our R session.

    suppressPackageStartupMessages({
        library(SummarizedExperiment)
    })
    
    getRPKM <- function(rse, length_var = "bp_length", mapped_var = NULL) {
        mapped <- if (!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)
        bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)
        len <- if (!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))
        wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
        return(assays(rse)$counts / (wid / 1000) / (bg / 1e6))
    }

Normalize data for PCA analysis.

    keepInd   <- which((colData(rse_gene)$Age > 13) &
    		   (colData(rse_gene)$Race %in% c("AA", "CAUC")))
    rse_gene  <- rse_gene[,keepInd]
    gene_rpkm <- getRPKM(rse_gene, "Length")
    rse_gene  <- rse_gene[rowMeans(gene_rpkm) > 0.2,]; rm(gene_rpkm)
    norm_df   <- getRPKM(rse_gene, 'Length')
    colnames(norm_df) <- gsub("\\_.*", "", colnames(norm_df))
    norm_df[1:3, 1:5]

                         R12864   R12865   R12866    R12867   R12868
    ENSG00000227232.5 2.1068042 2.135875 2.467500 2.3927082 3.263655
    ENSG00000278267.1 4.6081368 3.160039 8.591733 3.2047741 4.381159
    ENSG00000269981.1 0.6436247 0.864719 1.452124 0.8952303 1.049010

Load  MDS data.

    suppressPackageStartupMessages({library(dplyr)})
    gfile  <- "input/TOPMed_LIBD_AA_EA.mds"
    mds    <- data.table::fread(gfile) |>
      rename_at(.vars = vars(starts_with("C")),
    	    function(x){sub("C", "snpPC", x)}) |>
      mutate_if(is.character, as.factor)
    mds[1:2, 1:5]

                     FID    IID SOL     snpPC1     snpPC2
    1: 3998646007_R01C01 Br2585   0 -0.1168750 0.00105119
    2: 3998646007_R02C01 Br2565   0  0.0389051 0.00175104

Merge with phenotype information.

    new_pd <- colData(rse_gene) |> as.data.frame() |>
    	 inner_join(mds, by=c("BrNum"="IID"), multiple="all") |>
    	 distinct(RNum, .keep_all = TRUE) |> mutate(ids=RNum) |>
    	 tibble::column_to_rownames("ids")
    new_pd[1:2, 1:5]

            BrNum   RNum  Region RIN   Age
    R12864 Br1303 R12864 Caudate 9.6 42.98
    R12865 Br1320 R12865 Caudate 9.5 53.12

Generate model and do PCA.

    mod <- model.matrix(~Sex + Dx + Age + snpPC1 + snpPC2 + snpPC3,
    		    data=new_pd)
    pca_df <- prcomp(t(log2(norm_df[, new_pd$RNum]+1)))
    k      <- sva::num.sv(log2(norm_df[, new_pd$RNum]+1), mod)
    pcs    <- pca_df$x[, 1:k]
    dim(pcs)

    [1] 435  26

    pcs[1:2, 1:5]

                  PC1       PC2       PC3         PC4       PC5
    R12864 -25.639291  14.28025 -10.27202  0.05150028 -3.808366
    R12865  -2.663915 -24.74939  18.96207 -7.36002876 -9.967997

Subset samples and generate covariates.

    sample_df <- data.table::fread("./data/sample_id_to_brnum.tsv")
    head(sample_df, 2)

         RNum  BrNum
    1: R12995 Br2585
    2: R13019 Br5073

    head( mod, 2)

           (Intercept) SexM DxControl DxSchizo   Age     snpPC1      snpPC2
    R12864           1    0         0        1 42.98 -0.0940914 0.000457512
    R12865           1    1         0        1 53.12 -0.1119230 0.001828400
               snpPC3
    R12864 -0.0010998
    R12865  0.0018495

    covs <- cbind(mod[,c(-1)], pcs) |> as.data.frame() |>
        tibble::rownames_to_column("RNum") |>
        inner_join(sample_df, by=c("RNum")) |> select(-"RNum") |>
        rename("ID"="BrNum") |> tibble::column_to_rownames("ID") |>
        t() |> as.data.frame() |> tibble::rownames_to_column("ID")
    covs <- covs[, c("ID", sample_df$BrNum)]
    dim(covs)

    [1]  33 436

    covs[1:5, 1:5]

             ID    Br2585     Br5073     Br5347    Br5179
    1      SexM  0.000000  1.0000000  0.0000000  1.000000
    2 DxControl  0.000000  0.0000000  1.0000000  0.000000
    3  DxSchizo  1.000000  1.0000000  0.0000000  1.000000
    4       Age 33.350000 62.6100000 80.5400000 53.040000
    5    snpPC1 -0.116875 -0.0734824  0.0363273 -0.083326

Now, we can write the covarites to a file.

    data.table::fwrite(covs, "data/genes.combined_covariates.txt",
    		   sep='\t')

Now, we are ready to run tensorQTL!

