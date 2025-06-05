# Class 11: AlphaFold
Linh Dang (PID:A16897764)

- [Custom analysis of resulting
  models](#custom-analysis-of-resulting-models)
- [Predicted Alignment Error for
  domains](#predicted-alignment-error-for-domains)
- [Score Residue Conservation from alignment
  file](#score-residue-conservation-from-alignment-file)

# Custom analysis of resulting models

Here we read the results from AlphaFold and try to interpret all the
models and quality score metrics:

``` r
library(bio3d)

pth <- "dimer_23119/"
pdb.files <- list.files(path=pth, full.names = TRUE, pattern = ".pdb")
```

Align and superpose all these models

``` r
file.exists(pdb.files)
```

    [1] TRUE TRUE TRUE TRUE TRUE

``` r
pdbs <- pdbaln(pdb.files, fit = TRUE, exefile="msa")
```

    Reading PDB files:
    dimer_23119//dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb
    dimer_23119//dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000.pdb
    dimer_23119//dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000.pdb
    dimer_23119//dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000.pdb
    dimer_23119//dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000.pdb
    .....

    Extracting sequences

    pdb/seq: 1   name: dimer_23119//dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb 
    pdb/seq: 2   name: dimer_23119//dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000.pdb 
    pdb/seq: 3   name: dimer_23119//dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000.pdb 
    pdb/seq: 4   name: dimer_23119//dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000.pdb 
    pdb/seq: 5   name: dimer_23119//dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000.pdb 

``` r
library(bio3dview)

#view.pdbs(pdbs)
```

``` r
# Read a reference PDB structure
pdb <- read.pdb("1hsg")
```

      Note: Accessing on-line PDB file

``` r
plot(pdbs$b[1,], typ = "l", ylim=c(0,100), ylab="pLDDT score")
lines(pdbs$b[2,], typ = "l", col="blue")
lines(pdbs$b[3,], typ = "l", col="green")
lines(pdbs$b[4,], typ = "l", col="orange")
lines(pdbs$b[5,], typ = "l", col="red")
```

![](class11_alphafold_files/figure-commonmark/unnamed-chunk-6-1.png)

# Predicted Alignment Error for domains

``` r
rd <- rmsd(pdbs)
```

    Warning in rmsd(pdbs): No indices provided, using the 198 non NA positions

``` r
rd
```

                                                                           dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000
    dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000                                                                  0.000
    dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000                                                                  0.157
    dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000                                                                  0.367
    dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000                                                                  0.308
    dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000                                                                 13.306
                                                                           dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000
    dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000                                                                  0.157
    dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000                                                                  0.000
    dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000                                                                  0.383
    dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000                                                                  0.318
    dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000                                                                 13.284
                                                                           dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000
    dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000                                                                  0.367
    dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000                                                                  0.383
    dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000                                                                  0.000
    dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000                                                                  0.483
    dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000                                                                 13.406
                                                                           dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000
    dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000                                                                  0.308
    dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000                                                                  0.318
    dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000                                                                  0.483
    dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000                                                                  0.000
    dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000                                                                 13.247
                                                                           dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000
    dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000                                                                 13.306
    dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_5_seed_000                                                                 13.284
    dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_4_seed_000                                                                 13.406
    dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_1_seed_000                                                                 13.247
    dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000                                                                  0.000

``` r
library(jsonlite)

# Listing of all PAE JSON files
pae_files <- list.files(path=pth,
                        pattern=".*model.*\\.json",
                        full.names = TRUE)
pae_files
```

    [1] "dimer_23119//dimer_23119_scores_rank_001_alphafold2_multimer_v3_model_2_seed_000.json"
    [2] "dimer_23119//dimer_23119_scores_rank_002_alphafold2_multimer_v3_model_5_seed_000.json"
    [3] "dimer_23119//dimer_23119_scores_rank_003_alphafold2_multimer_v3_model_4_seed_000.json"
    [4] "dimer_23119//dimer_23119_scores_rank_004_alphafold2_multimer_v3_model_1_seed_000.json"
    [5] "dimer_23119//dimer_23119_scores_rank_005_alphafold2_multimer_v3_model_3_seed_000.json"

``` r
pae1 <- read_json(pae_files[1],simplifyVector = TRUE)
pae5 <- read_json(pae_files[5],simplifyVector = TRUE)

attributes(pae1)
```

    $names
    [1] "plddt"   "max_pae" "pae"     "ptm"     "iptm"   

``` r
plot.dmat(pae5$pae, 
          xlab="Residue Position (i)",
          ylab="Residue Position (j)",
          grid.col = "black",
          zlim=c(0,30))
```

![](class11_alphafold_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
plot.dmat(pae1$pae, 
          xlab="Residue Position (i)",
          ylab="Residue Position (j)",
          grid.col = "black",
          zlim=c(0,30))
```

![](class11_alphafold_files/figure-commonmark/unnamed-chunk-11-1.png)

# Score Residue Conservation from alignment file

AlphaFold returns it’s large elignment file used for analysis. Here we
read this file and score conservation per position.

``` r
aln_file <- list.files(path=pth,
                       pattern=".a3m$",
                        full.names = TRUE)
aln_file
```

    [1] "dimer_23119//dimer_23119.a3m"

Read the alignment file

``` r
aln <- read.fasta(aln_file[1], to.upper = TRUE)
```

    [1] " ** Duplicated sequence id's: 101 **"
    [2] " ** Duplicated sequence id's: 101 **"

We can score residue conservation in the alignment with the conserv()
function.

``` r
sim <- conserv(aln)
```

``` r
plotb3(sim[1:99], sse=trim.pdb(pdb, chain="A"),
       ylab="Conservation Score")
```

![](class11_alphafold_files/figure-commonmark/unnamed-chunk-15-1.png)

Find the consensus sequence at a very high cut-off to find invariant
residues

``` r
con <- consensus(aln, cutoff = 0.9)
con$seq
```

      [1] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [19] "-" "-" "-" "-" "-" "-" "D" "T" "G" "A" "-" "-" "-" "-" "-" "-" "-" "-"
     [37] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [55] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [73] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [91] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [109] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [127] "-" "-" "-" "-" "-" "-"
