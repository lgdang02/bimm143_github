# Class 6: R functions
Linh Dang (PID: A16897764)

- [1. Function basics](#1-function-basics)
- [2. Generate DNA sequence](#2-generate-dna-sequence)
- [3. Generate Protein function](#3-generate-protein-function)

## 1. Function basics

Let’s start writing our first silly function to add some numbers:

Every R function has 3 things:

- name (we get to pick this)
- input arguments (there can loads of these separated by a comma)
- the body (the R code that does the work)

``` r
add <- function(x, y=10, z=0){
  x + y + z
}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

``` r
add( x=c(1,2,3,4), y=100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=10`) in the function defination.

``` r
add(1, 100, 10)
```

    [1] 111

> Q. Write a function to returna DNA sequence of a user specified
> length? Call it `generate_dna()`

The `sample()` function can help here

``` r
#generate_dna <- function(size=5) { }

students <- c("jeff", "jeremy", "peter")

sample(students, size = 5, replace = TRUE)
```

    [1] "peter"  "jeremy" "jeff"   "jeff"   "jeremy"

## 2. Generate DNA sequence

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")
sample(bases, size = 10, replace = TRUE)
```

     [1] "C" "A" "A" "G" "A" "A" "C" "T" "C" "T"

Now I have a working ‘snippet’ of code I can use this as the body of my
first function version here:

``` r
generate_dna <- function(size=5) {
  bases <- c("A", "C", "G", "T")
  sample(bases, size = size, replace = TRUE)
  }
```

``` r
generate_dna(100)
```

      [1] "G" "T" "C" "C" "A" "T" "A" "G" "A" "G" "T" "C" "C" "G" "A" "C" "C" "A"
     [19] "C" "G" "A" "A" "T" "T" "A" "A" "T" "A" "C" "G" "T" "C" "C" "G" "G" "T"
     [37] "G" "A" "A" "A" "C" "A" "T" "C" "G" "C" "C" "G" "T" "T" "T" "A" "A" "A"
     [55] "A" "G" "T" "A" "C" "A" "A" "G" "G" "A" "A" "A" "G" "A" "G" "G" "A" "G"
     [73] "A" "C" "G" "G" "G" "G" "G" "G" "C" "A" "C" "A" "C" "A" "T" "T" "A" "C"
     [91] "G" "T" "A" "T" "C" "A" "T" "G" "T" "T"

``` r
generate_dna()
```

    [1] "C" "T" "G" "T" "A"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(size = 5, together = TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size = size, replace = TRUE)
  
  if(together) {
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

``` r
generate_dna(together=F)
```

    [1] "A" "T" "A" "C" "C"

## 3. Generate Protein function

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length?

We can get the set of 20 natural amino-acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

``` r
generate_protein <- function(size=6, together=TRUE ) {
  
  ## Get the 20 amino-acids as a vector
  aa <- bio3d::aa.table$aa1[1:20]
  sequence <- sample(aa, size, replace=TRUE)
  
  ## Optionally return a single element string
  if(together){
    sequence <- paste(sequence, collapse="")
  }
  return(sequence)
}
```

``` r
generate_protein(together=F)
```

    [1] "A" "Q" "M" "Y" "N" "R"

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
generate_protein(7)
```

    [1] "GRNLYFA"

``` r
generate_protein(8)
```

    [1] "EGMKVGHF"

``` r
generate_protein(9)
```

    [1] "ASWMIPLSG"

``` r
#generate_protein(size=6:12)
```

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. for a loop) or by
using the R **apply** family of utility functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "KQQCDA"       "CTCTGYC"      "HCGQKTQF"     "VGLKTGIMT"    "IFDLWYRPNQ"  
    [6] "AIYVNKVTYTA"  "MKSGSLVVAAPK"

It would be cool and useful if I could get FASTA format output

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "KNQFVY"       "PDPPARY"      "CFTKFQPH"     "GACPIPDIQ"    "DFQNCRIAKH"  
    [6] "PYYSEHTAVED"  "QNVTWDKALVKH"

``` r
cat(ans, sep="\n")
```

    KNQFVY
    PDPPARY
    CFTKFQPH
    GACPIPDIQ
    DFQNCRIAKH
    PYYSEHTAVED
    QNVTWDKALVKH

I want this to look like FASTA format with an ID line e.g.

    >ID.6
    MNHTPI
    >ID.7
    HDLPFIW
    >ID.8
    LDAYVHKD

The functions `paste()` and `cat()` can help ups here…

``` r
cat( paste(">ID.", 6:12, "\n", ans, sep=""), sep="\n" )
```

    >ID.6
    KNQFVY
    >ID.7
    PDPPARY
    >ID.8
    CFTKFQPH
    >ID.9
    GACPIPDIQ
    >ID.10
    DFQNCRIAKH
    >ID.11
    PYYSEHTAVED
    >ID.12
    QNVTWDKALVKH

``` r
id.line <- paste(">ID.", 6:12, sep="")
id.line
```

    [1] ">ID.6"  ">ID.7"  ">ID.8"  ">ID.9"  ">ID.10" ">ID.11" ">ID.12"

``` r
id.line <- paste(">ID.", 6:12, sep="")
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n", file="myseq.fa")
```

> Q. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

I BLASTp searched my FASTA format sequences against NR and found that
length 6, 7 are not unique and can be found in the databases with 100%
coverage and 100% identity.

Random sequences of length 8 and above are unique and can’t be found in
the databases.
