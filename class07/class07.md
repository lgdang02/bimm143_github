# Class 7: Machine Learning 1
Linh Dang (A16897764)

- [Clustering](#clustering)
  - [K-means](#k-means)
  - [Hierarchical Clustering](#hierarchical-clustering)
- [Principal Compoenent Analysis
  (PCA)](#principal-compoenent-analysis-pca)
  - [Data import](#data-import)
  - [PCA to the rescue](#pca-to-the-rescue)

Today we will explore unsupervised machine learning methods starting
with clustering and dimensionality reduction.

## Clustering

To start let’s make up some data to cluster where we know what the
answer should be. The `rnorm()` function will help up here.

``` r
hist( rnorm(10000, mean=3) )
```

![](class07_files/figure-commonmark/unnamed-chunk-1-1.png)

Return 30 numbers centered on -3

``` r
tmp <- c( rnorm(30, mean=-3),
          rnorm(30, mean=3) )

x <- cbind(x=tmp, y=rev(tmp))

x
```

                   x          y
     [1,] -2.2816674  2.5228197
     [2,] -4.6975852  2.2708137
     [3,] -2.4853465  3.5560243
     [4,] -1.5784031  3.2811384
     [5,] -3.1250424  3.6293584
     [6,] -3.7853902  3.8424528
     [7,] -2.7823375  2.6407128
     [8,] -2.3348644  3.4052784
     [9,] -2.8614229  4.3250347
    [10,] -3.2167484  3.5886997
    [11,] -3.4350718  3.7916438
    [12,] -3.0208166  2.5680062
    [13,] -2.2294688  3.4257272
    [14,] -2.6052761  2.5788086
    [15,] -2.3075138  4.2500018
    [16,] -3.2955165  4.8595054
    [17,] -3.5512875  1.8377747
    [18,] -4.2530931  5.8765206
    [19,] -3.7686349  2.4817848
    [20,] -4.3804821  3.4510825
    [21,] -2.1799600  4.4677435
    [22,] -1.9227042  4.9253785
    [23,]  0.1585469  3.8674478
    [24,] -2.9508642  4.6851877
    [25,] -4.2177419  1.7183510
    [26,] -4.2479295  1.8861837
    [27,] -1.7877578  4.4616056
    [28,] -3.2211432  4.6149257
    [29,] -3.1961360  2.3724072
    [30,] -2.4658708  2.2650041
    [31,]  2.2650041 -2.4658708
    [32,]  2.3724072 -3.1961360
    [33,]  4.6149257 -3.2211432
    [34,]  4.4616056 -1.7877578
    [35,]  1.8861837 -4.2479295
    [36,]  1.7183510 -4.2177419
    [37,]  4.6851877 -2.9508642
    [38,]  3.8674478  0.1585469
    [39,]  4.9253785 -1.9227042
    [40,]  4.4677435 -2.1799600
    [41,]  3.4510825 -4.3804821
    [42,]  2.4817848 -3.7686349
    [43,]  5.8765206 -4.2530931
    [44,]  1.8377747 -3.5512875
    [45,]  4.8595054 -3.2955165
    [46,]  4.2500018 -2.3075138
    [47,]  2.5788086 -2.6052761
    [48,]  3.4257272 -2.2294688
    [49,]  2.5680062 -3.0208166
    [50,]  3.7916438 -3.4350718
    [51,]  3.5886997 -3.2167484
    [52,]  4.3250347 -2.8614229
    [53,]  3.4052784 -2.3348644
    [54,]  2.6407128 -2.7823375
    [55,]  3.8424528 -3.7853902
    [56,]  3.6293584 -3.1250424
    [57,]  3.2811384 -1.5784031
    [58,]  3.5560243 -2.4853465
    [59,]  2.2708137 -4.6975852
    [60,]  2.5228197 -2.2816674

Make a plot of `x`

``` r
plot(x)
```

![](class07_files/figure-commonmark/unnamed-chunk-3-1.png)

### K-means

The main function in “base” R for K-means clustering is called `kmeans`:

``` r
km <- kmeans(x, centers =2)
km
```

    K-means clustering with 2 clusters of sizes 30, 30

    Cluster means:
              x         y
    1 -2.934251  3.448247
    2  3.448247 -2.934251

    Clustering vector:
     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

    Within cluster sum of squares by cluster:
    [1] 62.23937 62.23937
     (between_SS / total_SS =  90.8 %)

    Available components:

    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

The `kmeans()` function return a “list” with 9 components. You can see
the named components of any list with the `attributes()` function.

``` r
attributes(km)
```

    $names
    [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    [6] "betweenss"    "size"         "iter"         "ifault"      

    $class
    [1] "kmeans"

> Q. How mnay points are in each cluster?

``` r
km$size
```

    [1] 30 30

> Q. Cluster assignment/membership vector?

``` r
km$cluster
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

> Q. Cluster centers?

``` r
km$centers
```

              x         y
    1 -2.934251  3.448247
    2  3.448247 -2.934251

> Q. Make a plot of our `kmeans()` results showing cluster assignment
> using different colors for each cluster/group of points and cluster
> centers.

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-9-1.png)

> Q. Run `kmeans()` again on `x` and this cluster in to 4
> groups/clusters and plot the same result figure as above.

``` r
km4 <- kmeans(x, centers=4)
plot(x, col=km4$cluster)
points(km4$centers, col="blue", pch=15, cex=2)
```

![](class07_files/figure-commonmark/unnamed-chunk-10-1.png)

> **key-point**: K-means clustering is super popular but can be
> miss-used. One big limitation is that it can impose a clustering
> pattern on your data even if clear natural grouping don’t exist -
> i.e. it does what you tell it to do in terms of `centers`.

### Hierarchical Clustering

The main function in “base” R for Hierarchical Clustering is called
`hclust()`.

You can’t just pass our dataset as is into `hclust()` you must give
“distance matrix” as input. We can get this from the `dist()` function
in R.

``` r
d <- dist(x)
hc <-hclust(d)
hc
```


    Call:
    hclust(d = d)

    Cluster method   : complete 
    Distance         : euclidean 
    Number of objects: 60 

The result of `hclust()` don’t have a useful `print()` method but do
have a special `plot()` method.

``` r
plot(hc)
abline(h=8, col="red")
```

![](class07_files/figure-commonmark/unnamed-chunk-12-1.png)

To get our main cluster assignment (membership vector) we need to “cut”
the tree at the big goal posts…

``` r
grps <- cutree(hc, h=8)
grps
```

     [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
table(grps)
```

    grps
     1  2 
    30 30 

``` r
plot(x, col=grps)
```

![](class07_files/figure-commonmark/unnamed-chunk-15-1.png)

Hierarchical Clustering is distinct in that the dendogram (tree figure)
can reveal the potential grouping in your data (unlike K-means)

## Principal Compoenent Analysis (PCA)

PCA is a common and highly useful dimensionality reduction technique
used in many field - particularly bioinformatics.

Here we will analyze some data from the UK on the food consumption.

### Data import

``` r
url <- "https://tinyurl.com/UK-foods"
x <- read.csv(url)

head(x)
```

                   X England Wales Scotland N.Ireland
    1         Cheese     105   103      103        66
    2  Carcass_meat      245   227      242       267
    3    Other_meat      685   803      750       586
    4           Fish     147   160      122        93
    5 Fats_and_oils      193   235      184       209
    6         Sugars     156   175      147       139

``` r
rownames(x) <- x[,1]
x <- x[,-1]
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
x <- read.csv(url, row.names =1)
head(x)
```

                   England Wales Scotland N.Ireland
    Cheese             105   103      103        66
    Carcass_meat       245   227      242       267
    Other_meat         685   803      750       586
    Fish               147   160      122        93
    Fats_and_oils      193   235      184       209
    Sugars             156   175      147       139

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](class07_files/figure-commonmark/unnamed-chunk-19-1.png)

One conventional plot that can be useful is callled a “paris” plot.

``` r
pairs(x, col=rainbow(nrow(x)), pch=16)
```

![](class07_files/figure-commonmark/unnamed-chunk-20-1.png)

### PCA to the rescue

The main function in base R for PCA `prcomp()`.

``` r
pca <- prcomp( t(x) )
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3       PC4
    Standard deviation     324.1502 212.7478 73.87622 2.921e-14
    Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

The `prcomp()` function returns a list object of our results with five
attributes/components

``` r
attributes(pca)
```

    $names
    [1] "sdev"     "rotation" "center"   "scale"    "x"       

    $class
    [1] "prcomp"

The two main “result” in here are `pca$x` and `pca$rotation`. The first
of these (`pca$x`) contains the scores of the data on the new PC axis -
we use these to make our “PCA plot”.

``` r
pca$x
```

                     PC1         PC2        PC3           PC4
    England   -144.99315   -2.532999 105.768945 -9.152022e-15
    Wales     -240.52915 -224.646925 -56.475555  5.560040e-13
    Scotland   -91.86934  286.081786 -44.415495 -6.638419e-13
    N.Ireland  477.39164  -58.901862  -4.877895  1.329771e-13

``` r
library(ggplot2)
library(ggrepel)

# Make a plot of pca$x with PC1 vs PC2
ggplot(pca$x) +
  aes(PC1, PC2, label=rownames(pca$x)) +
  geom_point() +
  geom_text_repel()
```

![](class07_files/figure-commonmark/unnamed-chunk-24-1.png)

This show the variation of food between the 4 countries. In this case,
N.Ireland is most dissimilar to the 3 other countries by intepreting the
relationship on the PC (N.Ireland has positive score compared to
England, Wales, and Scotland).

The second major result is contained in the `pca$rotation` object or
component. Let’s plot this to see what PCA is picking up…

``` r
ggplot(pca$rotation) +
  aes(PC1, rownames(pca$rotation)) +
  geom_col()
```

![](class07_files/figure-commonmark/unnamed-chunk-25-1.png)

The plot shows the contribution of variance going in both negative
(England, Wale, and Scotland) and positive direction (N.Ireland). We can
see that N.Ireland contribute to the right side of the plot for soft
drinks and fresh potatoes, while the other 3 countries contribute to the
lest side of the plot for fresh fruit and alcoholic drinks, etc.
