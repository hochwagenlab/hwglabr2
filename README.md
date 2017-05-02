# hwglabr2
### Hochwagen lab R package 2

Compilation of functions used frequently by Hochwagen lab members. This is a rewrite
of the original packege `hwglabr` 
(GitHub [repo](https://github.com/hochwagenlab/hwglabr)), focused on using 
Bioconductor's `GRanges` object class.

#### Installation

You can install the package directly from the source code on GitHub. For that
you will need Hadley Wickham's `devtools` R package:

``` r
install.packages("devtools")
```

Once you have `devtools` you can install and load `hwglabr2`:

``` r
devtools::install_github("hochwagenlab/hwglabr2")
library(hwglabr2)
```

#### Documentation

Use the package GitHub [repo](https://github.com/hochwagenlab/hwglabr2) and the
[documentation website](http://www.nyu.edu/projects/hochwagen/hwglabr2_docs/).

Function documentation is also accessible within R in the standard way, by
typing one of the following:

``` r
help("function_name")

?function_name
```


