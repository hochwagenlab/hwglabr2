[![Build Status](https://travis-ci.org/hochwagenlab/hwglabr2.svg?branch=master)](https://travis-ci.org/hochwagenlab/hwglabr2)

# hwglabr2
### Hochwagen lab R package 2

Compilation of functions used frequently by Hochwagen lab members. This is a
rewrite of the original package `hwglabr` 
(GitHub [repo](https://github.com/hochwagenlab/hwglabr)), designed to leverage 
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

See the
[Analysis recipes](https://github.com/hochwagenlab/hwglabr2/wiki/Analysis-recipes)
GitHub wiki page for example usage of `hwglabr2` to perform some common analysis.

For function documentation, use the package GitHub [repo](https://github.com/hochwagenlab/hwglabr2)
and the [documentation website](http://www.nyu.edu/projects/hochwagen/hwglabr2_docs/)
built with `pkgdown`.

Function documentation is also accessible within R in the standard way, by typing
one of the following:

``` r
help("function_name")

?function_name
```

## License

This project is licensed under the terms of the MIT license.
See [LICENSE](LICENSE) file for details.
