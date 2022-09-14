## Comments from initial CRAN submission

* If there are references describing the methods in your package, ...

> Added one reference.


* Please rather use the Authors@R field ...

> Done


* Please add \value to .Rd files ...

> Done


* Please ensure that you do not use more than 2 cores in your examples, vignettes, etc. 

> We use parallel::mclapply(...,mc.cores=cores) and double-checked that cores = 2.


* \dontrun{} in examples:

> The example for ReadGRAND takes > 5 sec to execute and is therefore wrapped into \dontrun; it only shows the usage, but does not really produce informative outputs
> The example for PairwiseDESeq2 also takes > 5 sec to execute and is therefore wrapped into \dontrun


* You are using installed.packages() in your code.

> Removed.




## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:
* checking CRAN incoming feasibility ... NOTE                                                                                                                                 
Maintainer: ‘Florian Erhard <Florian.Erhard@uni-wuerzburg.de>’                                                                                                                
                                                                                                                                                                              
New submission

Found the following (possibly) invalid URLs:                                                                                                                                  
  URL: GO:BP                                                                                                                                                                  
    From: inst/doc/differential-expression.html                                                                                                                               
    Message: Invalid URI scheme                                                                                                                                               
  URL: GO:CC                                                                                                                                                                  
    From: inst/doc/differential-expression.html                                                                                                                               
    Message: Invalid URI scheme                                                                                                                                               
  URL: GO:MF                                                                                                                                                                  
    From: inst/doc/differential-expression.html                                                                                                                               
    Message: Invalid URI scheme                                                                                                                                               

> these are not URLs but database identifiers
    

* checking R code for possible problems ... NOTE                                                                                                                            
<many "no visible binding for global variable">

> This package implements many functions that use ggplot to generate plots. Symbols
> used in the call to aes() are recognized as global variable by the check, but indeed
> are evaluated in an environment having the data frame used in ggplot.


