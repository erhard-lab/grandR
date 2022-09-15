## Relevant comments from earlier CRAN submission

* Please ensure that you do not use more than 2 cores in your examples, vignettes, etc. 

> We use parallel::mclapply(...,mc.cores=cores) and double-checked that cores = 2.


* \dontrun{} in examples:

> The example for ReadGRAND takes > 5 sec to execute and is therefore wrapped into \dontrun
> The example for PairwiseDESeq2 also takes > 5 sec to execute and is therefore wrapped into \dontrun




## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking CRAN incoming feasibility ... NOTE                                                                                                                                 
Maintainer: ‘Florian Erhard <Florian.Erhard@uni-wuerzburg.de>’                                                                                                                
                                                                                                                                                                              
New submission


