# Comments for submission of grandR 0.2.6
## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

## CRAN Package Check Results
There were no ERRORs or WARNINGs.

There was 1 additional NOTE:
  Package suggested but not available for checking: ‘monocle’
  
This is a Bioconductor package.



# Comments for submission of grandR 0.2.5
## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:
* checking CRAN incoming feasibility ... NOTE                                                                             
Maintainer: ‘Florian Erhard <Florian.Erhard@informatik.uni-regensburg.de>’                                                
                                                                                                                          
New maintainer:                                                                                                           
  Florian Erhard <Florian.Erhard@informatik.uni-regensburg.de>                                                            
Old maintainer(s):                                                                                                        
  Florian Erhard <Florian.Erhard@uni-wuerzburg.de>                                                                        
  
> I double-checked my new email address


## CRAN Package Check Results
There were no ERRORs or WARNINGs.

There was 1 additional NOTE:
  Package suggested but not available for checking: ‘monocle’
  
This is a Bioconductor package.


# Comments for submission of grandR 0.2.2
## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 


# Comments for submission of grandR 0.2.1
## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 




# Relevant comments from earlier CRAN submission

* Please wrap examples that are generally executable but need > 5 sec in \donttest{}. \dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.

> Changed to \donttest in both cases


* Please ensure that you do not use more than 2 cores in your examples, vignettes, etc. 

> We use parallel::mclapply(...,mc.cores=cores) and double-checked that cores = 2.


* \dontrun{} in examples:

> The example for ReadGRAND takes > 5 sec to execute and is therefore wrapped into \dontrun > now changed to \donttest
> The example for PairwiseDESeq2 also takes > 5 sec to execute and is therefore wrapped into \dontrun > now changed to \donttest




## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:
* checking CRAN incoming feasibility ... NOTE                                                                                                                                 
Maintainer: ‘Florian Erhard <Florian.Erhard@uni-wuerzburg.de>’                                                                                                                
                                                                                                                                                            
New submission

* checking examples (25.4s)
   Examples with CPU (user + system) or elapsed time > 5s
                    user system elapsed
   PairwiseDESeq2 7.242  0.212  7.464

> This is wrapped in \donttest
