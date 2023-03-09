

#' Density estimation in 2d
#'
#' Estimate point densities on a regular grid for.
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param facet factor: estimate for each unique factor; can be NULL
#' @param n size of the grid
#' @param margin one of 'n','x' or 'y'; should the density be computed along both axes ('n'), or along 'x' or 'y' axis only
#'
#' @return a density value for each point
#' @export
#'
#' @concept helper
density2d=function(x, y, facet=NULL, n=100, margin='n') {
  bandwidth.nrd.ex=function (x)
  {
    r <- range(x)
    h <- (r[2L] - r[1L])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
  }

  if (is.null(facet)) {
        use=is.finite(x+y)
        if (min(x[use])==max(x[use])) {
          x[use]=0:1
        }
        if (min(y[use])==max(y[use])) {
          y[use]=0:1
        }
        bw.x=MASS::bandwidth.nrd(x[use])
        if (bw.x==0) bw.x=bandwidth.nrd.ex(x[use])
        bw.y=MASS::bandwidth.nrd(y[use])
        if (bw.y==0) bw.y=bandwidth.nrd.ex(y[use])
        d=MASS::kde2d(x[use], y[use], h=c(bw.x,bw.y), n=n)
        if (margin=='x') d$z=d$z/apply(d$z,1,max)
        else if (margin=='y') d$z=t(t(d$z)/apply(d$z,2,max))
        else d$z=d$z/max(d$z)
        r=rep(NA,length(x))
        r[use]=d$z[cbind(findInterval(x[use], d$x),findInterval(y[use], d$y))]
	r=r/max(r,na.rm=T)
        return(r)
    }
    re<-rep(NA,length(x))
    for (f in unique(facet)) re[f==facet]=density2d(x[f==facet],y[f==facet],n=n)
    re
}


#' Make a PCA plot
#' @param data the grandR object that contains the data to plot
#' @param mode.slot the mode and slot of data to plot; slot in the grandr object (eg "count")
#' @param ntop how many genes to use
#' @param aest parameter to set the visual attributes
#' @param x number of principal component to show on the x axis (numeric)
#' @param y number of principal component to show on the y axis (numeric)
#' @param columns which columns (i.e. samples or cells) to perform PCA on (see details)
#' @param do.vst perform a variance stabilizing transformation for count data?
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a PCA plot
#' @export
#'
#' @concept globalplot
PlotPCA=function(data, mode.slot=DefaultSlot(data), ntop=500,aest=NULL,x=1,y=2,columns=NULL,do.vst=TRUE) {

  aest=setup.default.aes(data,aest)

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  mat=as.matrix(GetTable(data,type=mode.slot,columns=columns,ntr.na = FALSE))

	cd=do.call("rbind",lapply(1:length(mode.slot), function(i) cbind(data$coldata,data.frame(type=mode.slot[i]))))


	mat=mat[,columns]
	cd=cd[columns,]

	rm.na=!apply(is.na(mat),2,sum)==nrow(mat)
	mat=mat[,rm.na]
	cd=cd[rm.na,]
	if (do.vst) {
	  checkPackages("DESeq2")
	  mat <- DESeq2::vst(cnt(mat))
	}

	checkPackages("matrixStats")
	rv <- matrixStats::rowVars(mat)

  	select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  	pca <- prcomp(t(mat[select,]))
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
	d <- as.data.frame(pca$x)
	names(d)=paste0("PC",1:dim(d)[2])
	d=cbind(d, cd)
	ggplot(d,utils::modifyList(aes(!!sym(paste0("PC",x)),!!sym(paste0("PC",y))),aest))+cowplot::theme_cowplot()+
	  geom_point(size=3)+xlab(paste0("PC",x,": ",round(percentVar[x] * 100),"% variance"))+ylab(paste0("PC",y,": ",round(percentVar[y] * 100),"% variance"))+coord_fixed()
}

Transform=function(name,label=NULL){
  FUN=get(paste0("Transform.",name))()
  function(mat) {
    re=FUN(mat)
    if (!is.null(label)) attr(re,"label")=label
    re
  }
}

#' Transformations for PlotHeatmap
#'
#' Functions to perform transformations on the matrix used for \link{PlotHeatmap}.
#'
#' @param label label that is used for the heatmap legend
#' @param center perform centering when computing Z scores (see \link{scale})
#' @param scale perform scaling when computing Z scores (see \link{scale})
#' @param LFC.fun function to compute log fold changes (default: \link[lfc]{PsiLFC}, other viable option: \link[lfc]{NormLFC})
#' @param columns which columns (i.e. samples or cells) to use as reference when computing log fold changes (see details)
#' @param ... further parameters passed down to LFC.fun
#'
#' @details These functions should be used as transform parameter to \link{PlotHeatmap}. Available data transformations are
#' \itemize{
#'  \item{transform=Transform.Z(): compute z scores for each row; you can omit the usual centering or scaling by setting the respective parameters to false; see \link{scale}}
#'  \item{transform=Transform.VST(): do a variance stabilizing transformation using \link[DESeq2]{vst}}
#'  \item{transform=Transform.logFC(): compute log2 fold changes to one or several reference columns; see below how to define them; fold changes are computed using the lfc package)}
#'  \item{transform=Transform.no(): do not transform}
#' }
#'
#' @details The label to be used in the heatmap legend can be changed by specifying the label parameter.
#'
#' @details For Transform.logFC, columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#'
#' @return A function that transforms a matrix.
#'
#' @export
#'
#' @concept globalplot
Transform.no=function(label=" ") function(m) {re=m; attr(re,"label")=label; re}
#' @rdname Transform.no
#' @export
Transform.Z=function(label="z score",center=TRUE,scale=TRUE) function(m) {re=t(scale(t(m),center=center,scale=scale)); attr(re,"label")=label; re}
#' @rdname Transform.no
#' @export
Transform.VST=function(label="VST") function(m) {
  if (nrow(m)<1000) re = DESeq2::varianceStabilizingTransformation(cnt(m)) else re=DESeq2::vst(cnt(m))
  attr(re,"label")=label
  re
  }
#' @rdname Transform.no
#' @export
Transform.logFC=function(label="log2 FC",LFC.fun=NULL,columns=NULL,...) {
  function(m) {
    if (is.null(columns)) columns=1:ncol(m)
    ref=rowMeans(m[,columns,drop=F])
    if (is.null(LFC.fun)) LFC.fun=lfc::PsiLFC
    re=apply(m,2,function(v) LFC.fun(v,ref,normalizeFun=function(vv) vv))
    attr(re,"label")=label
    re
  }
}


make.continuous.colors=function(values,colors=NULL,breaks=NULL) {
  if (quantile(values,0.25,na.rm=TRUE)<0) {
    quant=c(50,95) #c(seq(0,1,length.out=nq+1)[c(-1,-nq-1)]*100,95)
    if (length(breaks)==1) {
      quant=c(seq(0,1,length.out=breaks+1)[c(-1,-breaks-1)]*100,95)
      breaks=NULL
    }
    if (is.null(breaks)) {
      upper=quantile(values[values>0],quant/100,na.rm=TRUE)
      lower=quantile(-values[values<0],quant/100,na.rm=TRUE)
      breaks=c(-rev(pmax(upper,lower)),0,pmax(upper,lower))
    }
    if (is.null(colors)) colors="RdBu"
  } else {
    quant=c(5,25,50,75,95)
    if (length(breaks)==1) {
        quant=c(5,seq(0,1,length.out=breaks)[c(-1,-breaks)]*100,95)
        breaks=NULL
    }
    if (is.null(breaks)) {
      breaks=quantile(values,quant/100,na.rm=TRUE)
    }
    if (is.null(colors)) colors="YlOrRd"
  }

  if (length(colors)==1) checkPackages(c("RColorBrewer","viridisLite"))

  if (length(colors)==1 && colors %in% rownames(RColorBrewer::brewer.pal.info)) {
    colors = RColorBrewer::brewer.pal(length(breaks),colors)
  } else if (length(colors)==1) {
    colors = viridisLite::viridis(length(breaks),option = colors)
  }
  list(breaks = breaks,colors=colors)
}


#' Create heatmaps from grandR objects
#'
#' Convenience method to compare among more two variables (slot data or analyses results).
#'
#' @param data the grandR object that contains the data to plot
#' @param type Either a mode.slot (see details) or a regex to be matched against analysis names. Can also be a vector
#' @param columns a vector of columns (either condition/cell names if the type is a mode.slot, or names in the output table from an analysis; use \link{Columns}(data,<analysis>) to learn which columns are available); all condition/cell names if NULL
#' @param genes the genes to be included in the plot (default: all genes)
#' @param summarize Should replicates by summarized? Can only be specified if columns is NULL; either a summarization matrix (\link{GetSummarizeMatrix}) or TRUE (in which case \link{GetSummarizeMatrix}(data) is called)
#' @param transform apply a transformation to the selected data; can be a function, or a character (see details)
#' @param cluster.genes should genes be clustered?
#' @param cluster.columns should samples (or cells) be clustered?
#' @param label.genes should genes be labeled?
#' @param xlab The names to show at the x axis (only works if type is a single slot)
#' @param breaks vector of color breaks; can be NULL (see details)
#' @param colors an RColorBrewer palette name; can be NULL (see details)
#' @param title the title for the plot; can be NULL
#' @param return.matrix if TRUE, return a list containing the data matrix and the heatmap instead of the heatmap alone
#' @param na.to convert NA values in the matrix to this value immediately before computing the heatmap
#' @param ... additional parameters forwarded to \link[ComplexHeatmap]{Heatmap}
#'
#' @details This is just a convenience function which
#' \enumerate{
#' \item{Calls \link{GetTable} with the parameter \code{type,columns,summarize,genes}}
#' \item{Transforms the returned table using the \code{transform} parameter}
#' \item{Determines reasonable colors using \code{breaks} and \code{colors}}
#' \item{and then calls ComplexHeatmap::Heatmap}
#' }
#'
#' @details \code{type} and \code{columns} can refer to values from data slots values from analyses (and can be mixed).
#' If there are types from both data and analyses, columns must be NULL.
#' Otherwise columns must either be condition/cell names (if type refers to one or several data slots), or regular expressions
#' to match against the names in the analysis tables.
#'
#' @details Columns definitions for data slots can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total)
#' followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @details The transform parameter either is a function that transforms a matrix (which can conveniently be done using the Transform.XXX functions described next), or
#' a character (which must be the XXX to find such a function). Available data transformations are
#' \itemize{
#'  \item{transform=Transform.Z() or transform="Z": compute z scores for each row (see \link{Transform.Z})}
#'  \item{transform=Transform.VST() or transform="VST": do a variance stabilizing transformation (see \link{Transform.VST})}
#'  \item{transform=Transform.logFC() or transform="logFC": compute log2 fold changes to one or several reference columns; which must be defined via parameters (see \link{Transform.logFC})}
#'  \item{transform=Transform.no() or transform="no": do not transform (see \link{Transform.no})}
#' }
#'
#' @details Reasonable coloring is chosen depending on the value distribution in the matrix. If the values are zero centered (e.g. z scores or most often log fold changes), then
#' by default the 50% and 95% quantiles of all positive and all negative values are determined. Let q95 be the 95% quantile with the larger absolute value, and q50 likewise the 50%
#' quantile with the larger value. The breaks are -q90,q50,0,q50,q90, and, by default, the red to blue "RdBu" palette from RColorBrewer is taken. If the values are not zero centered,
#' the 5%,25%,50%,75%, and 95% quantiles are used as breaks and the yellow-orange-red ("YlOrRd") palette is taken. Breaks can also be specified (as absolute values).
#'
#' @details xlab can be given as a character vector or an expression that evaluates into a character vector.
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables.
#'
#' @return a ComplexHeatmap object
#'
#' @seealso \link{GetTable},\link[ComplexHeatmap]{Heatmap}
#'
#' @export
#'
#' @concept globalplot
PlotHeatmap=function(data,
                     type=DefaultSlot(data),
                     columns=NULL,
                     genes=NULL,
                     summarize=NULL,
                     transform="Z",
                     cluster.genes=TRUE,
                     cluster.columns=FALSE,
                     label.genes=NULL,
                     xlab=NULL,
                     breaks=NULL,
                     colors=NULL,
                     title=NULL,return.matrix=FALSE,
                     na.to=NA,...) {

  checkPackages(c("ComplexHeatmap","circlize"))


  mode.slot=check.mode.slot(data,type)
  if (any(mode.slot)) {
    columns=substitute(columns)
    columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
    columns=Columns(data,columns)

    if (length(type)==1) {
      xlab=substitute(xlab)
      xlab=if (is.null(xlab)) NULL else rep(eval(xlab,Coldata(data),parent.frame()),length.out=length(columns))
    }
  }

  mat=as.matrix(GetTable(data,type=type,genes = genes,columns=columns,summarize = summarize,ntr.na = FALSE))
  if (is.character(transform)) transform=Transform(transform)
  mat=transform(mat)


  name=attr(mat,"label")
  col=make.continuous.colors(mat,colors = colors,breaks=breaks)
  col=circlize::colorRamp2(breaks=col$breaks,colors=col$colors)

  if (length(xlab)==ncol(mat)) colnames(mat)=xlab
  if (is.null(label.genes)) label.genes=nrow(mat)<=50

  mat[is.na(mat)]=na.to
  hm=ComplexHeatmap::Heatmap(mat,name=name,
                          cluster_rows = cluster.genes,
                          cluster_columns = cluster.columns,
                          show_row_names = label.genes,
                          col = col,
                          column_title=if (is.null(title)) "" else title,
                          row_title=sprintf("n=%d",nrow(mat)),
                          ...)

  if (return.matrix) return(list(Matrix=mat,Heatmap=hm))
  hm
}

#PlotTestOverlap=function(data,names=NULL,alpha=0.05,type=c("venn","euler")) {
#  # R CMD check guard for non-standard evaluation
#  name <- NULL#

#  mat=GetAnalysisTable(data,gene.info=FALSE,genes=names,columns='^Q$')
#	df=setNames(as.data.frame(mat<alpha & !is.na(mat)),gsub(".Q$","",names(mat)))
#	pl=switch(type[1],euler=eulerr::euler(df),venn=eulerr::venn(df))
#	plot(pl,main=name)
#}

#' Formatting function for correlations
#'
#' Returns a function that takes x and y and returns a formatted output to describe the correlation of x and y
#'
#' @param method how to compute correlation coefficients (can be pearson, spearman or kendall)
#' @param n.format format string for the number of data points (see \link{sprintf}); can be NULL (don't output the number of data points)
#' @param coeff.format format string for the correlation coefficient (see \link{sprintf}); can be NULL (don't output the correlation coefficient)
#' @param p.format  format string for the P value (see \link{sprintf}); can be NULL (don't output the P value)
#'
#' @details Use this for the \code{correlation} parameter of \link{PlotScatter}
#'
#' @return a function
#' @export
#'
#' @examples
#'
#' set.seed(42)
#' data <- data.frame(u=runif(500))  # generate some correlated data
#' data$x <- rnorm(500,mean=data$u)
#' data$y <- rnorm(500,mean=data$u)
#'
#' fun <- FormatCorrelation()
#' fun(data$x,data$y)
#'
#' fun <- FormatCorrelation(method="spearman",p.format="%.4g")
#' fun(data$x,data$y)
#'
#' @concept globalplot
FormatCorrelation=function(method="pearson",n.format=NULL,coeff.format="%.2f",p.format="%.2g") {
  function(x,y) {
    if (length(x)!=length(y)) stop("Cannot compute correlation, unequal lengths!")
    use=is.finite(x)&is.finite(y)
    if (sum(use)<length(x)) warning(sprintf("Removed %d/%d non finite values while computing correlation!",sum(!use),length(x)))
    x=x[use]
    y=y[use]
    cc=cor.test(x,y,method=method)
    p.name=switch(method,pearson="R",spearman="\U03C1",kendall="\U03C4")
    formatted.n=if (!is.null(n.format)) sprintf(sprintf("n=%s",n.format),length(x))
    formatted.p=if (!is.null(p.format)) sprintf(sprintf("p%s%s",if (cc$p.value<2.2e-16) "<" else "=",p.format),if (cc$p.value<2.2e-16) 2.2e-16 else cc$p.value)
    formatted.coeff=if (!is.null(coeff.format)) sprintf(sprintf("%s=%s",p.name,coeff.format),cc$estimate)
    paste(c(formatted.n,formatted.coeff,formatted.p),collapse="\n")
  }
}


#' Make a scatter plot
#'
#' Convenience method to compare two variables (slot data or analyses results).
#'
#' @param data the grandR object (can also be a plain data frame)
#' @param x an expression to compute the x value or a character corresponding to a sample (or cell) name or a fully qualified analysis result name (see details)
#' @param y an expression to compute the y value or a character corresponding to a sample (or cell) name or a fully qualified analysis result name (see details)
#' @param analysis the name of an analysis table (can be NULL; see details)
#' @param xcol a character corresponding to a sample (or cell) name or a fully qualified analysis result name (see details)
#' @param ycol a character corresponding to a sample (or cell) name or a fully qualified analysis result name (see details)
#' @param xlab the label for x (can be NULL, then the x parameter is used)
#' @param ylab the label for y (can be NULL, then the y parameter is used)
#' @param log if TRUE, use log scales for x and y axis
#' @param log.x  if TRUE, use log scale for the x axis
#' @param log.y  if TRUE, use log scale for the y axis
#' @param remove.outlier configure how outliers are selected (is the coef parameter to \link[grDevices]{boxplot.stats}); can be FALSE, in which case no points are considered outliers (see details)
#' @param lim define the both x and y axis limits (vector of length 2 defining the lower and upper bound, respectively)
#' @param xlim define the x axis limits (vector of length 2 defining the lower and upper bound, respectively)
#' @param ylim define the y axis limits (vector of length 2 defining the lower and upper bound, respectively)
#' @param size the point size to use
#' @param cross add horizontal and vertical lines through the origin?
#' @param diag if TRUE, add main diagonal; if numeric vector, add these diagonals
#' @param genes restrict to these genes; can be either numeric indices, gene names, gene symbols or a logical vector
#' @param highlight highlight these genes; can be either numeric indices, gene names, gene symbols or a logical vector (see details)
#' @param label label these genes; can be either numeric indices, gene names, gene symbols or a logical vector (see details)
#' @param label.repel force to repel labels from points and each other (increase if labels overlap)
#' @param facet an expression (evaluated in the same environment as x and y); for each unique value a panel (facet) is created; can be NULL
#' @param color either NULL (use point density colors), or a name of the \link{GeneInfo} table (use scale_color_xxx to define colors), or a color for all points
#' @param colorpalette either NULL (use default colors), or a palette name from color brewer or viridis
#' @param density.margin for density colors, one of 'n','x' or 'y'; should the density be computed along both axes ('n'), or along 'x' or 'y' axis only
#' @param density.n how many bins to use for density calculation (see \link[MASS]{kde2d})
#' @param rasterize use ggrastr to rasterize points? (can be NULL, see details)
#' @param correlation a function to format correlation statistics to be annotated (see details)
#' @param correlation.x x coordinate to put the correlation annotation in the plot (see details)
#' @param correlation.y y coordinate to put the correlation annotation in the plot (see details)
#' @param correlation.hjust x adjustment to put the correlation annotation in the plot (see details)
#' @param correlation.vjust y adjustment to put the correlation annotation in the plot (see details)
#' @param layers.below list of ggplot geoms to add before adding the layer containing the points
#'
#' @details Both the x and y parameter are either expressions or names. Names are either sample (or cell, in case of single cell experiments) names or
#' fully qualified analysis results (analysis name followed by a dot and the analysis result table column). If the analysis parameter is given, the analysis
#' name must be omitted from x and y. These names can be used within expressions using non-standard evaluation.
#' Defining by names only works with character literals like "kinetics.Synthesis", but if you give an expression (e.g. a variable name that contains a character),
#' the situation is more complicated, since PlotScatter will try to evaluate this for defining the values, not the name of the column. If the expression evaluates
#' into a single character string that is equal to a name (see above!), PlotScatter knows what to do. For more complicated situations that cannot be resolved by this,
#' you can use the xcol and ycol parameters instead of the x and y parameters!
#'
#' @details By default the limits of x and y axis are chosen after removing outliers (using the same algorithm used for \link{boxplot}). Thus, larger numbers filter
#' less stringently. remove.outlier can also be set to FALSE (no outlier filtering). If xlim or ylim are set, this overrides outlier filtering. Points outside of the limits
#' (i.e. outliers or points outside of xlim or ylim) are set to infinity (such that they are shown at the border of the plot in gray)
#'
#' @details By default, all genes are shown. This can be restricted using the \code{genes} parameter (see \link{ToIndex}). It is also possible to highlight a subset of the genes
#' using \code{highlight}. This parameter either describes a subset of the genes (either numeric indices, gene names, gene symbols or a logical vector), in which case these genes
#' are plotted in red and with larger points size, or it can be a list of such vectors. The names of this list must be valid colors. Genes can also be labeled (make sure that this
#' is really only a small subset of the genes).
#'
#' @details When rendering to vector based devices (such as svg or pds), a genome-wide scatterplot often is painfully big (and rendering therefore slow). The \code{rasterize}
#' parameter can be used to automatically rasterize the points only (via the ggrastr package). If this parameter is NULL, ggrastr is used if more than 1000 points are plotted!
#'
#' @details Often scatter plots show that x and y coordinates are correlated. Correlations can be annotated using the \link{FormatCorrelation} function. Most often you will use
#' \code{PlotScatter(data,x,y,correlation=FormatCorrelation())}. To use a different correlation measure, other formats for correlation coefficient and P values or omit one of these
#' statistics, parametrize \code{FormatCorrelation}. Use correlation.x and correlation.y to place the annotation in the plot, and correlation.hjust/correlation.vjust to align the
#' annotation at the given x,y coordinates. Infinite values for correlation.x/correlation.y will put the annotation at the border of the plot.
#'
#' @return a ggplot object with the data frame used as the df attribute
#' @export
#'
#' @concept globalplot
PlotScatter=function(data,
                     x=NULL, y=NULL, analysis=NULL,xcol=NULL,ycol=NULL, xlab=NULL, ylab=NULL,
                     log=FALSE, log.x=log, log.y=log,
                     remove.outlier=1.5, lim=NULL,xlim=lim, ylim=lim,
                     size=0.3,
                     cross=NULL,diag=NULL,
                     genes=NULL,highlight=NULL, label=NULL, label.repel=1,
                     facet=NULL,
                     color=NULL, colorpalette=NULL,
                     density.margin = 'n', density.n = 100,
                     rasterize=NULL,
                     correlation=NULL,correlation.x=-Inf,correlation.y=Inf,correlation.hjust=0.5,correlation.vjust=0.5,
                     layers.below=NULL) {
  if (is.grandR(data)) {
    if (!is.null(analysis) || IsSparse(data)) {
      df=cbind(GetAnalysisTable(data,analyses = analysis,regex=FALSE,prefix.by.analysis = FALSE,gene.info = FALSE),GeneInfo(data))
    }else{
      df=cbind(GetAnalysisTable(data,gene.info = FALSE),GetTable(data,type=DefaultSlot(data)),GeneInfo(data))
    }
    if (!is.null(genes)) df=df[ToIndex(data,genes),]
  } else {
    df=as.data.frame(data)
  }

  facet=substitute(facet)
  if (!is.null(facet)) {
    df$facet=eval(facet,df,parent.frame())
  }

  if (!is.data.frame(df)) stop("df must be a data frame (or at least coercable into a data frame)")
  adaptInf=function(df,rx,ry) {
    # workaround to also "brush" infinity points at the border of the plane
    if (log.x) rx=log10(rx)
    rx=c(rx[1]-0.04*(rx[2]-rx[1]),rx[2]+0.04*(rx[2]-rx[1]))
    if (log.x) rx=10^rx
    df$A[is.infinite(df$A) & df$A>0]=rx[2]
    df$A[is.infinite(df$A) & df$A<0]=rx[1]

    if (log.y) ry=log10(ry)
    ry=c(ry[1]-0.04*(ry[2]-ry[1]),ry[2]+0.04*(ry[2]-ry[1]))
    if (log.y) ry=10^ry
    df$B[is.infinite(df$B) & df$B>0]=ry[2]
    df$B[is.infinite(df$B) & df$B<0]=ry[1]

    df
  }


  rn=rownames(df)

  x=substitute(x)
  if (is.null(x) && is.null(xcol))  {
    A=df[[1]]
    x=1
  } else if (is.null(x)) {
    A=df[[xcol]]
    x=xcol
  } else if (is.null(xcol)) {
    A=if (is.character(x) || is.numeric(x)) df[[x]] else eval(x,df,parent.frame())
    if (length(A)==1 && is.character(A)) {
      if (A %in% names(df)) A=df[[A]] else stop("Cannot make heads or tails out of the x parameter. Try to use xcol to access a specific column!")
    }
    if (length(A)==1) stop("Cannot make heads or tails out of the x parameter. Try to use xcol to access a specific column!")
  } else {
    stop("You must not specify both x and xcol!")
  }

  y=substitute(y)
  if (is.null(y) && is.null(ycol))  {
    B=df[[2]]
    y=2
  } else if (is.null(y)) {
    B=df[[ycol]]
    y=ycol
  } else if (is.null(ycol)) {
    B=if (is.character(y) || is.numeric(y)) df[[y]] else eval(y,df,parent.frame())
    if (length(B)==1 && is.character(B)) {
      if (B %in% names(df)) B=df[[B]] else stop("Cannot make heads or tails out of the y parameter. Try to use ycol to access a specific column!")
    }
    if (length(B)!=nrow(df)) stop("Cannot make heads or tails out of the y parameter. Try to use xcol to access a specific column!")
  } else {
    stop("You must not specify both y and ycol!")
  }

  if (is.null(xlab)) { if (is.character(x)) xlab=x else if (is.numeric(x)) xlab=names(df)[x] else xlab=deparse(x)}
  if (is.null(ylab)) { if (is.character(y)) ylab=y else if (is.numeric(y)) ylab=names(df)[y] else ylab=deparse(y)}
  if (is.na(xlab) || xlab=="") xlab=NULL
  if (is.na(ylab) || ylab=="") ylab=NULL

  df$A=A
  df$B=B
  rownames(df)=rn


  df$A.trans=if(log.x) log10(df$A) else df$A
  df$B.trans=if(log.y) log10(df$B) else df$B
  df$A.trans.unclipped=df$A.trans
  df$B.trans.unclipped=df$B.trans

  set.coord=remove.outlier!=FALSE || !is.null(xlim) || !is.null(ylim)
  if (set.coord) {
    if (is.null(xlim) && remove.outlier) {
      xlim=grDevices::boxplot.stats(df$A.trans[is.finite(df$A.trans)],coef=remove.outlier)$stats[c(1,5)] #(quantile(df[,1]-xmean,pnorm(c(-2,2)))*1.5)+xmean
      xlim=c(df$A[which(df$A.trans==xlim[1])[1]],df$A[which(df$A.trans==xlim[2])[1]])
    }
    if (is.null(ylim) && remove.outlier) {
      ylim=grDevices::boxplot.stats(df$B.trans[is.finite(df$B.trans)],coef=remove.outlier)$stats[c(1,5)] #(quantile(df[,2]-ymean,pnorm(c(-2,2)))*1.5)+ymean
      ylim=c(df$B[which(df$B.trans==ylim[1])[1]],df$B[which(df$B.trans==ylim[2])[1]])
    }
    if (is.null(xlim)) xlim=range(df$A[is.finite(df$A)])
    if (is.null(ylim)) ylim=range(df$B[is.finite(df$B)])
    clip=function(v,ch,lim,minus) ifelse(!is.finite(v),v,ifelse(ch<lim[1],minus,ifelse(ch>lim[2],Inf,v)))
    df$A.trans=clip(df$A.trans,df$A,xlim,-Inf)
    df$B.trans=clip(df$B.trans,df$B,ylim,-Inf)
    df$A=clip(df$A,df$A,xlim,if (log.x) 0 else -Inf)
    df$B=clip(df$B,df$B,ylim,if (log.y) 0 else -Inf)
  }
  else {
    xlim=range(df$A[!is.infinite(df$A)])
    ylim=range(df$B[!is.infinite(df$B)])
  }

  color=substitute(color)
  if (is.null(color)) {
    if (is.null(df$facet)) {
      df$color=density2d(df$A.trans, df$B.trans, n = density.n,margin = density.margin)
    } else {
      df$color=density2d(df$A.trans, df$B.trans, n = density.n,margin = density.margin,facet = df$facet)
    }
    if (is.null(colorpalette)) {
      colorscale=scale_color_viridis_c(name = "Density",guide="none")
    } else {
      col=make.continuous.colors(values=df$color,colors = colorpalette)
      colorscale = scale_color_gradientn(name="Density",guide="none",colors=col$colors)
    }
#  } else if (length(color)==1 && as.character(color) %in% names(GeneInfo(data))) {
#    df$color=GeneInfo(data,color)[ToIndex(data,genes)]
#    colorscale=scale_color_discrete(deparse(color))
  } else {
    ccol=if (length(color)==1 && as.character(color) %in% names(df)) df[[as.character(color)]] else eval(color,df,parent.frame())
    if (length(ccol)==1 && is.character(ccol) && ccol %in% names(df)) {
      ccol=df[[ccol]]
    }
    df$color=ccol
    if (is.character(df$color)) {
      colorscale=scale_color_identity(guide="none")
    } else if (is.factor(df$color)) {
      colorscale=scale_color_discrete(deparse(color))
    } else {
      col=make.continuous.colors(values = df$color,colors=colorpalette)
      df$color=pmin(pmax(df$color,min(col$breaks)),max(col$breaks))
      colorscale = scale_color_gradientn(deparse(color),colors=col$colors[c(1,3,5)],breaks=scales::breaks_pretty(n=5)(df$color),limits=col$breaks[c(1,5)])
    }
  }

  point.fun=function(...) geom_point(...)
  if ((is.null(rasterize) && nrow(df)>1000) || rlang::is_true(rasterize)) {
    checkPackages("ggrastr")
    point.fun=function(...) ggrastr::geom_point_rast(...)
  }

  g=ggplot(df,aes(A,B,color=color))+cowplot::theme_cowplot()
  if (!is.null(cross)) g=g+geom_vline(xintercept = 0,linetype=2,color='gray')+geom_hline(yintercept = 0,linetype=2,color='gray')
  if (!is.null(diag)) {
    if (is.logical(diag)) diag=0
      g=g+geom_abline(linetype=2,color='gray',intercept = diag)
  }

  if (!is.null(layers.below)) for (e in layers.below) g=g+e
  g=g+
    point.fun(size=size)+
    xlab(xlab)+ylab(ylab)

  if (!is.null(df$facet)) g=g+facet_wrap(~facet)

  if (!is.null(colorscale)) g=g+colorscale

  if (!is.null(highlight)) {
    if (is.list(highlight)){
      for (col in names(highlight)) {
        g=g+geom_point(data=df[highlight[[col]],],color=col,size=size*3)
      }
    } else {
      g=g+geom_point(data=df[highlight,],color='red',size=size*3)
    }
  }
  if (!is.null(label)) {
    df2=df[if (is.grandR(data)) ToIndex(data,label) else label,]
    df2$label=if (is.character(label)) label else rownames(df2)[label]
    if (checkPackages("ggrepel",error = FALSE,warn = FALSE)) {
      g=g+ggrepel::geom_label_repel(data=df2,mapping=aes(label=label),show.legend = FALSE,max.overlaps = Inf,min.segment.length = 0,force=label.repel)
    } else {
      singleMessage("Install the ggrepel package to have nicer labels!")
      g=g+ggplot2::geom_label(data=df2,mapping=aes(label=label),show.legend = FALSE)
    }
  }
  if (!is.null(correlation)) {
    if (correlation.x==-Inf) correlation.hjust=-0.05
    if (correlation.x==Inf) correlation.hjust=1.05
    if (correlation.y==-Inf) correlation.vjust=-0.05
    if (correlation.y==Inf) correlation.vjust=1.05

    if(log.x && correlation.x==-Inf) correlation.x=0
    if(log.y && correlation.y==-Inf) correlation.y=0

    if (is.null(df$facet)) {
      cor.format=correlation(df$A.trans.unclipped,df$B.trans.unclipped)
      g=g+annotate("text",x=correlation.x,y=correlation.y,label=cor.format,hjust=correlation.hjust,vjust=correlation.vjust)
    } else {
      cor.format = plyr::ddply(df,plyr::.(facet),function(sub) data.frame(label=correlation(sub$A.trans.unclipped,sub$B.trans.unclipped)))
      cor.format$correlation.x=correlation.x
      cor.format$correlation.y=correlation.y
      g=g+geom_text(data=cor.format,aes(x=correlation.x,y=correlation.y,label=label),inherit.aes=FALSE, parse=FALSE, hjust=correlation.hjust,vjust=correlation.vjust)
    }
  }
  if (set.coord) g=g+coord_cartesian(ylim=ylim,xlim=xlim)
  if (log.x) g=g+scale_x_log10()
  if (log.y) g=g+scale_y_log10()

  attr(g, 'df') <- adaptInf(df,xlim,ylim)
  g
}

#PlotToxicityTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=PsiLFC,hl.quantile=0.8) {
#	w=GetData(data,"count",conditions=w4sU,table=T)[,1]
#	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetData(data,"count",conditions=no4sU,table=T)[,1]
#	ntr=GetData(data,"ntr",conditions=w4sU,table=T)[,1]
#	use=!is.na(w+n+ntr)
#	w=w[use]
#	n=n[use]
#	ntr=ntr[use]
#
#	phl=comp.hl(ntr)
#	df=data.frame(lfc=LFC.fun(w,n),PHL=phl)[ntr<1,]
#	df=df[df$PHL<quantile(df$PHL[is.finite(df$PHL)],hl.quantile),]
#	ggplot(df,aes(PHL,lfc,color=density2d(PHL, lfc, n = 100)))+
#			scale_color_viridis_c(name = "Density",guide='none')+
#			geom_point(alpha=1)+
#			geom_hline(yintercept=0)+
##			geom_smooth(method="loess")+
#			xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
#			scale_x_continuous(breaks=c())+
#			coord_cartesian(ylim=ylim)
#}



PlotExpressionTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=lfc::PsiLFC,hl.quantile=0.8) {
  # R CMD check guard for non-standard evaluation
  M <- lfc <- NULL

  w=GetTable(data,type="count",columns=w4sU)[,1]
	nn=no4sU
	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetTable(data,type="count",columns=nn)[,1]
	use=!is.na(w+n)
	w=w[use]
	n=n[use]

	df=data.frame(lfc=LFC.fun(w,n),M=(log10(w+1)+log10(n+1))/2)
	ggplot(df,aes(M,lfc,color=density2d(M, lfc, n = 100)))+
	  cowplot::theme_cowplot()+
	  scale_color_viridis_c(name = "Density",guide='none')+
			geom_point(alpha=1)+
			geom_hline(yintercept=0)+
#			geom_smooth(method="loess")+
			xlab("Mean expression")+ylab("log FC 4sU/no4sU")+
			coord_cartesian(ylim=ylim)
}

#' Convenience function to make the same type of plot for multple analyses.
#' @param data the grandR object that contains the data to be plotted
#' @param plot.fun the plottinf function to apply
#' @param analyses the analyses to plot (default: all)
#' @param add additional ggplot (e.g., geoms) objects to add
#' @param ... passed further to plot.fun
#' @return ggplot objects
#' @export
#' @concept globalplot
PlotAnalyses=function(data,plot.fun,analyses=Analyses(data),add=NULL,...) {
  lapply(analyses,function(analysis) {
    re=plot.fun(data,analysis=analysis,...)
    if (!is.null(add)) for (e in if (is.list(add)) add else list(add)) re=re+e
    re
  })
}

#' Make a Vulcano plot
#'
#' Plot log2 fold changes against -log10 multiple testing adjusted P values
#'
#' @param data the grandR object that contains the data to be plotted
#' @param analysis the analysis to plot (default: first analysis)
#' @param p.cutoff p-value cutoff (default: 0.05)
#' @param lfc.cutoff log fold change cutoff (default: 1)
#' @param annotate.numbers if TRUE, label the number of genes
#' @param ... further parameters passed to \link{PlotScatter}
#' @return a ggplot object
#' @export
#' @concept globalplot
VulcanoPlot=function(data,analysis=Analyses(data)[1],p.cutoff=0.05,lfc.cutoff=1,
                     annotate.numbers=TRUE,...) {
  # R CMD check guard for non-standard evaluation
  Q <- NULL

  df=GetAnalysisTable(data,analyses=analysis,regex=FALSE,columns=c("LFC|Q"),gene.info = FALSE)
  names(df)=gsub(".*.Q","Q",gsub(".*.LFC","LFC",names(df)))
  g=PlotScatter(df,x=LFC,y=-log10(Q),remove.outlier = FALSE,...)+
    xlab(bquote(log[2]~FC))+
    ylab(bquote("-"~log[10]~FDR))+
    geom_hline(yintercept=-log10(p.cutoff),linetype=2)+
    geom_vline(xintercept=if (lfc.cutoff==0) 0 else c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(analysis)

  if (annotate.numbers) {
    if (lfc.cutoff!=0) {
      n=table(cut(df$LFC,breaks=c(-Inf,-lfc.cutoff,lfc.cutoff,Inf)),factor(df$Q>p.cutoff,levels=c("FALSE","TRUE")))
      g=g+annotate("label",x=c(-Inf,0,Inf,-Inf,0,Inf),y=c(Inf,Inf,Inf,-Inf,-Inf,-Inf),label=paste0("n=",as.numeric(n)),hjust=c(-0.1,0.5,1.1,-0.1,0.5,1.1),vjust=c(1.1,1.1,1.1,-0.1,-0.1,-0.1))
    } else {
      n=table(cut(df$LFC,breaks=c(-Inf,lfc.cutoff,Inf)),factor(df$Q>p.cutoff,levels=c("FALSE","TRUE")))
      g=g+annotate("label",x=c(-Inf,Inf,-Inf,Inf),y=c(Inf,Inf,-Inf,-Inf),label=paste0("n=",as.numeric(n)),hjust=c(-0.1,1.1,-0.1,1.1),vjust=c(1.1,1.1,-0.1,-0.1))
    }
  }
  g
}


#' Make an MA plot
#'
#' Plot average expression vs. log2 fold changes
#'
#' @param data the grandR object that contains the data to be plotted
#' @param analysis the analysis to plot (default: first analysis)
#' @param aest parameter to set visual attributes of the plot
#' @param p.cutoff p-value cutoff (default: 0.05)
#' @param lfc.cutoff log fold change cutoff (default: 1)
#' @param label.numbers if TRUE, label the number of genes
#' @param highlight highlight these genes; can be either numeric indices, gene names, gene symbols or a logical vector (see details)
#' @param label label these genes; can be either numeric indices, gene names, gene symbols or a logical vector (see details)
#' @param label.repel force to repel labels from points and each other (increase if labels overlap)
#' @return a ggplot object
#' @export
#' @concept globalplot
MAPlot=function(data,analysis=Analyses(data)[1],aest=aes(),p.cutoff=0.05,
                lfc.cutoff=1,
                label.numbers=TRUE,
                highlight=NULL,
                label=NULL,
                label.repel=1) {
  # R CMD check guard for non-standard evaluation
  M <- Q <- NULL

  df=GetAnalysisTable(data,analyses=analysis,regex=FALSE,columns=c("M|LFC|Q"),gene.info = FALSE)
  if (is.numeric(analysis)) analysis=Analyses(data)[analysis]
  names(df)=gsub(".*.Q","Q",gsub(".*.LFC","LFC",gsub(".*.M","M",names(df))))
  aes=utils::modifyList(aes(M+1,LFC,color=ifelse(Q<p.cutoff,"Sig.","NS")),aest)
  g=ggplot(df,mapping=aes)+
    cowplot::theme_cowplot()+
    geom_point(size=1)+
    scale_x_log10()+
    scale_color_manual(values=c(Sig.="black",NS="grey50"),guide='none')+
    ylab(bquote(log[2]~FC))+
    xlab("Total expression")+
    geom_hline(yintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(analysis)

  if (!is.null(highlight)) {
    if (is.list(highlight)){
      for (col in names(highlight)) {
        g=g+geom_point(data=df[ToIndex(data,highlight[[col]]),],color=col,size=1.5)
      }
    } else {
      g=g+geom_point(data=df[ToIndex(data,highlight),],color='red',size=1.5)
    }
  }
  if (!is.null(label)) {
    if (label=="auto") label=abs(df$LFC)>lfc.cutoff & df$Q<p.cutoff & !is.na(df$LFC) & !is.na(df$Q)
    df2=df
    df2$label=""
    df2[label,"label"]=rownames(df2)[label]
    g=g+ggrepel::geom_label_repel(data=df2,mapping=aes(label=label),show.legend = FALSE,max.overlaps = Inf,min.segment.length = 0,force=label.repel)
  }
  if (label.numbers) {
    n=c(sum(df$LFC>lfc.cutoff & df$Q<p.cutoff,na.rm = TRUE),sum(df$LFC< -lfc.cutoff & df$Q<p.cutoff,na.rm = TRUE))
    g=g+annotate("label",x=c(Inf,Inf),y=c(Inf,-Inf),label=paste0("n=",n),hjust=c(1.1,1.1),vjust=c(1.1,-0.1))
  }
  g
}

#' Plot the distribution of gene types
#' @param data the grandR object to get the data to be plotted from
#' @param mode.slot which mode and slot to use
#' @param relative show percentage values?
#'
#' @return a ggplot object
#'
#' @export
#' @concept globalplot
PlotTypeDistribution=function(data,mode.slot=DefaultSlot(data),relative=FALSE) {
  # R CMD check guard for non-standard evaluation
  value <- Type <- NULL

  df=GetTable(data,type=mode.slot)
	df=sapply(levels(data$gene.info$Type),function(type) colSums(df[ data$gene.info$Type==type,]))
	df=df[,colSums(df)>0]
	if (relative) {
		df=df/rowSums(df)*100
		mode.slot=sprintf("%s [%%]",mode.slot)
	}
	df=reshape2::melt(df,varnames=c("Condition","Type"))
	ggplot(df,aes(Condition,value,fill=Type))+cowplot::theme_cowplot()+
	  geom_bar(stat="Identity")+scale_fill_brewer(palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab(mode.slot)+xlab(NULL)
}


setup.default.aes=function(data,aest) {
  # R CMD check guard for non-standard evaluation
  Replicate <- NULL

  if (is.null(aest)) aest=aes()
  if (!is.null(Condition(data)) && !any(c("color","colour") %in% names(aest))) aest=utils::modifyList(ggplot2::aes(color=Condition),aest)
  if (!is.null(Coldata(data)$Replicate)&& !any(c("shape") %in% names(aest))) aest=utils::modifyList(ggplot2::aes(shape=Replicate),aest)
  aest
}

#' Gene plot comparing old vs new RNA
#'
#' Plot the old vs new RNA values of a gene
#'
#' @param data the grandR object to get the data to be plotted from
#' @param gene the gene to plot
#' @param slot the slot of the grandR object to get the data from
#' @param columns which columns (i.e. samples or cells) to show (see details)
#' @param log show both axes in log scale
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param aest parameter to set the visual attributes of the plot
#' @param size the point size used for plotting; overridden if size is defined via aest
#'
#' @details The value of the aest parameter must be an \emph{Aesthetic mapping} as generated by \link[ggplot2]{aes}.
#'
#' @details The table used for plotting is the table returned by \link{GetData} with coldata set to TRUE, i.e. you can use all names from the \link{Coldata} table for aest.
#'
#' @details By default, aest is set to aes(color=Condition,shape=Replicate) (if both Condition and Replicate are names in the Coldata table).
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a ggplot object.
#'
#' @seealso \link{GetData}, \link{PlotGeneTotalVsNtr},\link{PlotGeneGroupsPoints},\link{PlotGeneGroupsBars}
#'
#' @export
#' @concept geneplot
PlotGeneOldVsNew=function(data,gene,slot=DefaultSlot(data),columns=NULL,log=TRUE,show.CI=FALSE,
                          aest=NULL,size=2) {
  # R CMD check guard for non-standard evaluation
  lower <- upper <- NULL

  aest=setup.default.aes(data,aest)
  if (length(slot)!=1) stop("Provide a single slot name!")
  slot=get.mode.slot(data,slot,allow.ntr = FALSE)$slot
  new=paste0("new.",slot)
  old=paste0("old.",slot)

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  df=GetData(data,genes=gene,mode.slot=c(old,new),columns=columns,by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)
  g=ggplot(df,utils::modifyList(aes(!!sym(old),!!sym(new)),aest))+cowplot::theme_cowplot()
  if (!is.null(aest$size)) g=g+geom_point() else g=g+geom_point(size=size)
  if (log) {
    g=g+scale_x_log10()+
      scale_y_log10()
  } else {
    g=g+scale_x_continuous()+
      scale_y_continuous()
  }
  g=g+xlab(paste0("Old RNA (",slot,")"))
  g=g+ylab(paste0("New RNA (",slot,")"))

  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper",slot),by.rows=FALSE,coldata=T,ntr.na = FALSE)[,c("lower","upper",slot)])
    g=g+geom_errorbar(data=df,mapping=aes(ymin=lower*!!sym(slot),ymax=upper*!!sym(slot)))
    g=g+geom_errorbarh(data=df,mapping=aes(xmin=(1-upper)*!!sym(slot),xmax=(1-lower)*!!sym(slot)))
  }
  g
}

#' Gene plot comparing total RNA vs the NTR
#'
#' Plot the total RNA expression vs the new-to-total RNA ratio for a gene
#'
#' @param data the grandR object to get the data to be plotted from
#' @param gene the gene to plot
#' @param slot the slot of the grandR object to get the data from
#' @param columns which columns (i.e. samples or cells) to show (see details)
#' @param log show the x axis (total RNA) in log scale
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param aest parameter to set the visual attributes of the plot
#' @param size the point size used for plotting; overridden if size is defined via aest
#'
#' @details The value of the aest parameter must be an \emph{Aesthetic mapping} as generated by \link[ggplot2]{aes}.
#'
#' @details The table used for plotting is the table returned by \link{GetData} with coldata set to TRUE, i.e. you can use all names from the \link{Coldata} table for aest.
#'
#' @details By default, aest is set to aes(color=Condition,shape=Replicate) (if both Condition and Replicate are names in the Coldata table).
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a ggplot object.
#'
#' @seealso \link{GetData}, \link{PlotGeneOldVsNew},\link{PlotGeneGroupsPoints},\link{PlotGeneGroupsBars}
#'
#' @export
#' @concept geneplot
PlotGeneTotalVsNtr=function(data,gene,slot=DefaultSlot(data),columns=NULL,log=TRUE,show.CI=FALSE,
                            aest=NULL,size=2) {
  # R CMD check guard for non-standard evaluation
  lower <- upper <- ntr <- NULL

  aest=setup.default.aes(data,aest)
  if (length(slot)!=1) stop("Provide a single slot name!")
  slot=get.mode.slot(data,slot,allow.ntr = FALSE)$slot

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  df=GetData(data,genes=gene,mode.slot=c("ntr",slot),columns=columns,by.rows=FALSE,coldata=T,ntr.na = FALSE)
  g=ggplot(df,utils::modifyList(aes(!!sym(slot),ntr),aest))+cowplot::theme_cowplot()
  if (!is.null(aest$size)) g=g+geom_point() else g=g+geom_point(size=size)
  if (log) {
    g=g+scale_x_log10()
  } else {
    g=g+scale_x_continuous()
  }
  g=g+scale_y_continuous("NTR")
  g=g+xlab(paste0("Total RNA (",slot,")"))

  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper"),by.rows=FALSE,coldata=T,ntr.na = FALSE)[,c("lower","upper")])
    g=g+geom_errorbar(data=df,mapping=aes(ymin=lower,ymax=upper))
  }
  g
}

#' Plot gene groups as points
#'
#' Plot either old, new or total RNA of a gene in a row, per condition.
#'
#' @param data the grandR object to get the data to be plotted from
#' @param gene the gene to plot
#' @param group how to group the genes (default: Condition)
#' @param mode.slot the  mode.slot of the grandR object to get the data from
#' @param columns which columns (i.e. samples or cells) to show (see details)
#' @param log show the y axis in log scale
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param aest parameter to set the visual attributes of the plot
#' @param size the point size used for plotting; overridden if size is defined via aest
#' @param transform function that is called on the data frame directly before plotting (can be NULL)
#'
#' @details The value of the aest parameter must be an \emph{Aesthetic mapping} as generated by \link[ggplot2]{aes}.
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a
#' dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @details The table used for plotting is the table returned by \link{GetData} with coldata set to TRUE, i.e. you can use all names from the \link{Coldata} table for aest.
#'
#' @details By default, aest is set to aes(color=Condition,shape=Replicate) (if both Condition and Replicate are names in the Coldata table).
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a ggplot object.
#'
#' @seealso \link{GetData}, \link{PlotGeneTotalVsNtr},\link{PlotGeneOldVsNew},\link{PlotGeneGroupsBars}
#'
#' @export
#' @concept geneplot
PlotGeneGroupsPoints=function(data,gene,group="Condition",mode.slot=DefaultSlot(data),
                              columns=NULL,
                              log=TRUE,
                              show.CI=FALSE,
                              aest=NULL,size=2,transform=NULL) {
  # R CMD check guard for non-standard evaluation
  lower <- upper <- value <- NULL

  aest=setup.default.aes(data,aest)

  if (length(group)!=1 && !group %in% names(Coldata(data))) stop("Group must be a name in the Coldata table!")

  if (length(mode.slot)!=1) stop("Provide a single slot name!")
  mode.slot=get.mode.slot(data,mode.slot)
  slot=mode.slot$slot
  mode=mode.slot$mode

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  if (slot=="ntr") {
    df=GetData(data,genes=gene,mode.slot="ntr",columns=columns,by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)
    df$ntr=df$Value
    df$value=df[[slot]]
    log=FALSE
  } else {
    df=GetData(data,genes=gene,mode.slot=c(slot,"ntr"),columns=columns,by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)
    df$value=switch(mode[1],total=df[[slot]],new=df[[slot]]*df[["ntr"]],old=df[[slot]]*(1-df[["ntr"]]),stop(paste0(mode," unknown!")))
  }

  if (!is.null(transform)) df=transform(df)

  g=ggplot(df,utils::modifyList(aes(!!sym(group),value),aest))+cowplot::theme_cowplot()
  if (!is.null(aest$size)) g=g+geom_point(position=if(show.CI) position_dodge(width=0.4) else "identity") else g=g+geom_point(size=size,position=if(show.CI) position_dodge(width=0.4) else "identity")
  g=g+xlab(NULL)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (log) g=g+scale_y_log10()

  if (slot=="ntr")g=g+ylab("NTR") else g=g+ylab(paste0(toupper(substr(mode,1,1)),substr(mode,2,nchar(mode))," RNA (",slot,")"))

  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper"),by.rows=FALSE,coldata=FALSE,ntr.na = FALSE)[,c("lower","upper")])
    if (slot=="ntr") {
      g=g+geom_errorbar(data=df,mapping=aes(ymin=lower,ymax=upper),width=0,position=position_dodge(width=0.4))
    } else {
      g=switch(mode[1],
               total=g,
               new=g+geom_errorbar(data=df,mapping=aes(ymin=lower*!!sym(slot),ymax=upper*!!sym(slot)),width=0,position=position_dodge(width=0.4)),
               old=g+geom_errorbar(data=df,mapping=aes(ymin=(1-upper)*!!sym(slot),ymax=(1-lower)*!!sym(slot)),width=0,position=position_dodge(width=0.4)),
               stop(paste0(mode," unknown!")))
    }
  }
  g
}

#' Plot gene values as bars
#'
#' Plot old and new RNA of a gene in a row.
#'
#' @param data the grandR object to get the data to be plotted from
#' @param gene the gene to plot
#' @param slot the slot of the grandR object to get the data from
#' @param columns which columns (i.e. samples or cells) to show (see details)
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param xlab The names to show at the x axis;
#' @param transform function that is called on the data frame directly before plotting (can be NULL)
#'
#' @details xlab can be given as a character vector or an expression that evaluates into a character vector.
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently it.
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a ggplot object.
#'
#' @seealso \link{GetData}, \link{PlotGeneTotalVsNtr},\link{PlotGeneOldVsNew},\link{PlotGeneGroupsBars}
#'
#' @export
#' @concept geneplot
PlotGeneGroupsBars=function(data,gene,slot=DefaultSlot(data),columns=NULL,show.CI=FALSE,xlab=NULL,transform=NULL) {
  # R CMD check guard for non-standard evaluation
  Name <- Value <- mode.slot <- upper <- lower <- NULL

  if (length(slot)!=1) stop("Provide a single slot name!")
  slot=get.mode.slot(data,slot,allow.ntr = FALSE)$slot

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  df=GetData(data,genes=gene,mode.slot=paste0(c("new.","old."),slot),columns=columns,by.rows=TRUE,coldata=TRUE,ntr.na = FALSE)

  xlab=substitute(xlab)
  xlab=if (is.null(xlab)) df$Name else eval(xlab,df,parent.frame())
  df$xlab=xlab

  if (!is.null(transform)) df=transform(df)

  g=ggplot(df,aes(Name,Value,fill=mode.slot))+
    cowplot::theme_cowplot()+
    geom_bar(stat="identity",position=position_stack())+
    scale_fill_manual(NULL,values = c('red','gray'),guide="none")+
    xlab(NULL)+
    scale_x_discrete(breaks=df$Name,labels=df$xlab)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  g=g+ylab(paste0("Total RNA (",slot,")"))
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df2=GetData(data,genes=gene,mode.slot=c("lower","upper",slot),by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)
    g=g+geom_errorbar(data=df2,mapping=aes(y=!!sym(slot),fill=NULL,ymin=(1-upper)*!!sym(slot),ymax=(1-lower)*!!sym(slot)),width=0)
  }
  g
}

#' Gene plot for snapshot timecourse data
#'
#' Plot the total RNA expression vs the new-to-total RNA ratio for a gene
#'
#' @param data the grandR object to get the data to be plotted from
#' @param gene the gene to plot
#' @param time the times to show on the x axis (see details)
#' @param mode.slot the mode.slot of the grandR object to get the data from
#' @param columns which columns (i.e. samples or cells) to show (see details)
#' @param average.lines add average lines?
#' @param exact.tics use axis labels directly corresponding to the available temporal values?
#' @param log show the y axis in log scale
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param aest parameter to set the visual attributes of the plot
#' @param size the point size used for plotting; overridden if size is defined via aest
#'
#' @details The x axis of this plot will show a temporal dimension. The time parameter defines a name in the \link{Coldata} table containing the temporal values for each sample.
#'
#' @details The value of the aest parameter must be an \emph{Aesthetic mapping} as generated by \link[ggplot2]{aes}.
#'
#' @details The table used for plotting is the table returned by \link{GetData} with coldata set to TRUE, i.e. you can use all names from the \link{Coldata} table for aest.
#'
#' @details By default, aest is set to aes(color=Condition,shape=Replicate) (if both Condition and Replicate are names in the Coldata table).
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a ggplot object.
#'
#' @seealso \link{GetData}, \link{PlotGeneOldVsNew},\link{PlotGeneGroupsPoints},\link{PlotGeneGroupsBars}
#'
#' @export
#' @concept geneplot
PlotGeneSnapshotTimecourse=function(data,gene,time=Design$dur.4sU,
                            mode.slot=DefaultSlot(data),
                            columns=NULL,
                            average.lines=TRUE,
                            exact.tics=TRUE,
                            log=TRUE,
                            show.CI=FALSE,aest=NULL,size=2) {
  # R CMD check guard for non-standard evaluation
  x <- colour <- group <- Value <- ntr <- lower <- upper <- linetype <- NULL

  aest=setup.default.aes(data,aest)
  if (length(slot)!=1) stop("Provide a single slot name!")

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  df=GetData(data,genes=gene,mode.slot=mode.slot,columns=columns,by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)

  aes=utils::modifyList(aes(x=!!sym(time),y=Value),aest)
  if (!is.null(Condition(data))) aes=utils::modifyList(aes,aes(group=Condition))
  breaks=if (exact.tics) sort(unique(df[[time]])) else scales::breaks_extended(5)(df[[time]])

  oslot=mode.slot
  mode.slot=get.mode.slot(data,mode.slot)
  slot=mode.slot$slot
  mode=mode.slot$mode
  if (slot=="ntr") log=FALSE

  g=ggplot(df,mapping=aes)+cowplot::theme_cowplot()
  if (!is.null(aest$size)) g=g+geom_point(position=if(show.CI) position_dodge(width=0.4) else "identity") else g=g+geom_point(size=size,position=if(show.CI) position_dodge(width=0.4) else "identity")
  g=g+scale_x_continuous(NULL,labels = scales::number_format(accuracy = max(0.01,my.precision(breaks)),suffix="h"),breaks=breaks)
  if (slot=="ntr")g=g+ylab("NTR") else g=g+ylab(paste0(toupper(substr(mode,1,1)),substr(mode,2,nchar(mode))," RNA (",slot,")"))
  if (log) g=g+scale_y_log10()
  if (average.lines) {
    # compute average line:
    if (!is.null(Condition(data))) {
      ddf=as.data.frame(lapply(aes,function(col) rlang::eval_tidy(col,data=df)))
      ddf=plyr::ddply(ddf,intersect(names(ddf),c("x","colour","linetype","group")),function(s) c(Value=mean(s$y,na.rm=TRUE)))
      daes=aes(x,Value,group=group)
      if ("colour" %in% names(ddf)) daes=utils::modifyList(daes,aes(color=colour))
      if ("linetype" %in% names(ddf)) daes=utils::modifyList(daes,aes(linetype=linetype))
      g=g+geom_line(data=ddf,mapping=daes,inherit.aes=FALSE)
    } else {
      ddf=as.data.frame(lapply(aes,function(col) rlang::eval_tidy(col,data=df)))
      ddf=plyr::ddply(ddf,plyr::.(x),function(s) c(Value=mean(s$y,na.rm=TRUE)))
      g=g+geom_line(data=ddf,mapping=aes(x,Value),inherit.aes=F)
    }
  }
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df2=GetData(data,genes=gene,mode.slot=unique(c("lower","upper",slot,oslot)),by.rows=FALSE,coldata=TRUE,ntr.na = FALSE)
    if (slot=="ntr") {
      g=g+geom_errorbar(data=df2,mapping=aes(y=ntr,ymin=lower,ymax=upper),width=0,position=position_dodge(width=0.4))
    } else {
      g=switch(mode[1],
               total=g,
               new=g+geom_errorbar(data=df2,mapping=aes(y=!!sym(oslot),ymin=lower*!!sym(slot),ymax=upper*!!sym(slot)),width=0,position=position_dodge(width=0.4)),
               old=g+geom_errorbar(data=df2,mapping=aes(y=oslot,ymin=(1-upper)*!!sym(slot),ymax=(1-lower)*!!sym(slot)),width=0,position=position_dodge(width=0.4)),
               stop(paste0(mode," unknown!")))
    }
  }
  g
}

#' Rotate x axis labels
#'
#' Add this to a ggplot object to rotate the x axis labels
#'
#' @param angle the angle by which to rotate
#'
#' @return a ggplot theme object
#' @export
#'
#' @concept helper
RotatateAxisLabels=function(angle=90) {
  theme(axis.text.x = element_text(angle = angle, vjust = if (abs(angle-90)<10) 0.5 else 1, hjust=1))
}



