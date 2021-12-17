
read.tsv=function(t,...) {
	t=read.delim(t,stringsAsFactors=FALSE,check.names=FALSE,...)
	t=as.data.frame(lapply(t,function(c) if (is.character(c)) factor(c,levels=unique(c)) else c),check.names=FALSE)
	t
}

combine=function(...,sep=".") {
	l=list(...)
	re=l[[1]]
	if (length(l)>1) for (i in 2:length(l)) {
		re=paste(rep(re,each=length(l[[i]])),rep(l[[i]],length(re)),sep=sep)
	}
	re
}

