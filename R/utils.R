
read.tsv=function(t) {
	t=read.delim(t,stringsAsFactors=FALSE,check.names=FALSE)
	t=as.data.frame(lapply(t,function(c) if (is.character(c)) factor(c,levels=unique(c)) else c),check.names=FALSE)
	t
}

