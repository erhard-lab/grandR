

Simulate=function(s=100*d,d=log(2)/hl,hl=2,l0=s/d,min.time=-1,max.time=10,length.out = 1000,by = ((max.time - 0)/(length.out - 1)),name=NULL) {
    ode.new=function(t,s,d) ifelse(t<0,0,s/d*(1-exp(-t*d)))
    ode.old=function(t,f0,s,d) ifelse(t<0,s/d,f0*exp(-t*d))
    t=seq(min.time,max.time,by=by)
    old=ode.old(t,l0,s,d)
    new=ode.new(t,s,d)
    re=data.frame(
        Time=t,
        Value=c(old,new,old+new,new/(old+new)),
        Type=factor(rep(c("Old","New","Total","NTR"),each=length(t)),levels=c("Old","New","Total","NTR"))
    )
    if (!is.null(name)) re$Name=name
    re
}


PlotSimulation=function(sim.df,ntr=TRUE,old=TRUE,new=TRUE,total=TRUE) {
    if (!ntr) sim.df=sim.df[sim.df$Type!="NTR",]
    if (!old) sim.df=sim.df[sim.df$Type!="Old",]
    if (!new) sim.df=sim.df[sim.df$Type!="New",]
    if (!total) sim.df=sim.df[sim.df$Type!="Total",]
    ggplot(sim.df,aes(Time,Value,color=Type))+
        geom_line(size=1)+
        scale_color_manual(NULL,values=c(Old="#54668d",New="#953f36",Total="#373737",NTR="#e4c534"))+
        facet_wrap(~ifelse(Type=="NTR","NTR","Timecourse"),scales="free_y",ncol=1)+
        ylab(NULL)+
        scale_x_continuous(breaks=scales::pretty_breaks())+
        theme(
	  strip.background = element_blank(),
	  strip.text.x = element_blank()
	)
}

PlotCompareNTRs=function(...) {
    dfs=list(...)
    df=do.call("rbind",dfs)
    ggplot(df[df$Type=="NTR",],aes(Time,Value,color=Name))+
        geom_line(size=1)
}

#sim.hl8.steady=Simulate(hl=8,name="Steady-state, 8h")
#sim.hl4.steady=Simulate(hl=4,name="Steady-state, 4h")
#sim.hl2.steady=Simulate(hl=2,name="Steady-state, 2h")
#sim.hl1.steady=Simulate(hl=1,name="Steady-state, 1h")
#sim.hl8.2x=Simulate(l0=10,hl=4,s=10*log(2)/8,name="2x down, 8h")
#sim.hl4.2x=Simulate(l0=10,hl=2,s=10*log(2)/4,name="2x down, 4h")
#sim.hl2.2x=Simulate(l0=10,hl=1,s=10*log(2)/2,name="2x down, 2h")
#sim.hl1.2x=Simulate(l0=10,hl=0.5,s=10*log(2)/1,name="2x down, 1h")
#sim.hl8.10x=Simulate(l0=10,hl=0.8,s=10*log(2)/8,name="10x down, 8h")
#sim.hl4.10x=Simulate(l0=10,hl=0.4,s=10*log(2)/4,name="10x down, 4h")
#sim.hl2.10x=Simulate(l0=10,hl=0.2,s=10*log(2)/2,name="10x down, 2h")
#sim.hl1.10x=Simulate(l0=10,hl=0.1,s=10*log(2)/1,name="10x down, 1h")



#pdf("destabilized.pdf",width=10,height=10)
#plot_grid(
#    PlotSimulation(sim.hl8.steady),
#    PlotSimulation(sim.hl8.2x),
#    PlotSimulation(sim.hl8.10x),
#    PlotCompareNTRs(sim.hl8.steady,sim.hl8.2x,sim.hl8.10x),
#    ncol=2
#)
    

#plot_grid(
#    PlotSimulation(sim.hl4.steady),
#    PlotSimulation(sim.hl4.2x),
#    PlotSimulation(sim.hl4.10x),
#    PlotCompareNTRs(sim.hl4.steady,sim.hl4.2x,sim.hl4.10x),
#    ncol=2
#)

#plot_grid(
#    PlotSimulation(sim.hl2.steady),
#    PlotSimulation(sim.hl2.2x),
#    PlotSimulation(sim.hl2.10x),
#    PlotCompareNTRs(sim.hl2.steady,sim.hl2.2x,sim.hl2.10x),
#    ncol=2
#)

#plot_grid(
#    PlotSimulation(sim.hl1.steady),
#    PlotSimulation(sim.hl1.2x),
#    PlotSimulation(sim.hl1.10x),
#    PlotCompareNTRs(sim.hl1.steady,sim.hl1.2x,sim.hl1.10x),
#    ncol=2
#)

#dev.off()
