
ServeData=function(data,...,df=GetDiffExpTable(data,cols=c("LFC","Q")),sizes=12,height=400,plots=NULL,title=Title(data),show.sessionInfo=FALSE,help=list(".Q: multiple testing corrected p values",".LFC: log2 fold changes") ) {

  
  plot.funs=list(...)
  
  if (length(plot.funs)==0) plot.funs=list(PlotGeneOldVsNew)
  if (length(sizes)==1) sizes=rep(floor(12/min(4,length(plot.funs))),min(4,length(plot.funs)))
  if (length(sizes)!=length(plot.funs)) stop("sizes need to be length 1 or same length as plots!")
  sizes=c(sizes,rep(1,8))
  
  if (!is.null(plots)) {
    plots=lapply(plots, function(p) if (is.function(p)) p(data) else p)
  }
  
  if (!is.null(help) && is.list(help)) help=sprintf("<span style='padding-top:25px;'><span class='help-block well'>Table columns:%s</span></span>", paste(sapply(help,function(s) sprintf("<li><span>%s</span></li>",s)),collapse="\n"))
  
	server=function(input, output,session) {
	  options(DT.options = list(pageLength = 12))
	  output$tab <- DT::renderDataTable(DT::datatable(df, selection = 'single',rownames = FALSE, escape=-1,filter = "top",extensions = 'Buttons', options = list(
	    dom = 'Bfrtip',
	    buttons = c('copy', 'csv', 'excel')
	  )) %>%formatRound(names(df)[grepl("\\.LFC$",names(df))], 2)%>%formatSignif(names(df)[grepl("\\.Q$",names(df))], 2))

	  output$plot1=renderPlot({ if (length(input$tab_rows_selected)==1) plot.funs[[1]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot2=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=2) plot.funs[[2]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot3=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=3) plot.funs[[3]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot4=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=4) plot.funs[[4]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot5=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=5) plot.funs[[5]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot6=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=6) plot.funs[[6]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot7=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=7) plot.funs[[7]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$plot8=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.funs)>=8) plot.funs[[8]](data=data,gene=df$Symbol[input$tab_rows_selected])  })
	  output$helpText=renderText({ if (length(input$tab_rows_selected)==0 && !is.null(help)) help  })
	  
	  if (!is.null(plots)) {
	    for (n in names(plots)) {
	      create=function(n) {
  	      env=new.env()
  	      env$n=n
  	      getwidth=function() {
  	        w=attr(plots[[n]][[input[[paste0(n,"list")]]]],"width")
  	        if (is.null(w)) w=7*100
  	        w
  	      }
  	      getheight=function() {
  	        w=attr(plots[[n]][[input[[paste0(n,"list")]]]],"height")
  	        if (is.null(w)) w=7*100
  	        w
  	      }
  	      renderPlot({plots[[n]][[input[[paste0(n,"list")]]]](data)},width=getwidth,height=getheight,env=env)
	      }
	      output[[paste0(n,"plot")]]=create(n)
	    }
	  }
	  
	  if (show.sessionInfo) output$sessionInfo <- renderPrint({
	    capture.output(sessionInfo())
	  })
	  
	  
	}


	plot.ui=NULL
	if (!is.null(plots)) {
	  
	  plist=c(lapply(names(plots),function(n) tabPanel(n,
	                                                   selectInput(paste0(n,"list"),n,names(plots[[n]]),selectize=FALSE,size=10),
	                                                   plotOutput(paste0(n,"plot"))
	                                                   )),list(title="Plots"))
	  
	  plot.ui=do.call("navbarMenu",plist)
	}
	
	more=NULL
	if (show.sessionInfo)
	  more=navbarMenu("More",
	                  tabPanel("Info",verbatimTextOutput("sessionInfo"))
	  )
	
	
	ui=list(
        	  tabPanel("Data",
        	  fluidPage(
        	  fluidRow(
        	    column(12, DT::dataTableOutput('tab'))
        	  ),
        	  conditionalPanel(
        	    condition = "helpText",
        	    fluidRow(column(10, htmlOutput("helpText")))
        	  ),
        	  fluidRow(
        	    column(sizes[1], plotOutput("plot1",height = height)),
        	    conditionalPanel(
        	      condition = "plot2",
        	      column(sizes[2], plotOutput("plot2",height = height))
        	    ),
        	    conditionalPanel(
        	      condition = "plot3",
        	      column(sizes[3], plotOutput("plot3",height = height))
        	    ),
        	    conditionalPanel(
        	      condition = "plot4",
        	      column(sizes[4], plotOutput("plot4",height = height))
        	    )
        	  ),
        	  fluidRow(
        	    conditionalPanel(
        	      condition = "plot5",
        	      column(sizes[5], plotOutput("plot5",height = height))
        	    ),
        	    conditionalPanel(
        	      condition = "plot6",
        	      column(sizes[6], plotOutput("plot6",height = height))
        	    ),
        	    conditionalPanel(
        	      condition = "plot7",
        	      column(sizes[7], plotOutput("plot7",height = height))
        	    ),
        	    conditionalPanel(
        	      condition = "plot8",
        	      column(sizes[8], plotOutput("plot8",height = height))
        	    )
        	  )
        	  #do.call("fluidRow",lapply(1:length(plot.funs),function(i) column(sizes[i],plotOutput(paste0("plot.funs",i),height=height))))
        	)),
        	
        	plot.ui,
        	
        	more,
        	  
        	tags$script(HTML(sprintf("
        	var header = $('.navbar> .container-fluid');
          header.append('<div class=\"nav navbar-nav\" style=\"float:right\"><span class=\"navbar-brand\">%s</span></div>')",
        	                         VersionString(data)
        	))),
        	
        	tags$script(HTML(
        	  "$(document).ready(function(){
              
              $('[data-toggle=tab]').on('click', function(e){
                if (($(this).attr('data-value'))=='Download table') {
                  $('#download')[0].click()
                  $('.dropdown').removeClass('open');
                  e.stopPropagation();
                  e.preventDefault();
                }
              });
            });"
        	))
	)
	
	ui=ui[!sapply(ui,is.null)]
	myui=function(...) navbarPage(title,...)
	ui=do.call("myui",ui)
	
	shinyApp(ui = ui, server = server)
}



ServeData.legacy=function(data,df=GetDiffExpTable(data,cols=c("LFC","Q")),aest=aes(color=Sample),aest.ts=aes(color=Sample),time="hpi",ggparam=NULL,ggparam.ts=NULL,average.lines=FALSE, plot2=NA, plot3=NA) {
  
  server=function(input, output,session) {
    options(DT.options = list(pageLength = 12))
    output$tab <- DT::renderDataTable(DT::datatable(df, selection = 'single',rownames = FALSE, escape=-1)) # %>%formatRound(names(df)[grepl("q$",names(df))], 2)
    
    output$plot1=renderPlot({
      if (length(input$tab_rows_selected)==1) {
        gene=df$Symbol[input$tab_rows_selected]
        g=ggplot(GetData(data,gene=gene,type=c("tpm","ntr"),melt=F,coldata=T),modifyList(aes(tpm,ntr),aest))+geom_point(size=3)+scale_x_log10()
        if (!is.null(ggparam)) g=g+ggparam
        g
      }
    })
    
    output$plot2=renderPlot({
      if (length(input$tab_rows_selected)==1) {
        gene=df$Symbol[input$tab_rows_selected]
        if (is.null(plot2)) {
          NULL			
        } else if (is.na(plot2)){
          df=GetData(data,gene=gene,type=c("tpm"),melt=T,coldata=T)
          aes=modifyList(aes_string(time,"Value"),aest.ts)
          g=ggplot(df,mapping=aes)+geom_point()+scale_y_log10()+xlab(NULL)+ylab("Total TPM")
          if (!is.null(ggparam.ts)) g=g+ggparam.ts
          if (average.lines) {
            # compute average line:
            ddf=as.data.frame(lapply(aes,function(col) rlang::eval_tidy(col,data=df)))
            ddf=ddply(ddf,.(x,colour),function(s) c(Value=mean(s$y,na.rm=TRUE)))
            g=g+geom_line(data=ddf,mapping=aes(x,Value,colour=colour,group=colour),inherit.aes=F)
          }
          g
        } else {
          plot2(gene)
        }
      }
    })
    
    output$plot3=renderPlot({
      if (length(input$tab_rows_selected)==1) {
        gene=df$Symbol[input$tab_rows_selected]
        if (is.null(plot3)) {
          NULL			
        } else if (is.na(plot3)){
          df=GetData(data,gene=gene,type=c("new.tpm"),melt=T,coldata=T)
          aes=modifyList(aes_string(time,"Value"),aest.ts)
          g=ggplot(df,mapping=aes)+geom_point()+scale_y_log10()+xlab(NULL)+ylab("New TPM")
          if (!is.null(ggparam.ts)) g=g+ggparam.ts
          if (average.lines) {
            # compute average line:
            ddf=as.data.frame(lapply(aes,function(col) rlang::eval_tidy(col,data=df)))
            ddf=ddply(ddf,.(x,colour),function(s) c(Value=mean(s$y,na.rm=TRUE)))
            g=g+geom_line(data=ddf,mapping=aes(x,Value,colour=colour,group=colour),inherit.aes=F)
          }
          g
        } else {
          plot3(gene)
        }
      }
    })
    
  }
  
  
  ui=fluidPage(
    fluidRow(
      column(12, DT::dataTableOutput('tab'))
    ),
    fluidRow(
      column(4, plotOutput("plot1",height = 400)),
      column(if (is.null(plot3)) 8 else 4, plotOutput("plot2",height = 400)),
      column(4, plotOutput("plot3",height = 400))
    )
  )
  
  shinyApp(ui = ui, server = server)
}


