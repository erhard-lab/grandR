
ServeData=function(data,
                   df=GetAnalysisTable(data,columns="Synthesis|Half-life|LFC|Q"),
                   df.set=df,
                   sizes=NA,height=400,
                   plot.single=list(), plot.set=list(), plot.static=list(),
                   df.identifier="Symbol",
                   title=Title(data),
                   show.sessionInfo=FALSE,
                   help=list(".Q: multiple testing corrected p values",".LFC: log2 fold changes") ) {


  if (length(plot.single)==0) plot.single=list(PlotGeneOldVsNew)
  if (length(sizes)==1 & is.na(sizes)) sizes=rep(floor(12/min(4,length(plot.single))),length(plot.single))
  if (length(sizes)!=length(plot.single)) stop("sizes need to be length 1 or same length as plots!")
  sizes=c(sizes,rep(1,8))

#    plot.static=lapply(plot.static, function(p) if (is.function(p)) p(data) else p)

  if (!is.null(help) && is.list(help)) help=sprintf("<span style='padding-top:25px;'><span class='help-block well'>Table columns:%s</span></span>", paste(sapply(help,function(s) sprintf("<li><span>%s</span></li>",s)),collapse="\n"))
	server=function(input, output,session) {
	  dttab=DT::datatable(df,
	                      callback = DT::JS("$('div#buttons').css('float','left').css('margin-right','50px'); $('div#clip').css('float','left'); $('div#buttons').append($('#download1')); $('div#buttons').append($('#clip')); "),
	                      selection = 'single',
	                      rownames = FALSE,
	                      escape=-1,
	                      filter = "top",
	                      options = list(
	                        pageLength = 10,
	                        dom = '<"#buttons">lfrtip'
	                      ))
	  dttab=DT::formatRound(dttab,names(df)[sapply(df,class)=="numeric"], 2)
	  dttab=DT::formatRound(dttab,names(df)[grepl("\\.LFC$",names(df))], 2)
	  dttab=DT::formatSignif(dttab,names(df)[grepl("\\.Q$",names(df))], 2)
	  dttab=DT::formatSignif(dttab,names(df)[grepl("\\.P$",names(df))], 2)
	  output$tab <- DT::renderDataTable(dttab)
	  output$download1 <- downloadHandler(
	    filename = function() {
	      paste0(title,"-", Sys.Date(), ".tsv")
	    },
	    content = function(file) {
	      write.table(df[input$tab_rows_all,], file,row.names=F,col.names=T,quote=F,sep="\t")
	    }
	  )
	  output$clip <- renderUI({
	    nn=if(.row_names_info(df)<0) df[input$tab_rows_all,1] else rownames(df)[input$tab_rows_all]
	    rclipboard::rclipButton("clipbtn", "Copy", paste(nn,collapse="\n"), icon("clipboard"))
	  })
	  observeEvent(input$clipbtn, {showNotification(
	    sprintf("Copied %d names",length(input$tab_rows_all)),
	    duration = 2,
	    type = "message"
	  )})


	  output$plot1=renderPlot({ if (length(input$tab_rows_selected)==1) plot.single[[1]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot2=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=2) plot.single[[2]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot3=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=3) plot.single[[3]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot4=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=4) plot.single[[4]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot5=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=5) plot.single[[5]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot6=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=6) plot.single[[6]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot7=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=7) plot.single[[7]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot8=renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.single)>=8) plot.single[[8]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$helpText=renderText({ if (length(input$tab_rows_selected)==0 && !is.null(help)) help  })

	  for (n in names(plot.static)) {
	    create=function(n) {
	      env=new.env()
	      env$n=n
	      getwidth=function() {
	        w=attr(plot.static[[n]][[input[[paste0(n,"list")]]]],"width")
	        if (is.null(w)) w=7*100
	        w
	      }
	      getheight=function() {
	        w=attr(plot.static[[n]][[input[[paste0(n,"list")]]]],"height")
	        if (is.null(w)) w=7*100
	        w
	      }
	      renderPlot({plot.static[[n]][[input[[paste0(n,"list")]]]](data)},width=getwidth,height=getheight,env=env)
	    }
	    output[[paste0(n,"plot")]]=create(n)
	  }

	  for (n in names(plot.set)) {
	    create=function(n) {
	      env=new.env()
	      env$n=n
	      getwidth=function() {
	        w=attr(plot.set[[n]],"width")
	        if (is.null(w)) w=7*100
	        w
	      }
	      getheight=function() {
	        w=attr(plot.set[[n]],"height")
	        if (is.null(w)) w=7*100
	        w
	      }
	      renderPlot({plot.set[[n]](df.set)},width=getwidth,height=getheight,env=env)
	    }
	    output[[make.names(paste0(n,"plotset"))]]=create(n)
	    e=new.env()
	    e$ddf=attr(plot.set[[n]](df.set),"df")
	    e$n=n

	    observe({
	      brushgenes=rownames(brushedPoints(ddf, input[[make.names(paste0(n,"plotsetbrush"))]]))
	      updateTextAreaInput(session, make.names(paste0(n,"plotsetgenes")), value = paste(brushgenes,collapse="\n"), label=sprintf("Selected genes (n=%d)",length(brushgenes)))
	    },env=e)

	  }

	  if (show.sessionInfo) output$sessionInfo <- renderPrint({
	    capture.output(sessionInfo())
	  })


	} # end server


	plot.static.ui=NULL
	if (length(plot.static)>0) {

	  plist=c(lapply(names(plot.static),function(n) shiny::tabPanel(n,
	                                                                shiny::selectInput(paste0(n,"list"),n,names(plot.static[[n]]),selectize=FALSE,size=10),
	                                                                shiny::plotOutput(paste0(n,"plot"))
	  )),list(title="Plots"))

	  plot.static.ui=do.call(shiny::"navbarMenu",plist)
	}

	plot.set.ui=NULL
	if (length(plot.set)>0) {

	  plist=c(lapply(names(plot.set),function(n) shiny::tabPanel(n,
	                                                       shiny::fluidRow(
	                                                         shiny::column(8,shiny::plotOutput(make.names(paste0(n,"plotset")),brush = shiny::brushOpts(id = make.names(paste0(n,"plotsetbrush"))))),
	                                                         shiny::column(4,shiny::textAreaInput(make.names(paste0(n,"plotsetgenes")), label="Selected genes",height = 300,cols=40))
	                                                      )
	  )),list(title="Global level"))

	  plot.set.ui=do.call(shiny::"navbarMenu",plist)
	}

	more=NULL
	if (show.sessionInfo)
	  more=shiny::navbarMenu("More",
	                  shiny::tabPanel("Info",shiny::verbatimTextOutput("sessionInfo"))
	  )


	ui=list(
	  shiny::tabPanel("Gene level",
	                  shiny::fluidPage(
	                    shiny::fluidRow(
        	    rclipboard::rclipboardSetup(),
        	    shiny::uiOutput("clip"),
        	    shiny::downloadButton("download1","Download"),
        	    shiny::column(12, DT::dataTableOutput('tab'))
        	  ),
        	  shiny::conditionalPanel(
        	    condition = "helpText",
        	    shiny::fluidRow(shiny::column(10, shiny::htmlOutput("helpText")))
        	  ),
        	  shiny::fluidRow(
        	    shiny::column(sizes[1], shiny::plotOutput("plot1",height = height)),
        	    shiny::conditionalPanel(
        	      condition = "plot2",
        	      shiny::column(sizes[2], shiny::plotOutput("plot2",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot3",
        	      shiny::column(sizes[3], shiny::plotOutput("plot3",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot4",
        	      shiny::column(sizes[4], shiny::plotOutput("plot4",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot5",
        	      shiny::column(sizes[5], shiny::plotOutput("plot5",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot6",
        	      shiny::column(sizes[6], shiny::plotOutput("plot6",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot7",
        	      shiny::column(sizes[7], shiny::plotOutput("plot7",height = height))
        	    ),
        	    shiny::conditionalPanel(
        	      condition = "plot8",
        	      shiny::column(sizes[8], shiny::plotOutput("plot8",height = height))
        	    )
        	  )
        	  #do.call("fluidRow",lapply(1:length(plot.funs),function(i) column(sizes[i],plotOutput(paste0("plot.funs",i),height=height))))
        	)),

        	plot.set.ui,
        	plot.static.ui,

        	more,

        	htmltools::tags$head(
        	  htmltools::tags$style(
        	    htmltools::HTML("#shiny-notification-panel {
                              top: 0;
                              bottom: unset;
                              left: 0;
                              right: 0;
                              margin-left: auto;
                              margin-right: auto;
                              width: 100%;
                              max-width: 450px;
                            }"
        	    )
        	  )
        	),

	  htmltools::tags$script(htmltools::HTML(sprintf("
        	var header = $('.navbar> .container-fluid');
          header.append('<div class=\"nav navbar-nav\" style=\"float:right\"><span class=\"navbar-brand\">%s</span></div>')",
        	                         VersionString(data)
        	)))
	)

	ui=ui[!sapply(ui,is.null)]
	myui=function(...) shiny::navbarPage(title,...)
	ui=do.call("myui",ui)

	shiny::shinyApp(ui = ui, server = server)
}



