
#' Serve a shiny web interface
#'
#' Fire up a shiny web server for exploratory analysis of grandR data.
#'
#' @param data the grandR object (or a file name to an rds file containing a grandR object)
#' @param table the table to display (can be NULL; see details)
#' @param sizes the widths for the gene plots to show (12 is full screen with); must be a vector as long as there are gene plots
#' @param height the height for the gene plots in pixel
#' @param plot.gene a list of gene plots; can be NULL, then the stored gene plots are used (see \link{Plots})
#' @param plot.global  a list of global plots; can be NULL, then the stored global plots are used (see \link{Plots})
#' @param df.identifier the main identifier (column name) from the table; this is used when calling the gene plot functions;
#' @param title the title to show in the header of the website
#' @param show.sessionInfo whether to show session info
#' @param help a list of characters that is shown as help text at the beginning (when no gene plot is shown); should describe the contents of your table
#'
#' @details If the table parameter is NULL, either an analysis table named "ServeGrandR" is
#' used (if it exists), otherwise the columns "Q", "LFC", "Synthesis" and "Half-life" of all analysis tables are used.
#'
#' @details The gene plots must be functions that accept two parameters: the grandR object and a gene identifier. You can either use
#' functions directly (e.g. \code{plot.gene=list(PlotGeneOldVsNew)}), or use \link{Defer} in cases you need to specify additional parameters,
#' e.g. \code{plot.gene=list(Defer(PlotGeneOldVsNew,log=FALSE))}. The global plots are functions accepting a single parameter (the grandR object). Here
#' the use of \link{Defer} is encouraged due to its caching mechanism.
#'
#' @return a shiny web server
#' @export
#' @concept shiny
ServeGrandR=function(data,
                   table=NULL,
                   sizes=NA,height=400,
                   plot.gene=NULL,
                   plot.global=NULL,
                   df.identifier="Symbol",
                   title=Title(data),
                   show.sessionInfo=FALSE,
                   help=list(".Q: multiple testing corrected p values",".LFC: log2 fold changes") ) {

  checkPackages(c("shiny","rclipboard","DT","htmltools"))

  if (is.character(data) && file.exists(data)) data = readRDS(data)

  plot.static=list()
  df=table
  if (is.null(df)) {
    df=if ("ServeGrandR" %in% Analyses(data)) GetAnalysisTable(data,analyses = "ServeGrandR",regex = FALSE,gene.info = FALSE,prefix.by.analysis=FALSE) else GetAnalysisTable(data,columns="Synthesis|Half-life|LFC|Q|log2FC|ROPE",gene.info = FALSE)
  }

  if (ncol(df)==0) stop("Empty table given!")

  if (is.null(plot.gene)) plot.gene=data$plots$gene
  if (is.null(plot.global)) plot.global=data$plots$global
  if (is.null(plot.gene)) plot.gene=PlotGeneGroupsBars

  if (!is.list(plot.gene)) plot.gene=list(plot.gene)
  if (length(plot.gene)==0) plot.gene=list(PlotGeneGroupsBars)
  if (length(sizes)==1 && is.na(sizes)) sizes=rep(floor(12/min(4,length(plot.gene))),length(plot.gene))
  if (length(sizes)!=length(plot.gene)) stop("sizes need to be length 1 or same length as plots!")
  sizes=c(sizes,rep(1,8))

  if (!df.identifier %in% names(df) && !is.null(rownames(df))) df=cbind(setNames(data.frame(rownames(df),stringsAsFactors = FALSE),df.identifier),df,stringsAsFactors = FALSE)

#    plot.static=lapply(plot.static, function(p) if (is.function(p)) p(data) else p)

  if (!is.null(help) && is.list(help)) help=sprintf("<span style='padding-top:25px;'><span class='help-block well'>Table columns:%s</span></span>", paste(sapply(help,function(s) sprintf("<li><span>%s</span></li>",s)),collapse="\n"))
	server=function(input, output,session) {
	  # R CMD check guard for non-standard evaluation
	  ddf <- NULL

	  df.rounded = as.data.frame(lapply(df,function(v) if(is.numeric(v)) round(v,5) else v),check.names=FALSE,stringsAsFactors = FALSE)

	  dttab=DT::datatable(df.rounded,
	                      callback = DT::JS("$('div#buttons').css('float','left').css('margin-right','50px'); $('div#clip').css('float','left'); $('div#buttons').append($('#downloadraw')); $('div#buttons').append($('#download1')); $('div#buttons').append($('#clip')); "),
	                      selection = 'single',
	                      rownames = FALSE,
	                      filter = "top",
	                      escape = FALSE,
	                      options = list(
	                        pageLength = 10,
	                        lengthMenu =list(c(5, 10, 25, 50, 100,-1),c(5, 10, 25, 50, 100,"All")),
	                        dom = '<"#buttons">lfrtip'
	                      ))
	  dttab=DT::formatRound(dttab,names(df)[sapply(df,class)=="numeric"], 2)
	  if (any(grepl("LFC$",names(df)))) dttab=DT::formatRound(dttab,names(df)[grepl("\\.LFC$",names(df))], 2)
	  if (any(grepl("\\.Q$",names(df)))) dttab=DT::formatSignif(dttab,names(df)[grepl("\\.Q$",names(df))], 2)
	  if (any(grepl("\\.P$",names(df)))) dttab=DT::formatSignif(dttab,names(df)[grepl("\\.P$",names(df))], 2)
	  if (any(grepl("\\.Half-life$",names(df)))) dttab=DT::formatRound(dttab,names(df)[grepl("\\.Half-life$",names(df))], 2)
	  if (any(grepl("\\.Synthesis$",names(df)))) dttab=DT::formatRound(dttab,names(df)[grepl("\\.Synthesis$",names(df))], 2)
	  output$tab <- DT::renderDataTable(dttab)
	  output$download1 <- shiny::downloadHandler(
	    filename = function() {
	      paste0(title,"-", Sys.Date(), ".tsv")
	    },
	    content = function(file) {
	      write.table(df[input$tab_rows_all,], file,row.names=F,col.names=T,quote=F,sep="\t")
	    }
	  )

	  mods=list(`Raw data`=unlist(lapply(Slots(data),function(sl) if(sl %in% c("ntr","alpha","beta")) sl else paste0(c("","new.","old."),sl))),Analyses=Analyses(data))
	  shiny::observeEvent(input$downloadraw, {
	    shiny::showModal(shiny::modalDialog(
	      shiny::selectInput("datamodality","Data modality",choices=mods,selected=DefaultSlot(data),selectize = FALSE),

	      title="Download data",easyClose=TRUE,footer=htmltools::tagList(
	        shiny::modalButton("Cancel"),
	        shiny::downloadButton("downloadrawdoit", "OK")
	      )
	    ))
	  })
	  output$downloadrawdoit <- shiny::downloadHandler(
	    filename = function() {
	      paste0(title,"-", Sys.Date(), ".tsv.gz")
	    },
	    content = function(file) {
	      on.exit(shiny::removeModal())
	      ggg=as.character(df[input$tab_rows_all,1])
	      tab=GetTable(data,type=input$datamodality,ntr.na = FALSE,gene.info = TRUE,genes = ggg)
	      write.table(tab, gzfile(file),row.names=F,col.names=T,quote=F,sep="\t")
	    }
	  )

	  output$clip <- shiny::renderUI({
	    nn=if(.row_names_info(df)<0) df[input$tab_rows_all,1] else rownames(df)[input$tab_rows_all]
	    rclipboard::rclipButton("clipbtn", "Copy", paste(nn,collapse="\n"), modal=TRUE,icon=shiny::icon("clipboard"))
	  })
	  shiny::observeEvent(input$clipbtn, {shiny::showNotification(
	    sprintf("Copied %d names",length(input$tab_rows_all)),
	    duration = 2,
	    type = "message"
	  )})


	  output$plot1=shiny::renderPlot({ if (length(input$tab_rows_selected)==1) plot.gene[[1]](data=data,gene=df[[df.identifier]][input$tab_rows_selected]) })
	  output$plot2=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=2) plot.gene[[2]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot3=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=3) plot.gene[[3]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot4=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=4) plot.gene[[4]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot5=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=5) plot.gene[[5]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot6=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=6) plot.gene[[6]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot7=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=7) plot.gene[[7]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$plot8=shiny::renderPlot({ if (length(input$tab_rows_selected)==1 && length(plot.gene)>=8) plot.gene[[8]](data=data,gene=df[[df.identifier]][input$tab_rows_selected])  })
	  output$helpText=shiny::renderText({ if (length(input$tab_rows_selected)==0 && !is.null(help)) help  })

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
	      shiny::renderPlot({plot.static[[n]][[input[[paste0(n,"list")]]]](data)},width=getwidth,height=getheight,env=env)
	    }
	    output[[paste0(n,"plot")]]=create(n)
	  }

	  highlighted.genes <- shiny::reactiveValues(genes = c())
	  for (n in names(plot.global)) {
	    create=function(n) {
	      env=new.env()
	      env$n=n
	      getwidth=function() {
	        w=attr(plot.global[[n]],"width")
	        if (is.null(w)) w=7*100
	        w
	      }
	      getheight=function() {
	        w=attr(plot.global[[n]],"height")
	        if (is.null(w)) w=7*100
	        w
	      }
	      shiny::renderPlot({plot.global[[n]](data,highlight=highlighted.genes$genes,label=df[[df.identifier]][input$tab_rows_selected])},width=getwidth,height=getheight,env=env)
	    }
	    output[[make.names(paste0(n,"plotset"))]]=create(n)
	    e=new.env()
	    e$ddf=attr(plot.global[[n]](data),"df")
	    e$n=n

	    shiny::observe({
	      brushgenes=rownames(shiny::brushedPoints(ddf, input[[make.names(paste0(n,"plotsetbrush"))]]))
	      shiny::updateTextAreaInput(session, make.names(paste0(n,"plotsetgenes")), value = paste(brushgenes,collapse="\n"), label=sprintf("Selected genes (n=%d)",length(brushgenes)))
	    },env=e)
	  }

	  lapply(names(plot.global),function(n) {
	    shiny::observeEvent(input$tab_rows_selected, {
	      session$resetBrush(make.names(paste0(n,"plotsetbrush")))
	    })

	    shiny::observeEvent(input[[make.names(paste0(n,"sethighlight"))]], {
	      highlighted.genes$genes <- strsplit(input[[make.names(paste0(n,"plotsetgenes"))]],"\n")[[1]]
	      session$resetBrush(make.names(paste0(n,"plotsetbrush")))
	    })

	    shiny::observeEvent(input[[make.names(paste0(n,"addhighlight"))]], {
	      highlighted.genes$genes <- c(highlighted.genes$genes,strsplit(input[[make.names(paste0(n,"plotsetgenes"))]],"\n")[[1]])
	      session$resetBrush(make.names(paste0(n,"plotsetbrush")))
	    })
	  })

	  if (show.sessionInfo) output$sessionInfo <- shiny::renderPrint({
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

	plot.global.ui=NULL
	if (length(plot.global)>0) {

	  plist=c(lapply(names(plot.global),function(n) shiny::tabPanel(n,
	                                                       shiny::fluidRow(
	                                                         shiny::column(8,shiny::plotOutput(make.names(paste0(n,"plotset")),brush = shiny::brushOpts(id = make.names(paste0(n,"plotsetbrush"))))),
	                                                         shiny::column(4,
	                                                                       shiny::textAreaInput(make.names(paste0(n,"plotsetgenes")), label="Selected genes",height = 300,cols=40),
	                                                                       shiny::actionButton(make.names(paste0(n,"sethighlight")), label="Set highlight"),
	                                                                       shiny::actionButton(make.names(paste0(n,"addhighlight")), label="Add highlight")
	                                                                       )
	                                                      )
	  )),list(title="Global level"))

	  plot.global.ui=do.call(shiny::"navbarMenu",plist)
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
        	    shiny::downloadButton("download1","Table"),
        	    shiny::actionButton("downloadraw","Data",icon = shiny::icon("download")),
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

        	plot.global.ui,
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
          header.append('<div class=\"nav navbar-nav\" style=\"float:right\"><span class=\"navbar-brand\">grandR v%s</span></div>')",
	                                                 packageVersion("grandR")
        	)))
	)

	ui=ui[!sapply(ui,is.null)]
	myui=function(...) shiny::navbarPage(title,...)
	ui=do.call("myui",ui)

	shiny::shinyApp(ui = ui, server = server)
}



