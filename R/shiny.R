
#' Serve a shiny web interface
#'
#' Fire up a shiny web server for exploratory analysis of grandR data.
#'
#' @param data the grandR object (or a file name to an rds file containing a grandR object)
#' @param table the table to display (can be NULL or a named list; see details)
#' @param sizes the widths for the gene plots to show (12 is full screen with); must be a vector as long as there are gene plots
#' @param height the height for the gene plots in pixel
#' @param plot.gene a list of gene plots; can be NULL, then the stored gene plots are used (see \link{Plots})
#' @param plot.global  a list of global plots; can be NULL, then the stored global plots are used (see \link{Plots})
#' @param plot.window a list of static plots to show in a floating window; see details
#' @param highlight a vector of gene names that are highlighted in the beginning
#' @param df.identifier the main identifier (column name) from the table; this is used when calling the gene plot functions;
#' @param title the title to show in the header of the website
#' @param show.sessionInfo whether to show session info
#' @param help a list of characters that is shown as help text at the beginning (when no gene plot is shown); should describe the contents of your table
#'
#' @details If the table parameter is NULL, either an analysis table named "ServeGrandR" is
#' used (if it exists), otherwise the columns "Q", "LFC", "Synthesis" and "Half-life" of all analysis tables are used. If it is a list, a menu is created in the navbar
#'
#' @details plot.window must be a list of plotting functions that receive the grandR object and return a plot. It can also be a list of list, in which case more than one plotting windows are generated.
#' Each plot will be rendered with a size of 350x350.
#'
#' @details The gene plots must be functions that accept two parameters: the grandR object and a gene identifier. You can either use
#' functions directly (e.g. \code{plot.gene=list(PlotGeneOldVsNew)}), or use \link{Defer} in cases you need to specify additional parameters,
#' e.g. \code{plot.gene=list(Defer(PlotGeneOldVsNew,log=FALSE))}. The global plots are functions accepting a single parameter (the grandR object). Here
#' the use of \link{Defer} is encouraged due to its caching mechanism.
#'
#' @return a shiny web server
#' @export
#'
#' @examples
#' \dontrun{
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' sars <- Pairwise(sars,contrasts = GetContrasts(sars,contrast = c("Condition","SARS","Mock")))
#' sars <- AddGenePlot(sars,"timecourse",
#'                     Defer(PlotGeneProgressiveTimecourse,steady.state=c(Mock=TRUE,SARS=FALSE)))
#' sars <- AddGlobalPlot(sars,"Vulcano",VulcanoPlot)
#' ServeGrandR(sars)
#'
#' }
#'
#' @concept shiny
ServeGrandR=function(data,
                     table=NULL,
                     sizes=NA,height=400,
                     plot.gene=NULL,
                     plot.global=NULL,
                     plot.window=NULL,
                     highlight=NULL,
                     df.identifier="Symbol",
                     title=Title(data),
                     show.sessionInfo=FALSE,
                     help=list(".Q: multiple testing corrected p values",".LFC: log2 fold changes") ) {

  checkPackages(c("shiny","rclipboard","DT","htmltools"))

  if (is.character(data) && file.exists(data)) data = readRDS(data)

  plot.static=list()
  if (is.null(table)) {
    table=if ("ServeGrandR" %in% Analyses(data)) GetAnalysisTable(data,analyses = "ServeGrandR",regex = FALSE,gene.info = FALSE,prefix.by.analysis=FALSE) else GetAnalysisTable(data,columns="Synthesis|Half-life|LFC|Q|log2FC|ROPE",gene.info = FALSE)
  }


  if (is.data.frame(table)) table = list(Table=table)

  for (i in 1:length(table)) {
    df=table[[i]]
    if (ncol(df)==0) stop("Empty table given!")
    if (!df.identifier %in% names(df) && !is.null(rownames(df))) df=cbind(setNames(data.frame(rownames(df),stringsAsFactors = FALSE),df.identifier),df,stringsAsFactors = FALSE)
    if (!df.identifier %in% names(df)) stop("Neither identifier nor row names found in table")
    table[[i]] = df
  }
  #gene.level.map=setNames(make.names(names(table)),names(table))
  #names(table) = make.names(names(table))

  if (is.null(plot.gene)) plot.gene=data$plots$gene
  if (is.null(plot.global)) plot.global=data$plots$global
  if (is.null(plot.gene)) plot.gene=PlotGeneGroupsBars

  if (!is.list(plot.gene)) plot.gene=list(plot.gene)
  if (length(plot.gene)==0) plot.gene=list(PlotGeneGroupsBars)
  if (length(sizes)==1 && is.na(sizes)) sizes=rep(floor(12/min(4,length(plot.gene))),length(plot.gene))
  if (length(sizes)!=length(plot.gene)) stop("sizes need to be length 1 or same length as plots!")
  sizes=c(sizes,rep(1,8))

  if (is.null(plot.window)) {
    if (!is.null(data$plots$floating)) {
      n=names(data$plots$floating)
      nodot = !grepl(".",n,fixed=TRUE)
      n[nodot]=paste0("Plots.",n[nodot])

      cats = unique(sapply(strsplit(n,".",fixed=TRUE),function(cc) cc[1]))
      plot.window = lapply(cats,function(cc) {
        use = substr(n,1,nchar(cc))==cc
        setNames(data$plots$floating[use],substr(n[use],nchar(cc)+2,nchar(n[use])))
      })
      names(plot.window) = cats

    } else plot.window = list()
  }
  if (length(plot.window)>0 && !is.list(plot.window[[1]])) plot.window=list(Plots=plot.window)
  # after this, plot.window is either an empty list, or a list of lists.
  if (is.null(names(plot.window)) && length(plot.window)>0) stop("plot.window must have names!")


  plot.wins = lapply(names(plot.window),function(n)CreateWindows(paste0("floating-",n),title=n,plots=plot.window[[n]]))


  #    plot.static=lapply(plot.static, function(p) if (is.function(p)) p(data) else p)

  htmls = list.files(pattern="html$")
  names(htmls) = gsub(".html","",htmls)
  shiny::addResourcePath("htmls", getwd())

  if (!is.null(help) && is.list(help)) help=sprintf("<span style='padding-top:25px;'><span class='help-block well'>Table columns:%s</span></span>", paste(sapply(help,function(s) sprintf("<li><span>%s</span></li>",s)),collapse="\n"))
  server=function(input, output,session) {


    highlighted.genes <- shiny::reactiveValues(genes = highlight,selected.gene=NULL,filtered.rows=NULL,active.table=NULL)

    win_ids = lapply(plot.wins,function(pw) pw$server(data))

    for (n in names(table)) {
      create.tab = function(ns,df) {
        # R CMD check guard for non-standard evaluation
        output=list()
        ddf <- NULL

        df.rounded = as.data.frame(lapply(df,function(v) if(is.numeric(v)) round(v,5) else v),check.names=FALSE,stringsAsFactors = FALSE)
        plot.window.adding = ""
        for (i in 1:length(plot.window)) {
          plot.window.adding = paste(plot.window.adding,sprintf("$('div[id=\"%s\"] > div > div#buttons').append($('[id=\"%s\"]'));",ns("tab"),ns(paste0("btn-plot-window-",names(plot.window)[i]))))
          obs=function(i) {
            btnid = ns(paste0("btn-plot-window-",names(plot.window)[i]))
            shiny::observeEvent(input[[btnid]], {
              shinyjs::js$toggleZ(win_ids[[i]])
              shinyjs::show(win_ids[[i]])
            })
          }
          obs(i)
        }
        dttab=DT::datatable(df.rounded,
                            callback = DT::JS(sprintf("$('div[id=\"%s\"] > div > div#buttons').css('float','left').css('margin-right','50px'); $('div[id=\"%s\"]').css('float','left'); $('div[id=\"%s\"] > div > div#buttons').append($('[id=\"%s\"]')); $('div[id=\"%s\"] > div > div#buttons').append($('[id=\"%s\"]')); $('div[id=\"%s\"] > div > div#buttons').append($('[id=\"%s\"]')); $('div[id=\"%s\"] > div > div#buttons').append($('[id=\"%s\"]')); %s",ns("tab"),ns("clip"),ns("tab"),ns("pdf"),ns("tab"),ns("downloadraw"),ns("tab"),ns("download1"),ns("tab"),ns("clip"),plot.window.adding)),
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
            utils::write.table(df[highlighted.genes$filtered.rows,], file,row.names=F,col.names=T,quote=F,sep="\t")
          }
        )

        shiny::observeEvent(input[[ns("pdf")]], {
          shiny::showModal(shiny::modalDialog(
            shiny::numericInput(ns("pdfwidth"),"Width",value=7,min=2,max=20,step=1),
            shiny::numericInput(ns("pdfheight"),"Height",value=4,min=2,max=20,step=1),

            title="Generate pdf",easyClose=TRUE,footer=htmltools::tagList(
              shiny::modalButton("Cancel"),
              shiny::downloadButton(ns("pdfdoit"), "OK")
            )
          ))
        })
        output$pdfdoit <- shiny::downloadHandler(
          filename = function() {
            if (!is.null(highlighted.genes$selected.gene)) paste0(highlighted.genes$selected.gene,"-", Sys.Date(), ".pdf") else paste0(Sys.Date(), ".pdf")
          },
          content = function(file) {
            on.exit(shiny::removeModal())
            grDevices::pdf(file,width=input[[ns("pdfwidth")]],height=input[[ns("pdfheight")]])
            for (i in 1:length(plot.gene)) print(plot.gene[[i]](data=data,gene=highlighted.genes$selected.gene))
            grDevices::dev.off()
          }
        )

        mods=list(`Raw data`=unlist(lapply(Slots(data),function(sl) if(sl %in% c("ntr","alpha","beta")) sl else paste0(c("","new.","old."),sl))),Analyses=Analyses(data))
        shiny::observeEvent(input[[ns("downloadraw")]], {
          shiny::showModal(shiny::modalDialog(
            shiny::selectInput(ns("datamodality"),"Data modality",choices=mods,selected=DefaultSlot(data),selectize = FALSE),

            title="Download data",easyClose=TRUE,footer=htmltools::tagList(
              shiny::modalButton("Cancel"),
              shiny::downloadButton(ns("downloadrawdoit"), "OK")
            )
          ))
        })
        output$downloadrawdoit <- shiny::downloadHandler(
          filename = function() {
            paste0(title,"-", Sys.Date(), ".tsv.gz")
          },
          content = function(file) {
            on.exit(shiny::removeModal())
            ggg=as.character(df[highlighted.genes$filtered.rows,1])
            tab=GetTable(data,type=input[[ns("datamodality")]],ntr.na = FALSE,gene.info = TRUE,genes = ggg)
            utils::write.table(tab, gzfile(file),row.names=F,col.names=T,quote=F,sep="\t")
          }
        )

        output$clip <- shiny::renderUI({
          nn=if(.row_names_info(df)<0) df[highlighted.genes$filtered.rows,1] else rownames(df)[highlighted.genes$filtered.rows]
          rclipboard::rclipButton(ns("clipbtn"), "Copy", paste(nn,collapse="\n"), modal=TRUE,icon=shiny::icon("clipboard"))
        })
        shiny::observeEvent(input[[ns("clipbtn")]], {shiny::showNotification(
          sprintf("Copied %d names",length(highlighted.genes$filtered.rows)),
          duration = 2,
          type = "message"
        )})


        shiny::observeEvent(input[[ns("tab_rows_selected")]], ignoreNULL = FALSE, {
          highlighted.genes$selected.gene = if (is.null(input[[ns("tab_rows_selected")]])) NULL else df[[df.identifier]][input[[ns("tab_rows_selected")]]]
        })
        shiny::observeEvent( input[[ns("tab_rows_all")]], ignoreNULL = FALSE, {
          highlighted.genes$filtered.rows = input[[ns("tab_rows_all")]]
        })

        output
      }

      ns = shiny::NS(paste0("table",n))
      elements = create.tab(ns,table[[n]])
      for (n in names(elements)) {
        output[[ns(n)]] = elements[[n]]
      }
    }

    shiny::observeEvent(input$mainnavbar, {
      if (input$mainnavbar %in% names(table)) {
        ns = shiny::NS(paste0("table",input$mainnavbar))
        highlighted.genes$selected.gene <- if (is.null(input[[ns("tab_rows_selected")]])) NULL else df[[df.identifier]][input[[ns("tab_rows_selected")]]]
        highlighted.genes$filtered.rows = input[[ns("tab_rows_all")]]
        highlighted.genes$active.table <- input$mainnavbar
      }
    })


    for (i in 1:(length(plot.gene))) {
      create.plotgene = function(i) {
        env=new.env()
        env$i=i
        shiny::renderPlot({ if (!is.null(highlighted.genes$selected.gene)) plot.gene[[i]](data=data,gene=highlighted.genes$selected.gene) },env=env)
      }
      output[[paste0("plot",i)]]=create.plotgene(i)
    }

    output$helpText=shiny::renderText({ if (is.null(highlighted.genes$selected.gene) && !is.null(help)) help  })

    for (n in names(plot.static)) {
      create.plotstatic=function(n) {
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
      output[[paste0(n,"plot")]]=create.plotstatic(n)
    }


    for (n in names(plot.global)) {
      create.plotglobal=function(n) {
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
        shiny::renderPlot({
          re=plot.global[[n]](data,highlight=highlighted.genes$genes,label=highlighted.genes$selected.gene)
          ddf$ddf=attr(re,"df")
          re
          },width=getwidth,height=getheight,env=env)
      }
      output[[make.names(paste0(n,"plotset"))]]=create.plotglobal(n)
      e=new.env()
      e$n=n
      ddf <- shiny::reactiveValues(ddf=NULL)
      shiny::observe({
        brushgenes=if (is.null(ddf$ddf)) NULL else rownames(shiny::brushedPoints(ddf$ddf, input[[make.names(paste0(n,"plotsetbrush"))]]))
        shiny::updateTextAreaInput(session, make.names(paste0(n,"plotsetgenes")), value = paste(brushgenes,collapse="\n"), label=sprintf("Selected genes (n=%d)",length(brushgenes)))
      },env=e)
    }
    shiny::observe({
      shiny::updateTextAreaInput(session, "highlightedgenes", value = paste0(highlighted.genes$genes,collapse="\n"), label=sprintf("Highlighted genes (n=%d)",length(highlighted.genes$genes)))
    })
    shiny::observeEvent(input[["updatehighlight"]], {
      highlighted.genes$genes <- strsplit(input[["highlightedgenes"]],"\n")[[1]]
      for (n in names(plot.global)) session$resetBrush(make.names(paste0(n,"plotsetbrush")))
    })


    lapply(names(plot.global),function(n) {
      shiny::observeEvent(highlighted.genes$selected.gene, {
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
      utils::capture.output(utils::sessionInfo())
    })

  } # end server


  html.ui=NULL
  html.list.ui = ""
  if (length(htmls)>0) {
    html.ui=shiny::navbarMenu(title = "Reports")
    appends = sapply(names(htmls),function(name) {
      sprintf("$('.dropdown-toggle[data-value=\"Reports\"] + .dropdown-menu').append('<li><a target=\"_blank\" href=\"%s\">%s</a></li>');",paste0("htmls/",htmls[name]),name)
    })
    html.list.ui = paste(appends,collapse="\n")
  }


  plot.static.ui=NULL
  if (length(plot.static)>0) {

    plist=c(lapply(names(plot.static),function(n) shiny::tabPanel(n,
                                                                  shiny::selectInput(paste0(n,"list"),n,names(plot.static[[n]]),selectize=FALSE,size=10),
                                                                  shiny::plotOutput(paste0(n,"plot"))
    )),list(title="Plots"))

    plot.static.ui=do.call(shiny::"navbarMenu",plist)
  }

  plot.gene.ui=NULL
  if (length(table)==1) {
    ns=shiny::NS(paste0("table",names(table)[1]))
    plot.window.buttons = lapply(names(plot.window),function(n)shiny::actionButton(ns(paste0("btn-plot-window-",n)),n,icon=shiny::icon("chart-column")))
    plot.gene.ui=shiny::tabPanel("Gene level",
                    shiny::fluidPage(
                      shiny::fluidRow(
                        rclipboard::rclipboardSetup(),
                        shiny::uiOutput(ns("clip")),
                        shiny::downloadButton(ns("download1"),"Table"),
                        shiny::actionButton(ns("downloadraw"),"Data",icon = shiny::icon("download")),
                        shiny::actionButton(ns("pdf"),"PDF",icon = shiny::icon("file-pdf")),
                        plot.window.buttons,
                        shiny::column(12, DT::dataTableOutput(ns('tab')))
                      )
                    ))
   } else {
     plist=c(lapply(names(table),function(n) {
       ns=shiny::NS(paste0("table",n))
       plot.window.buttons = lapply(names(plot.window),function(n)shiny::actionButton(ns(paste0("btn-plot-window-",n)),n,icon=shiny::icon("chart-column")))
       shiny::tabPanel(n,
                       shiny::fluidPage(
                         shiny::fluidRow(
                           rclipboard::rclipboardSetup(),
                           shiny::uiOutput(ns("clip")),
                           shiny::downloadButton(ns("download1"),"Table"),
                           shiny::actionButton(ns("downloadraw"),"Data",icon = shiny::icon("download")),
                           shiny::actionButton(ns("pdf"),"PDF",icon = shiny::icon("file-pdf")),
                           plot.window.buttons,
                           shiny::column(12, DT::dataTableOutput(ns('tab')))
                         )
                       ))
     }
     ),list(title="Gene level"))
     plot.gene.ui=do.call(shiny::"navbarMenu",plist)
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
    plist=c(list(shiny::tabPanel("Highlighted Genes",
                                 shiny::fluidRow(
                                   shiny::column(4,
                                                 shiny::textAreaInput("highlightedgenes", label="Highlighted genes",height = 300,cols=40),
                                                 shiny::actionButton("updatehighlight", label="Update highlight")
                                   )
                                 )
    )),plist)

    plot.global.ui=do.call(shiny::"navbarMenu",plist)
  }

  more=NULL
  if (show.sessionInfo)
    more=shiny::navbarMenu("More",
                           shiny::tabPanel("Info",shiny::verbatimTextOutput("sessionInfo"))
    )


  window.ui = NULL
  if (length(plot.window)>0) {
    window.ui = c(list(InitWindows()),lapply(plot.wins,function(pw) pw$ui()))
  }

  ui=list(
    plot.gene.ui,
    plot.global.ui,
    plot.static.ui,
    html.ui,
    more
  )


  jslist = list(
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
          %s
        	var header = $('.navbar> .container-fluid');
          header.append('<div class=\"nav navbar-nav\" style=\"float:right\"><span class=\"navbar-brand\">grandR v%s</span></div>')",
                                                 html.list.ui,
                                                 utils::packageVersion("grandR")
  )))
  )

  footer = list(
    shiny::conditionalPanel(
      condition=sprintf("[%s,'Gene level'].includes(input.mainnavbar)",paste0("'",names(table),"'",collapse=",")),
      shiny::conditionalPanel(
        condition = "helpText",
        shiny::fluidRow(shiny::column(10, shiny::htmlOutput("helpText")))
      ),
      do.call(shiny::"fluidRow",lapply(1:length(plot.gene),function(i) shiny::column(sizes[i], shiny::plotOutput(paste0("plot",i),height = height))))
    ),
    window.ui
  )

  ui=ui[!sapply(ui,is.null)]
  myui=function(...) shiny::navbarPage(title,id = "mainnavbar",...,footer=footer,header=jslist)
  ui=do.call("myui",ui)

  shiny::shinyApp(ui = ui, server = server)
}



InitWindows = function() {
  jscode = "shinyjs.toggleZ = function(params){ $('.plot-window').removeClass('plot-window-topmost'); $('#'+params).addClass('plot-window-topmost');}"
  list(shinyjs::useShinyjs(),
       shinyjs::extendShinyjs(text=jscode,functions="toggleZ"),
       htmltools::tags$head(htmltools::tags$style(htmltools::HTML("
      .plot-window {
        position: fixed;
        background-color: white;
        border: 1px solid #ddd;
        box-shadow: 0px 0px 10px #aaa;
        padding: 2px;
        display: none;
        height: auto;
      }
      .plot-window-shown {
        display: block;
        z-index: 100;
      }
      .plot-window-topmost {
        z-index: 101;
      }
      .window-header {
        padding: 0px;
        background-color: #f7f7f7;
        border-bottom: 1px solid #ddd;
        cursor: move;
      }
      .window-title {
        font-weight: bold;
      }
      .btn-close {
        float: right;
        cursor: pointer;
        padding: 0px 6px;
        font-size: 13px;
      }
      .window-body {
        padding: 1px;
      }
    "))))
}

CreateWindows = function(id,plots,title="Plots",width=350,height=350,nrow=NULL,ncol=NULL,x=50,y=-20,selection=NA) {
  if (is.na(selection)) selection=length(plots)>4
  xco = if (x<0) "bottom" else "top"
  yco = if (y<0) "right" else "left"
  x=abs(x)
  y=abs(y)
  if (is.null(nrow) && is.null(ncol)) {
    ncol=if(selection==TRUE) 1 else ceiling(sqrt(length(plots)))
  }
  if (is.null(ncol)) ncol = if(selection==TRUE) 1 else ceiling(length(plots)/nrow)

  list(ui=function() {
    ns <- shiny::NS(id)

    selectlist=NULL
    if (selection==TRUE) {
      if (is.null(names(plots))) stop("Use a named list of plots for selection mode!")
      selectlist = htmltools::tags$select(id=ns("plotSelect"),style="width: 120px; margin-left: 10px;")
      selectlist= htmltools::tagAppendChildren(selectlist,lapply(1:length(plots),function(i) {
        re = htmltools::tagAppendChildren(htmltools::tags$option(value=names(plots)[i],names(plots)[i]))
        if (i==1) re=htmltools::tagAppendAttributes(re , selected=NA  )
        re
      }))
      plotctrl = shiny::plotOutput(ns("plot"),height = sprintf("%.0fpx",height))
    } else {
      plotctrl = lapply(1:length(plots),function(i) shiny::plotOutput(ns(paste0("plot",i)),height = sprintf("%.0fpx",height)))
    }
    htmltools::tagList(
      shinyjqui::jqui_draggable(
        htmltools::div(id = ns("plotWindow"), class = "plot-window", style = sprintf("width: %.0fpx; %s: %.0fpx; %s: %.0fpx;",width*ncol,xco,x,yco,y),
                       htmltools::div(class = "window-header",
                                      htmltools::span(title, class = "window-title"),
                                      selectlist,
                                      shiny::actionButton(ns("closePlotWindow"), label = "x", class = "btn-close")
                       ),
                       htmltools::div(class = "window-body",
                                      htmltools::div(style = sprintf("display: grid; grid-template-columns: repeat(%.0f, 1fr); gap: 1px;",ncol),plotctrl)
                       )
        )
      )
    )
  },
  server=function(data) {
    shiny::moduleServer(id, function(input, output, session) {
      ns <- session$ns
      for (i in 1:length(plots)) {
        oo=function(i) output[[paste0("plot",i)]] = shiny::renderPlot({ plots[[i]](data) })
        oo(i)
      }
      shiny::observeEvent(input$closePlotWindow, {
        shinyjs::js$toggleZ(shiny::NS(id)("plotWindow"))
        shinyjs::hide("plotWindow")
      })
      shiny::observeEvent(input$plotSelect, {
        output[["plot"]] = shiny::renderPlot({ plots[[input$plotSelect]]() })
      })
    })
    shiny::NS(id)("plotWindow")
  })
}



