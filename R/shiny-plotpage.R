

PlotTest=function(data,aest=aes(color=trans), meta.only=FALSE) {
  g=NA
  if (!meta.only) {
    g=ggplot(data,utils::modifyList(aes(displ,cty),aest))+
      geom_point()+
      cowplot::theme_cowplot()
    attr(g,"table")=data.frame(x=data$displ,y=data$cty)
  }
  attr(g,"aest.value")=list(
    color=if (is.null(as.list(aest)$colour)) "-" else rlang::as_name(as.list(aest)$colour),
    shape=if (is.null(as.list(aest)$shape)) "-" else rlang::as_name(as.list(aest)$shape)
  )
  attr(g,"aest")=list(
                  color=c("model","manufacturer","drv","trans","-"),
                  shape=c("model","manufacturer","drv","trans","-")
                )
  attr(g,"description")="This is a test for plotPage"
  attr(g,"title")="Test plot"
  g
}


#' Title
#'
#' @param id
#' @param title
#' @param plotlist
#'
#' @return
#' @export
#'
#' @details plot$description: character
#' plot$aest: named list; names are the aesthetic attributes that can be configured, values are vectors of choices for this (i.e. column name of the table); can be a named vector (name=display label, value=column name)
#'
#' @examples
plotPage = function(id,title,...) {
  plotlist=list(...)
  if (length(plotlist)==0) stop("No plot metadata given!")
  if (is.null(names(plotlist))) stop("Plots must be named!")

  ns=function(...) {
    v=list(...)
    re=id
    for (n in v) re=shiny::NS(re,n)
    re
  }
  .dc=function(FUN) function(list,...) {
    list=c(list(...),list)
    do.call(FUN,list)
  }
  .sidebarPanel=.dc(shiny::sidebarPanel)
  .conditionalPanel=.dc(shiny::conditionalPanel)
  .mainPanel=.dc(shiny::mainPanel)

  # build a list of UI elements for a single plot
  buildSideList=function(plot,name) {
    re=list()
    if (!is.null(plot$description)) re = c(re,list(shiny::p(plot$description)))
    if (!is.null(plot$aest)) {
      for (aest.name in names(plot$aest)) {
        sel = if (!is.null(plot$aest.value)) plot$aest.value[[aest.name]] else NULL
        re = c(re,list(shiny::selectInput(inputId = ns(name,"aest",aest.name),label = aest.name,choices = as.list(plot$aest[[aest.name]]),selected = sel)))
      }
    }
    re
  }

  side=NULL
  if (length(plotlist)==1) {
    n=names(plotlist)[1]
    side=.sidebarPanel(buildSideList(plotlist[[1]],n),shiny::h3(n))
    main=mainPanel( plotOutput(outputId = ns(n,"plot")) )
  } else {
    side=list()
    main=list()

    side=lapply(names(plotlist), function(n) {
      cond=sprintf("input['%s']=='%s'",ns("selectplot"),n)
      .conditionalPanel(buildSideList(plotlist[[n]],n),condition=cond)
    })
    main=lapply(names(plotlist), function(n) {
      cond=sprintf("input['%s']=='%s'",ns("selectplot"),n)
      .conditionalPanel(list(condition=cond,plotOutput(outputId = ns(n,"plot"))))
    })
    side=.sidebarPanel(side,
                       shiny::selectInput(inputId = ns("selectplot"),label = "Plots:",choices = names(plotlist))
                       )
    main=.mainPanel(main)
  }

  fluidPage(
    titlePanel(title),
    sidebarLayout(sidebarPanel = side,mainPanel = main)
  )
}


renderPlotPage=function(input,output,id,title=id,data,...) {
  funs=list(...)
  if (length(funs)==0) stop("No plots given!")
  if (is.null(names(funs))) stop("Plots must be named!")

  ns=function(...) {
    v=list(...)
    re=id
    for (n in v) re=shiny::NS(re,n)
    re
  }

  # obtains meta data by first trying to call the function with meta.only=TRUE; omit this if it fails
  get.meta = function(FUN) {
    re=tryCatch({FUN(data=data,meta.only=TRUE)},error=function(e) FUN(data=data))
    attributes(re)
  }
  meta=setNames(lapply(funs,function(FUN) get.meta(FUN)),names(funs))
  renderList=setNames(lapply(names(funs),function(n) renderPlot({
      para=list(data=data)

      # set aest parameter
      if (!is.null(meta[[n]]$aest)) {
        all.aest=aes()
        for (aest.name in names(meta[[n]]$aest)) {
          val=input[[ns(n,"aest",aest.name)]]
          ae=if (val!="-") ggplot2::aes_string(val) else ggplot2::aes(NULL)
          names(ae)=ggplot2::standardise_aes_names(aest.name)
          all.aest=utils::modifyList(all.aest,ae)
        }
        para=c(para,list(aest=all.aest))
      }

      do.call(funs[[n]],para)
  })),names(funs))
  for (n in names(funs)) {
    output[[ns(n,"plot")]] = renderList[[n]]
  }

  renderUI({do.call("plotPage",c(list(id=id,title = title),meta))})
}
server2 <- function(input, output) {
  output$plotpage1 = renderPlotPage(input,output,"test",data=mpg,A=PlotTest,B=Defer(PlotTest,aest=aes(color=manufacturer)))
}
shinyApp(fluidPage(
  uiOutput("plotpage1")
), server2)



