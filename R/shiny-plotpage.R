
BindPlotTest=function(input,output,id,FUN,data,...) {
  output[[NS(id,"plot")]] = renderPlot({FUN(data,...)})
}

MetaPlotTest=function(id,data) {
  list(aest=list(
                color=c("model","manufacturer","drv")
              ),
      description="This is a test for plotPage",
      title="Test plot"
  )
}
PlotTest=function(data,aest=aes(color=trans)) {
  g=ggplot(data,utils::modifyList(aes(displ,cty),aest))+
    geom_point()+
    cowplot::theme_cowplot()
  attr(g,"table")=data.frame(x=data$displ,y=data$cty)
  attr(g,"aest")=list(
    color=rlang::as_name(as.list(aest)$colour)
  )
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
plotPage = function(id,title,plotlist=list()) {
  if (length(plotlist)==0) stop("No plot metadata given!")
  if (is.null(names(plotlist))) stop("Plots must be named!")

  # build a list of UI elements for a single plot
  buildSideList=function(plot,name) {
    pid=shiny::NS(id,name)
    re=list()
    if (!is.null(plot$description)) re = c(re,list(shiny::p(plot$description)))
    if (!is.null(plot$aest)) {
      for (aest.name in names(plot$aest)) {
        re = c(re,list(shiny::selectInput(inputId = shiny::NS(pid,"color"),label = aest.name,choices = as.list(plot$aest[[aest.name]]))))
      }
    }
    re
  }

  side=NULL
  if (length(plotlist)==1) {
    n=names(plotlist)[1]
    side=c(
      list(shiny::h3(n)),
      buildSideList(plotlist[[1]],n)
    )
    side=do.call(shiny::"sidebarPanel",side)

    main=mainPanel(
      plotOutput(outputId = shiny::NS(shiny::NS(id,n),"plot"))
    )
  } else {
    side=list(shiny::selectInput(inputId = shiny::NS(id,"selectplot"),label = "Plots:",choices = names(plotlist)))
    main=list()
    for(n in names(plotlist)) {
      cond=sprintf("input['%s']=='%s'",shiny::NS(id,"selectplot"),n)
      slist=c(
        list(condition=cond),
        buildSideList(plotlist[[n]],n)
      )
      mlist=list(condition=cond,plotOutput(outputId = shiny::NS(shiny::NS(id,n),"plot")))
      side=c(side,list(do.call(shiny::"conditionalPanel",slist)))
      main=c(main,list(do.call(shiny::"conditionalPanel",mlist)))
    }
    side=do.call(shiny::"sidebarPanel",side)
    main=do.call(shiny::"mainPanel",main)
  }

  fluidPage(
    titlePanel(title),
    sidebarLayout(sidebarPanel = side,mainPanel = main)
  )
}

# Server logic
server2 <- function(input, output) {
  BindPlotTest(input,output,shiny::NS("test","A"),PlotTest,mpg)
  BindPlotTest(input,output,shiny::NS("test","B"),PlotTest,mpg,aest=aes(color=manufacturer))

}
shinyApp(plotPage("test","Test title",list(A=MetaPlotTest("a",mpg),B=MetaPlotTest("b",mpg))), server2)

