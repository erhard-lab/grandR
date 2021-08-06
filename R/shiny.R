
ServeData=function(data,aest=aes(color=Sample),time="hpi",ggparam=NULL,ggparam.ts=NULL) {
	df=GetFdrTab(data)

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
			g=ggplot(GetData(data,gene=gene,type=c("tpm"),melt=T,coldata=T),modifyList(aes_string(time,"Value"),aest))+geom_point()+scale_y_log10()+xlab(NULL)+ylab("Total TPM")
			if (!is.null(ggparam)) g=g+ggparam.ts
			g
		}
	  })

	  output$plot3=renderPlot({
		if (length(input$tab_rows_selected)==1) {
			gene=df$Symbol[input$tab_rows_selected]
			g=ggplot(GetData(data,gene=gene,type=c("new.tpm"),melt=T,coldata=T),modifyList(aes_string(time,"Value"),aest))+geom_point()+scale_y_log10()+xlab(NULL)+ylab("New TPM")
			if (!is.null(ggparam)) g=g+ggparam.ts
			g
		}
	  })

	}


	ui=fluidPage(
	  fluidRow(
	    column(12, DT::dataTableOutput('tab'))
	  ),
	  fluidRow(
	    column(4, plotOutput("plot1",height = 400)),
	    column(4, plotOutput("plot2",height = 400)),
	    column(4, plotOutput("plot3",height = 400))
	  )
	)

	shinyApp(ui = ui, server = server)
}

