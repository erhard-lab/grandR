
ServeData=function(data) {
	df=GetFdrTab(data)

	server=function(input, output,session) {
	  options(DT.options = list(pageLength = 12))
	  output$tab <- DT::renderDataTable(DT::datatable(df, selection = 'single',rownames = FALSE, escape=-1)) # %>%formatRound(names(df)[grepl("q$",names(df))], 2)

	  output$plot1=renderPlot({
		if (length(input$tab_rows_selected)==1) {
			gene=df$Symbol[input$tab_rows_selected]
			ggplot(GetData(data,gene=gene,type=c("tpm","ntr"),melt=F,coldata=T),aes(tpm,ntr,color=hpi,shape=Condition))+geom_point(size=3)+scale_x_log10(breaks = scales::pretty_breaks(n = 10))
		}
	  })

	  output$plot2=renderPlot({
		if (length(input$tab_rows_selected)==1) {
			gene=df$Symbol[input$tab_rows_selected]
			ggplot(GetData(data,gene=gene,type=c("tpm"),melt=T,coldata=T),aes(hpi,Value,color=Condition))+geom_point()+scale_y_log10(breaks = scales::pretty_breaks(n = 10))+xlab(NULL)+ylab("Total TPM")
		}
	  })

	  output$plot3=renderPlot({
		if (length(input$tab_rows_selected)==1) {
			gene=df$Symbol[input$tab_rows_selected]
			ggplot(GetData(data,gene=gene,type=c("new.tpm"),melt=T,coldata=T),aes(hpi,Value,color=Condition))+geom_point()+scale_y_log10(breaks = scales::pretty_breaks(n = 10))+xlab(NULL)+ylab("New TPM")
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

