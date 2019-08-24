# setwd("C:/Users/User/Desktop") # comment this out if pwd() returns the Desktop directory,
rm(list=ls()) # remove data from previous execution
#### Insert data
evaluators = c("utility","delay","exectime")
for (evaluator in evaluators){
  filename1 = paste("C:/Users/User/Desktop/thesis_solver_results/solver_",evaluator,"_results.txt",sep='') # read 
  filename2 = paste("C:/Users/User/Desktop/thesis_figures/evaluation_results/",evaluator,"_fig.png",sep='') # write
  filename3 = paste("C:/Users/User/Desktop/thesis_figures/evaluation_results/",evaluator,"_results.tex",sep='') # write
  
  dataa = read.csv(filename1, header = TRUE)
  head(dataa)
  dev.off(dev.list()["RStudioGD"])
  
  # require(graphics) # preliminary for postcript command
  # postscript(file=file,horiz=TRUE,onefile=FALSE,width=8.5,height=11,paper=letter) 
  titlees = c("Satellites: 4 - Stations: 2","Satellites: 16 - Stations: 4",
              "Satellites: 32 - Stations: 5","Satellites: 64 - Stations: 6")
  
  
  par(mfrow=c(2,2))
  for (i in 1:4){
    data_ngroup_1 = dataa[which(dataa$Node.Group==i),] # node group 1
    epochs = data_ngroup_1$Total_epochs # x axis
    heuristic_1_closest =  data_ngroup_1$heuristic_1_closest # y axis
    heuristic_1_fastest =  data_ngroup_1$heuristic_1_fastest # y axis
    linprog =  data_ngroup_1$linprog # y axis
    
    max_ylim = max(max(data_ngroup_1$linprog),
                   max(data_ngroup_1$heuristic_1_closest),
                   max(data_ngroup_1$heuristic_1_closest))
    
    color_hclose = rgb(0,0,1) 
    color_hfast = rgb(0.2,0.4,0.1)
    color_linprog = rgb(1,0,0)
    
    # Make a basic graph
    plot( heuristic_1_closest ~ epochs , type = "b" , bty="l" , xlab="Total Epochs" , ylab=evaluator, main = titlees[i] , col=color_hclose, lwd=3 , pch=17, ylim=c(0,max_ylim) )
    lines( heuristic_1_fastest ~ epochs , col = color_hfast, lwd=3 , pch=19 , type="b" )
    lines( linprog ~ epochs, col=color_linprog , lwd=3 , pch=20 , type="b" ) # y~x
  }
  # plot.new()
  # par(xpd=TRUE)
  # Add a legend
  # par(mar=c(0,0,0,0))
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  # legend("bottom", c("IM", "IBD", "1R", "2R"), xpd = TRUE, horiz = TRUE,
  #        inset = c(0,0), bty = "n", pch = c(4, 2, 15, 19), col = 1:4, cex = 2)
  legend(x = "center",inset = 0,
         legend = c("Heur/ic: Closest", "Heur/ic: Fastest", "Optimal"),
         border = "black",
         col = c(color_hclose,
                 color_hfast,
                 color_linprog),
         pch = c(17,19,20),
         bty = "n",
         pt.cex = 1.3,
         cex = 1,
         text.col = "black",
         horiz = FALSE,
         lwd=5,xpd="NA")
  
  # Save plot as png:
  dev.copy(png,filename2)
  dev.off()
  
  #### Show results in tables:
  # install.packages('xtable')
  library('xtable')
  output_table = dataa[order(dataa$Node.Group),] # sort by node group because that's how data are shown in the plot above.
  latex_table = xtable(output_table,caption = paste(evaluator," results",sep=''),
                       label = paste("tbl:",evaluator," results table"), type = "latex")
  print(xtable(latex_table), file = filename3)
  ## factor(x = character(), levels, labels = levels,
  ##        exclude = NA, ordered = is.ordered(x), nmax = NA)
}
