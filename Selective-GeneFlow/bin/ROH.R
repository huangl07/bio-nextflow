library(optparse)
options(bitmapType='cairo')
option_list <- list(
  make_option(c("-p", "--pop"), type = "character", action = "store", default = NULL,help = "Input keys word by pop"),
  make_option(c("-o", "--threshold"), type = "integer", action = "store",  default = 0.7 , help = "output file")

)
opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"))
if(is.null(opt$pop)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
if(is.null(opt$threshold)){opt = parse_args(OptionParser(option_list = option_list,usage = "Usage: %prog [options] \nDescription: This Script is used to draw manhattan!"),args="--help")}
times<-Sys.time()

library(detectRUNS)
slidingRuns<-slidingRUNS.run(genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."))
summaryList <- summaryRuns(runs = slidingRuns,genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."))
table<-tableRuns(slidingRuns,,genotypeFile=paste(opt$pop,"ped",sep="."),mapFile=paste(opt$pop,"map",sep="."),threshold=opt$threshold)
plot_manhattanRuns(savePlots=T,outputName="ROH",runs = slidingRuns,genotypeFile = paste(opt$pop,"ped",sep="."), mapFile = paste(opt$pop,"map",sep="."))
write.table(file="ROH.sliding",slidingRuns)
write.table(file="ROH.table",table)
write.table(file="summary_ROH_count_chr",summaryList$summary_ROH_count_chr)
write.table(file="summary_ROH_percentage_chr",summaryList$summary_ROH_percentage_chr)
write.table(file="summary_ROH_count",summaryList$summary_ROH_count)
write.table(file="summary_ROH_percentage",summaryList$summary_ROH_percentage)
write.table(file="summary_ROH_mean_chr",summaryList$summary_ROH_mean_chr)
write.table(file="summary_ROH_mean_class",summaryList$summary_ROH_mean_class)
write.table(file="result_Froh_genome_wide",summaryList$result_Froh_genome_wide)
write.table(file="result_Froh_chromosome_wide",summaryList$result_Froh_chromosome_wide)
write.table(file="result_Froh_class",summaryList$result_Froh_class)
write.table(file="SNPinRun",summaryList$SNPinRun)


escaptime=Sys.time()-times
print("Done!")
print(escaptime)
