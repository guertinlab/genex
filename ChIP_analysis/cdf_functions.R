
cdf.deseq.df <- function(df, genes = gene.file, chip.peaks = chip.peaks) {
  bed.tss.activated = filter.deseq.into.bed(df, genes, cat = "Repressed")
  bed.tss.unchanged = filter.deseq.into.bed(df, genes, cat = "Matched to Repressed")
  act.distance = bedTools.closest(bed1 = bed.tss.activated, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')
  unreg.distance = bedTools.closest(bed1 = bed.tss.unchanged, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')

  df.up.can = cbind(act.distance[,c(4, 10)], "Repressed")
  df.un.can = cbind(unreg.distance[,c(4, 10)], "Matched to Repressed")

  colnames(df.up.can) = c(colnames(df.up.can)[1:2], 'status')
  colnames(df.un.can) = c(colnames(df.up.can)[1:2], 'status')

  df.all = rbind(df.up.can, df.un.can)
  df.all$status = factor(df.all$status, levels = c("Repressed", "Matched to Repressed"))
  return(df.all)
}

filter.deseq.into.bed <- function(deseq.df, gene.file, cat = 'R1881 Activated') {
  deseq.df = deseq.df[deseq.df$treatment == cat,]
  #print(dim(deseq.df))
  #scientific notation was messing this up occasionally
  x = gene.file$V4
  #print(length(x))
  y = gene.file[x %in% rownames(deseq.df),]
  #print(dim(y))
  z = get.tss(y)
  #print(dim(z))
  return(z)
}

bedTools.closest <- function(functionstring="/usr/local/bin/closestBed",bed1,bed2,opt.string="") {
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file= 'a.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file= 'b.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
  
  # create the command string and call the command using system()
  command1=paste('sort -k1,1 -k2,2n', 'a.file.bed', '> a.file.sorted.bed')
  cat(command1,"\n")
  try(system(command1))
  command2=paste('sort -k1,1 -k2,2n', 'b.file.bed', '> b.file.sorted.bed')
  cat(command2,"\n")
  try(system(command2))
  
  command=paste(functionstring, opt.string,"-a",'a.file.sorted.bed',"-b",'b.file.sorted.bed',">",'out.file.bed',sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table('out.file.bed',header=F, comment.char='')
  
  command3=paste('rm', 'a.file.bed', 'b.file.bed', 'a.file.sorted.bed', 'b.file.sorted.bed', 'out.file.bed')
  cat(command3,"\n")
  try(system(command3))
  
  colnames(res) = c(colnames(bed1), colnames(bed2), 'dis' )
  return(res)
}

filter.deseq.into.bed <- function(deseq.df, gene.file, cat = 'R1881 Activated') {
  deseq.df = deseq.df[deseq.df$treatment == cat,]
  print(dim(deseq.df))
  #scientific notation was messing this up occasionally
  x = gene.file$V4
  print(length(x))
  xy=sapply(strsplit(sapply(strsplit(rownames(deseq.df), ':'), '[', 2), "-"), "[", 2)
  gene = sapply(strsplit(xy, "_"), "[", 2)
  y = gene.file[x %in% gene,]
  #print(dim(y))
  z = get.tss(y)
  #print(dim(z))
  return(z)
}

plot_cdf <- function(df.all, tf="quantile", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), line.type = c(1)) {
pdf(paste0(tf, "FIG_cdf_compare_Reg_classes.pdf"), width=6.2, height=3.83) 
    print(ecdfplot(~log(abs(dis), base = 10), groups = status, data = df.all,
         auto.key = list(lines=TRUE, points=FALSE),
         col = col.lines,
         aspect = 1,
                                        #xlim = c(0, 50000),
         scales=list(relation="free",alternating=c(1,1,1,1)),
         ylab = 'Cumulative Distribution Function',
         xlab = expression('log'[10]~'ZNF143 Distance from TSS'),
                                        #index.cond = list(c(2,1)),
         between=list(y=1.0),
         type = 'a',
         xlim = c(0,8),
         lwd=2,
         lty=line.type,
         par.settings = list(superpose.line = list(col = col.lines, lwd=3), strip.background=list(col="grey85")),
         panel = function(...) {
             panel.abline(v= 200, lty =2)
             panel.ecdfplot(...)
         }))
    dev.off()
}
