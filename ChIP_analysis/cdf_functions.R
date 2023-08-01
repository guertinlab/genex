
cdf.deseq.df <- function(genes = gene.file, chip.peaks = chip.peaks, cat = "Repressed") {
  bed.tss.activated = get.tss(genes[genes$V5 == cat,])
  bed.tss.unchanged = get.tss(genes[genes$V5 == paste0("Matched to ", cat),])
  act.distance = bedTools.closest(bed1 = bed.tss.activated, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')
  unreg.distance = bedTools.closest(bed1 = bed.tss.unchanged, bed2 = chip.peaks[,c(1:3)], opt.string = '-D a')

  df.up.can = cbind(act.distance[,c(4, 10)], cat)
  df.un.can = cbind(unreg.distance[,c(4, 10)], paste0("Matched to ", cat))

  colnames(df.up.can) = c(colnames(df.up.can)[1:2], 'status')
  colnames(df.un.can) = c(colnames(df.up.can)[1:2], 'status')

  df.all = rbind(df.up.can, df.un.can)
  df.all$status = factor(df.all$status, levels = c(cat, paste0("Matched to ", cat)))
  return(df.all)
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

plot_cdf <- function(df.all, tf="quantile", cat = "Repressed", col.lines = c("#ce228e", "grey60", "#2290cf","grey90"), line.type = c(1)) {
pdf(paste0(tf, "_FIG_cdf_compare_Reg_classes_", cat, ".pdf"), width=6.2, height=3.83) 
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
