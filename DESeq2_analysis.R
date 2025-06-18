library(DESeq2)

lb_data <- read.table("Lb_all_counts_table.txt", header=1)
as_data <- read.table("As_all_counts_table.txt", header=1)[,c(-1,-28,-27,-26)]

full_data <- cbind(lb_data, as_data)
count_data <- full_data[,-1]
rownames(count_data) <- full_data[,1]

coldata <- data.frame(cas = c(rep("Lb",24),rep("As",24)))
coldata$cas <- factor(coldata$cas)

coldata$effector = c(rep("KRAB", 12), rep("SID", 12)) 
coldata$effector <- factor(coldata$effector)

coldata$day = c(rep(5,3), rep(7,3), rep(14,3), rep(21,3))
coldata$day = coldata$day - 5

lb_krab_count_data <- count_data[,1:12]
lb_sid_count_data <- count_data[,12+(1:12)]
as_krab_count_data <- count_data[,24+(1:12)]
as_sid_count_data <- count_data[,36+(1:12)]

lb_krab_col_data <- coldata[1:12,]
lb_sid_col_data <- coldata[12+1:12,]
as_krab_col_data <- coldata[24+1:12,]
as_sid_col_data <- coldata[36+1:12,]

full_results <- DESeq(DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ day*effector*cas))

lb_krab_dds <- DESeq(DESeqDataSetFromMatrix(countData = lb_krab_count_data, colData = lb_krab_col_data, design = ~ day))
lb_krab_dds <- DESeq(lb_krab_dds, test="LRT", reduced = ~1)
lb_krab_res <- results(lb_krab_dds)

lb_sid_dds <- DESeq(DESeqDataSetFromMatrix(countData = lb_sid_count_data, colData = lb_sid_col_data, design = ~ day))
lb_sid_dds <- DESeq(lb_sid_dds, test="LRT", reduced = ~1)
lb_sid_res <- results(lb_sid_dds)

as_krab_dds <- DESeq(DESeqDataSetFromMatrix(countData = as_krab_count_data, colData = as_krab_col_data, design = ~ day))
as_krab_dds <- DESeq(as_krab_dds, test="LRT", reduced = ~1)
as_krab_res <- results(as_krab_dds)

as_sid_dds <- DESeq(DESeqDataSetFromMatrix(countData = as_sid_count_data, colData = as_sid_col_data, design = ~ day))
as_sid_dds <- DESeq(as_sid_dds, test="LRT", reduced = ~1)
as_sid_res <- results(as_sid_dds)

write.csv(file="normalized_counts.csv", counts(full_results, normalized=TRUE))

write.csv(file="lb_krab_sig.csv", lb_krab_res[which(lb_krab_res$padj < 0.001),])
write.csv(file="lb_sid_sig.csv", lb_sid_res[which(lb_sid_res$padj < 0.001),])
write.csv(file="as_krab_sig.csv", as_krab_res[which(as_krab_res$padj < 0.001),])
write.csv(file="as_sid_sig.csv", as_sid_res[which(as_sid_res$padj < 0.001),])

write.csv(file="lb_krab_all.csv", lb_krab_res)
write.csv(file="lb_sid_all.csv", lb_sid_res)
write.csv(file="as_krab_all.csv", as_krab_res)
write.csv(file="as_sid_all.csv", as_sid_res)

plot_grna <- function(dds, grna_name) {
   ggplot(data = plotCounts(dds, gene=grna_name, intgroup=c("day","cas","effector"), returnData=TRUE), aes(x=day, y=count, group=day, color=cas, shape=effector)) + 			stat_summary(fun=mean, geom="line", aes(group = factor(effector:cas))) + 
	geom_point() + 
	labs(title=grna_name)
}