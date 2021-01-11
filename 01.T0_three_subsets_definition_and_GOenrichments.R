#Define the working directory based on which computer working from
EVO_IMMUNO_POP="S:/"
EVO_IMMUNO_POP="/pasteur/projets/policy01/evo_immuno_pop/"
EVO_IMMUNO_POP="/Volumes/evo_immuno_pop/"

#Specify library path if on server and load packages
.libPaths('/pasteur/projets/policy01/evo_immuno_pop/R_3.6_scLib')

library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)

#Set the working directory
setwd(EVO_IMMUNO_POP)

#Load the data
keepMeta=readRDS("IAV/data/meta_dataframe_88559monocytes.rds")
mat=readRDS("IAV/data/normalized_count_matrix_88559monocytes_22593hostgenes.rds")

facs <- read.table(paste0(EVO_IMMUNO_POP, "IAV/FACSdata_withIAVpct_and_cellStates_199ind.txt"), header=T, sep="\t")

#Extract T0
t0 <- filter(keepMeta, TP==0)

#Extract the 8 donors info for FACS
tmp <- filter(facs, IID %in% levels(keepMeta$ID1))   

#tmp <- filter(facs, IID %in% c("EUB058", "EUB078", "AFB108", "EUB061", "AFB022", "AFB131", "AFB176", "EUB060"))
tmp$perNC <- (tmp$CD16pos_CD14pos-tmp$CD16pos_CD14hi_Q5)/(tmp$CD16pos_CD14pos+tmp$CD16neg_CD14pos)
tmp$perINT <- (tmp$CD16pos_CD14hi_Q5)/(tmp$CD16pos_CD14pos+tmp$CD16neg_CD14pos)
tmp$perCL <- (tmp$CD16neg_CD14pos)/(tmp$CD16pos_CD14pos+tmp$CD16neg_CD14pos)

facs$perNC <- (facs$CD16pos_CD14pos-facs$CD16pos_CD14hi_Q5)/(facs$CD16pos_CD14pos+facs$CD16neg_CD14pos)
facs$perINT <- (facs$CD16pos_CD14hi_Q5)/(facs$CD16pos_CD14pos+facs$CD16neg_CD14pos)
facs$perCL <- (facs$CD16neg_CD14pos)/(facs$CD16pos_CD14pos+facs$CD16neg_CD14pos)
facs$perCD16pos <- (facs$CD16pos_CD14pos)/(facs$CD16pos_CD14pos+facs$CD16neg_CD14pos)
facs$ratio <- (facs$CD16pos_CD14pos)/(facs$CD16neg_CD14pos)

mean(facs$perNC, na.rm=TRUE)
mean(facs$perINT, na.rm=TRUE)
mean(facs$perCL, na.rm=TRUE)

cnts <- data.frame(table(t0$ID1))
tmp$cnt <- cnts[match(tmp$IID, cnts$Var1), "Freq"]

q1 <- weighted.mean(tmp$perNC, tmp$cnt, na.rm=TRUE)
q2 <- weighted.mean(tmp$perINT, tmp$cnt, na.rm=TRUE)
q3 <- weighted.mean(tmp$perCL, tmp$cnt, na.rm=TRUE)

quantile(t0$t0PC1, probs=c(q1, (q1+q2)))

th1 <- unname(quantile(t0$t0PC1, probs=c(q1))) 
th2 <- unname(quantile(t0$t0PC1, probs=c((q1+q2))))

t0$class <- as.factor(cut(t0$t0PC1, breaks=c((min(t0$t0PC1)-1), th1, th2, (max(t0$t0PC1)+1))))
t0$class <- factor(t0$class, levels=levels(t0$class), labels=c("nonclassical", "intermediate", "classical"))
table(t0$class)

p1 <- ggplot(t0, aes(t0PC1)) +
  geom_histogram(aes(fill=sxCD16), bins=50) +
  theme_light() +
  scale_fill_manual(values=c("#999999", "#333333")) +
  theme(panel.grid = element_blank(),
        legend.position="none") +
  xlab("") +
  geom_vline(aes(xintercept = th1), linetype="dotted", color="red") +
  geom_vline(aes(xintercept = th2), linetype="dotted", color="red") +
  geom_text(aes(x=-30, y=550), label="CD16")

p2 <- ggplot(t0, aes(t0PC1)) +
  geom_histogram(aes(fill=sxCD14), bins=50) +
  theme_light() +
  scale_fill_manual(values=c("#999999", "#333333")) +
  theme(panel.grid = element_blank(),
        legend.position="none") +
  geom_vline(aes(xintercept = th1), linetype="dotted", color="red") +
  geom_vline(aes(xintercept = th2), linetype="dotted", color="red") +
  geom_text(aes(x=-30, y=550), label="CD14")

partA1 <- plot_grid(p1, p2, align="hv", axis="tblr", nrow=2)

partA2 <- ggplot(t0, aes(t0umap1, t0umap2)) +
  geom_point(aes(color=class)) +
  theme_light() +
  theme(panel.grid = element_blank(),
        legend.position=c(0.1,0.25)) +
  scale_color_viridis(direction=-1, discrete=T)

partA <- plot_grid(partA1, partA2, nrow=1)


comp <- t0 %>% group_by(ID1) %>% dplyr::count(class) %>% data.frame()
comp$total <- cnts[match(comp$ID1, cnts$Var1),'Freq']
comp$per <- comp$n/comp$total

library(reshape2)
tmpM <- tmp[,c(2,23:25)]
names(tmpM) <- c("IID", "nonclassical", "intermediate", "classical")
tmpM <- melt(tmpM, id.vars=c("IID"))
tmpM$merge <- paste(tmpM$IID, tmpM$variable, sep="_")
comp$merge <- paste(comp$ID1, comp$class, sep="_")
comp$facs <- tmpM[match(paste0(comp$ID1, comp$class), paste0(tmpM$IID, tmpM$variable)), 'value']

cor.test(comp[comp$class=="nonclassical", 'facs'], comp[comp$class=="nonclassical", 'per'])
cor.test(comp[comp$class=="classical", 'facs'], comp[comp$class=="classical", 'per'])
cor.test(comp[comp$class=="intermediate", 'facs'], comp[comp$class=="intermediate", 'per'])

tt <- t0 %>% group_by(ID1) %>% dplyr::count(sxCD16) %>% dcast(ID1~sxCD16) %>% data.frame()
tt$per <- tt$FALSE./(tt$FALSE.+tt$TRUE.)
tt$ratio <- tt$FALSE./(tt$TRUE.)
tt$facs

partB <- ggplot(comp, aes(x=facs, y=per)) +
  facet_wrap(~class, scales="free") +
  geom_point(aes(color=ID1)) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_light() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.position = "none") +
  scale_color_viridis(direction=-1, option="C", discrete=T, end=0.9) +
  geom_abline(slope=1, intercept=0) +
  ylab("% Monocytes (scRNA-seq)") +
  xlab("% CD14+ Cells (FACS)")

plot_grid(partA, partB, nrow=2)


keepMeta$T0subset <- t0[match(keepMeta$Barcode, t0$Barcode), 'class']
pp <-  arrange(keepMeta, !is.na(T0subset), T0subset)
ggplot(pp, aes(tsne1, tsne2)) +
  geom_point(aes(color=T0subset)) +
  theme_light() +
  theme(panel.grid=element_blank())

#Plot some genes 
identical(t0$CD14, unname(mat["CD14",t0$Barcode]))
t0$IFITM3 = unname(mat["IFITM3",t0$Barcode])
t0$HLADRA = unname(mat["HLA-DRA",t0$Barcode])

t0Plot <- melt(t0[,c('Barcode', 'ID1', 'class', 'CD14', 'FCGR3A', 'HLADRA')], id.vars=c('Barcode', 'ID1', 'class'))

partC <- ggplot(t0Plot, aes(x=class, y=value)) +
  facet_wrap(~variable) +
  geom_violin(aes(fill=class), scale="width") +
  geom_boxplot(fill="white", alpha=0.3, notch=T, width=0.5) +
  theme_light() +
  theme(panel.grid=element_blank(),
        legend.position="none") +
  scale_color_viridis(direction=-1, option="C", discrete=T, end=0.9) +
  scale_fill_viridis(discrete=T, direction=-1)

intPlot <- plot_grid(partA, partB, partC, nrow=3)

#
compDF <- dcast(comp[,c(1, 2, 5)], ID1 ~ class)
flu <- read.table(paste0(EVO_IMMUNO_POP, "sc_resources/important_files/EvoImmunoPop_perFlu_estimates.txt"), header=T)
compDF$EIP <- flu[match(paste(compDF$ID1, "5", sep="-"), flu$sample), 'perUniq']

summary(lm(EIP ~ classical, compDF))              
summary(lm(EIP ~ nonclassical, compDF))              
summary(lm(EIP ~ intermediate, compDF))              

cor.test(compDF$classical, compDF$EIP)
cor.test(compDF$nonclassical, compDF$EIP)
cor.test(compDF$intermediate, compDF$EIP)

cor.test(tt$per, compDF$EIP)
cor.test(tt$ratio, compDF$EIP)

sx <- t0 %>% group_by(ID1) %>% dplyr::count(sxCD16) %>% data.frame() %>% dcast(ID1 ~ sxCD16)
sx$per <- sx$`TRUE` / rowSums(sx[,2:3])
sx$ratio <- sx$`TRUE` / sx$`FALSE`


facs$EIP <- flu[match(facs$sample, flu$sample), 'perUniq']

summary(lm(EIP ~ perCD16pos + population, facs))   
summary(lm(EIP ~ ratio + population, facs))   
summary(lm(EIP ~ perNC + population, facs))   
summary(lm(EIP ~ perCL + population, facs))   
summary(lm(EIP ~ perINT + population, facs))   

#######
######
######

####
library(scran)

#print the PC1-based calls
table(t0$class)

#Isolate t0 gene
t0mat <- mat[,t0$Barcode]
dim(t0mat)

# calculate average gene expression per subset
tmpList = list()
for (i in levels(t0$class)){
  avg <- unname(rowMeans(t0mat[, t0$class==i]))
  tmpList[[i]] = avg
}
avgSubsets <- data.frame(tmpList)
avgSubsets$gene <- rownames(t0mat)

#isolate the genes with average expression > 0.1 in at least one subset
t0keep <- avgSubsets[avgSubsets$nonclassical > 0.1 | avgSubsets$intermediate > 0.1 | avgSubsets$classical > 0.1, 'gene']
length(t0keep)

inc <- filter(t0, class %in% c("nonclassical", "intermediate"))
inc2 <- filter(t0, class %in% c("classical", "intermediate"))

exp <- findMarkers(t0mat[t0keep,inc$Barcode], as.factor(as.character(inc$class)), test.type=c("wilcox"), block=inc$ID1)
exp2 <- findMarkers(t0mat[t0keep,inc2$Barcode], as.factor(as.character(inc2$class)), test.type=c("wilcox"), block=inc2$ID1)
expT <- findMarkers(t0mat[t0keep,inc$Barcode], as.factor(as.character(inc$class)), test.type=c("t"), block=inc$ID1)
exp2T <- findMarkers(t0mat[t0keep,inc2$Barcode], as.factor(as.character(inc2$class)), test.type=c("t"), block=inc2$ID1)

exp[["intermediate"]]$gene <- rownames(exp[["intermediate"]])
exp2[["intermediate"]]$gene <- rownames(exp2[["intermediate"]])
expT[["intermediate"]]$gene <- rownames(expT[["intermediate"]])
exp2T[["intermediate"]]$gene <- rownames(exp2T[["intermediate"]])

intGeneDiff <- left_join(data.frame(exp[["intermediate"]][,c(5, 2:4)]), data.frame(expT[["intermediate"]][,4:5]), by="gene")
intGeneDiff <- left_join(intGeneDiff, data.frame(exp2[["intermediate"]][,2:5]), by="gene")
intGeneDiff <- left_join(intGeneDiff, data.frame(exp2T[["intermediate"]][,4:5]), by="gene")
intGeneDiff$lab <- abs(intGeneDiff$logFC.classical+intGeneDiff$logFC.nonclassical)
intGeneDiff$labText <- ifelse(intGeneDiff$lab > 1, intGeneDiff$gene, NA)


ggplot(intGeneDiff, aes(x=logFC.nonclassical, y=logFC.classical)) +
  geom_point() +
  theme_light() +
  theme(panel.grid=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  geom_text(aes(label=labText))

ot <- filter(t0, class %in% c("nonclassical", "classical"))
NCC <- findMarkers(t0mat[t0keep,ot$Barcode], as.factor(as.character(ot$class)), test.type=c("wilcox"), block=ot$ID1)
NCCt <- findMarkers(t0mat[t0keep,ot$Barcode], as.factor(as.character(ot$class)), test.type=c("t"), block=ot$ID1)

NCC <- data.frame(NCC[[1]])
NCCt <- data.frame(NCCt[[1]])
NCC$gene <- rownames(NCC)
NCCt$gene <- rownames(NCCt)

NCC <- left_join(NCC, NCCt, by="gene")
NCC$direction <- ifelse(NCC$FDR.x >= 0.01 | abs(NCC$logFC.nonclassical) < 0.2, 'nonsig',
                        ifelse(NCC$logFC.nonclassical > 0, "classical", "nonclassical"))
table(NCC$direction)

intGeneDiff$classic <- ifelse(intGeneDiff$FDR.y >= 0.01 | abs(intGeneDiff$logFC.classical) < 0.2, 'nonsig',
                              ifelse(intGeneDiff$logFC.classical < 0, "down", "up"))
intGeneDiff$nonclassic <- ifelse(intGeneDiff$FDR.x >= 0.01 | abs(intGeneDiff$logFC.nonclassical < 0.2), 'nonsig',
                                 ifelse(intGeneDiff$logFC.nonclassical < 0, "down", "up"))
table(intGeneDiff$classic)
table(intGeneDiff$nonclassic)

#load the gene annotation
GA=readRDS("IAV/data/gene_annotation_22593hostgenes.rds")
dim(GA)

intGeneDiff$ens <- GA[match(intGeneDiff$gene, GA$geneAlt), 'stable']
which( is.na( intGeneDiff$ens ) )

NCC$ens <- GA[match(NCC$gene, GA$geneAlt), 'stable']
which( is.na( NCC$ens ) )

nonclassicList <- NCC[NCC$direction=="nonclassical", 'ens']
classicList <- NCC[NCC$direction=="classical", 'ens']
intUP <- intGeneDiff[intGeneDiff$classic=="up" & intGeneDiff$nonclassic=="up", 'ens']
intDOWN <- intGeneDiff[intGeneDiff$classic=="down" & intGeneDiff$nonclassic=="down", 'ens']

#Alternatively from t0.master
nonclassicList <- t0.master[t0.master$upregulated_in=="nonclassical", 'ens']
classicList <- t0.master[t0.master$upregulated_in=="classical", 'ens']
intUP <- t0.master[t0.master$relative_classic=="up" & t0.master$relative_nonclassic=="up", 'ens']
intDOWN <- t0.master[t0.master$relative_classic=="down" & t0.master$relative_nonclassic=="down", 'ens']
background = t0.master$ens



source(paste0(EVO_IMMUNO_POP, "/sc_resources/GOSeq_hg38.R"))

background = NCC$ens

enrLists = list(nonclassicList, classicList, intUP)
names(enrLists) <- c("nonclassicList", "classicList", "intUP")

paste("Performing enrichment")
subGO=list()
for (y in 1:length(enrLists)){
  t=unique(enrLists[[y]])
  goRES <- GOSeq(t, background, addCI=TRUE)
  if (length(goRES$category) > 0) {
    goRES$genes <- NA
    tmpGA <- dplyr::filter(GA, stable %in% t)
    for (i in 1:length(goRES$category)) {
      category=goRES$category[i]
      genes <- tmpGA[grepl(category, tmpGA$GO),'Symbol']
      genes <- paste(genes, collapse=",")
      goRES[i, 'genes'] = genes
    }
  }
  nm=names(enrLists)[y]
  subGO[[nm]]=goRES
}

enr <- rbindlist(subGO, idcol=T, fill=T)
write.table(enr[,-10], paste0(EVO_IMMUNO_POP, "IAV/results/T0_three_subsets_GO.enrichments_4589backgroundgenes.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

t0.master <- left_join(NCC[c(11, 5, 2:4,9:10)], intGeneDiff[,c(1:9, 12:13)], by="gene")
names(t0.master) <- c('ens', 'gene', 'classical_nonclassical_pval', 'classical_nonclassical_fdr', 'classical_nonclassical_AUC', 'classical_nonclassical_logFC', 'upregulated_in', 'intermediate_nonclassical_pval', 'intermediate_nonclassical_fdr', 'intermediate_nonclassical_AUC', 'intermediate_nonclassical_logFC', 'intermediate_classical_pval', 'intermediate_classical_fdr', 'intermediate_classical_AUC', 'intermediate_classical_logFC', 'relative_classic', 'relative_nonclassic')

write.table(t0.master, file="IAV/results/T0_three_subsets_stats_4589genes.txt", row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))


intUP <- t0.master[t0.master$relative_classic=="up" & t0.master$relative_nonclassic=="up", 'ens']


write.table(enrLists[["nonclassicList"]], file="IAV/results/gene_lists/t0_three_subsets/nonclassical_347genes.txt", row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(enrLists[["classicList"]], file="IAV/results/gene_lists/t0_three_subsets/classical_501genes.txt", row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(enrLists[["intUP"]], file="IAV/results/gene_lists/t0_three_subsets/intermediate_18genes.txt", row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(background, file="IAV/results/gene_lists/t0_three_subsets/background_4589genes.txt", row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))

