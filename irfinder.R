library(DESeq2)
source("~/transcriptomes/splcing.detection.comparsion/soft/IRFinder/bin/DESeq2Constructor.R")  #Load IRFinder-related function
setwd('~/transcriptomes/splcing.detection.comparsion/data/irfinder.res/')

ref <- read.csv('~/transcriptomes/splcing.detection.comparsion/as.events.csv',stringsAsFactors = F)

results = read.table("filePaths.txt")
paths = as.vector(results$V1)                                            # File names must be saved in a vector
experiment = read.table("experiment.txt",header=T)                       
experiment$Condition=factor(experiment$Condition,levels=c("Tg1","Tg3"))    # Set WT as the baseline in the analysis
rownames(experiment)=NULL                                                # Force removing rownames

# WARNING: make sure the rownames of `experiment` is set to NULL. 
# WARNING: users MUST check if the order of files in the `path` matches the order of samples in `experiment` before continue  

metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)



# The above line generate a meta list contains three slots
# First slot is a DESeq2 Object that can be directly pass to DESeq2 analysis.  
# Second slot is a matrix for trimmed mean of intron depth
# Third slot  is a matrix for correct splicing depth flanking introns
# We build a “null” regression model on the interception only. 
# A “real” model can be assigned either here directly, or in the downstream step. See below

dds = metaList$DESeq2Object                       # Extract DESeq2 Object with normalization factors ready
colData(dds)                                      # Check design of matrix


# Please note that sample size has been doubled and one additional column "IRFinder" has been added.
# This is because IRFinder considers each sample has two sets of counts: one for reads inside intronic region and one for reads at splice site, indicating by "IR" and "Splice" respectively.
# "IRFinder" is considered as an additional variable in the GLM model.
# Please also be aware that size factors have been set to 1 for all samples. Re-estimation of size factors is NOT recommended and is going to bias the result.
# More details at the end of the instruction.

design(dds) = ~Condition + Condition:IRFinder     # Build a formula of GLM. Read below for more details. 
dds = DESeq(dds)                                  # Estimate parameters and fit to model


res.Tg1 = results(dds, name = "ConditionTg1.IRFinderIR")

# This tests if the number of IR reads are significantly different from normal spliced reads, in the WT samples.
# We might only be interested in the "log2FoldChange" column, instead of the significance.
# This is because "log2FoldChange" represents log2(number of intronic reads/number of normal spliced reads).
# So we the value of (intronic reads/normal spliced reads) by

Tg1.IR_vs_Splice=2^res.Tg1$log2FoldChange

# As IR ratio is calculated as (intronic reads/(intronic reads+normal spliced reads))
# We can easily convert the above value to IR ratio by

IRratio.Tg1 = Tg1.IR_vs_Splice/(1+Tg1.IR_vs_Splice)

# Similarly, we can get IR ratio in the KO samples
res.Tg3 = results(dds, name = "ConditionTg3.IRFinderIR")
Tg3.IR_vs_Splice=2^res.Tg3$log2FoldChange
IRratio.Tg3 = Tg3.IR_vs_Splice/(1+Tg3.IR_vs_Splice)

# Finally we can test the difference of (intronic reads/normal spliced reads) ratio between WT and KO
res.diff = results(dds, contrast=list("ConditionTg1.IRFinderIR","ConditionTg3.IRFinderIR"))  

# We can plot the changes of IR ratio with p values
# In this example we defined significant IR changes as
# 1) IR changes no less than 10% (both direction) and 
# 2) with adjusted p values less than 0.05

IR.change = IRratio.Tg3 - IRratio.Tg1
#plot(IR.change,col=ifelse(res.diff$padj < 0.05 & abs(IR.change)>=0.1, "red", "black"))

res.diff <- data.frame(res.diff)
res.diff <- data.frame(res.diff)

res.diff <- res.diff[complete.cases(res.diff),]
res.diff <- res.diff[res.diff$padj < 0.05,]
View(res.diff)
#res.diff <- res.diff[res.diff$pvalue < 0.05 & res.diff$padj < 0.1,]

write.table(res.diff, '~/transcriptomes/splcing.detection.comparsion/results/irfinder.txt')


intersect(unlist(lapply(strsplit(rownames(res.diff), '/'),function(x)x[[2]])), ref[ref$as.type == 'RI',]$gene_id)




