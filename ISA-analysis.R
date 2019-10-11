#!/usr/bin/Rscript
library(limma)
library(pheatmap)
library(edgeR)
library(reshape2)
library(DRIMSeq)
library(SummarizedExperiment)
library(ggsignif)
library(gridExtra)

library(IsoformSwitchAnalyzeR)
library("BSgenome.Mmusculus.UCSC.mm10")
library(org.Mm.eg.db)
library(GenomicAlignments)
library(SummarizedExperiment)
library(dplyr)


get.chr <- function(x){
  gr <- grep(x, names(transcripts.seqs),value = T)
  if(identical(gr, character(0))){
    return('no')
  } else {
    rt <- return(strsplit(strsplit(gr, ' ')[[1]][3],':')[[1]][3])
  }
}

get.spl <- function(event){
  spl.sign <- spl[spl[[event]] > 0,]
  chr <- spl.sign$chr
  spl.sign <- spl.sign[,grep(paste('^', event, sep = ''), colnames(spl.sign))]
  spl.sign$chr <- chr
  return(spl.sign)
}


getSeq.mod <- function(crd){
  return(as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, strsplit(crd, ':')[[1]][1], start=as.numeric(strsplit(strsplit(crd, ':')[[1]][2],'-')[[1]][1]), end=as.numeric(strsplit(strsplit(crd, ':')[[1]][2],'-')[[1]][2]))))  
}

calc.gc <- function(seq){
  letters <- unlist(strsplit(seq, ''))
  return(100*(length(letters[letters %in% c('C', 'G')]) / length(letters)))
}

get.gc.length.event <- function(event){
  spl.sign <- get.spl(event)
  crds <- data.frame()
  cat('Building coordinates', sep = '\n')
  for(f in 1:nrow(spl.sign)){
    tmp <- spl.sign[f,]
    if(length(unlist(strsplit(tmp[[paste(event, '_genomic_start', sep = '')]], ';'))) > 1){
      ret <- data.frame(chr = paste('chr', tmp$chr, sep = ''), start = unlist(strsplit(tmp[[paste(event, '_genomic_start', sep = '')]], ';')), stop = unlist(strsplit(tmp[[paste(event, '_genomic_end', sep = '')]], ';')))
    } else {
      ret <- data.frame(chr = paste('chr', tmp$chr, sep = ''), start = unlist(strsplit(tmp[[paste(event, '_genomic_start', sep = '')]], ';')), stop = unlist(strsplit(tmp[[paste(event, '_genomic_end', sep = '')]], ';')))
    }
    crds <- rbind(ret,crds)
  }
  
  crds <- paste(crds$chr, ':', crds$start, '-', crds$stop, sep = '')
  cat('Extracting sequences', sep = '\n')
  seqs <- unlist(lapply(crds, getSeq.mod))
  cat('Calculating length and GC content', sep = '\n')
  seq.length <- unlist(lapply(seqs, function(x){nchar(x)}))
  seq.gc <- unlist(lapply(seqs, calc.gc))
  seq.info <- data.frame(len = seq.length, gc = seq.gc, event = event)
  return(seq.info)
}



transcripts.seqs <- readDNAStringSet('/media/alexander/Data/stringtie.results/Mus_musculus.GRCm38.cds.all.fa')
transcript.names <- unlist(lapply(strsplit(names(transcripts.seqs), ' '),function(x)x[[1]]))

##STRINGTIE LOAD
myQantifications <- importIsoformExpression(parentDir = '/media/alexander/Data/stringtie.results/stringtie/',
                                            pattern = 'SRR',readLength = 75,showProgress = T)

design <- read.table('/media/alexander/Data/stringtie.results/samples.txt',header = F,stringsAsFactors = F)
colnames(design) <- c('sampleID', 'condition')

cond <- data.frame(condition_1 = c('Tg1','Wt1','Wt1','Tg1'),
                   condition_2 = c('Tg2','Tg1','Wt3','Tg3'))

mySwitchList <- importRdata(isoformCountMatrix = myQantifications$abundance,
                            showProgress = T,
                            designMatrix = design,
                            isoformExonAnnoation = '/media/alexander/Data/stringtie.results/merged.annotated.gtf',
                            comparisonsToMake = cond)


mySwitchList <- preFilter(mySwitchList,geneExpressionCutoff = 0.01, IFcutoff = 0,quiet = F,isoformExpressionCutoff = 3)
#isoformSwitchAnalysisCombined(switchAnalyzeRlist = mySwitchList,n = 100,pathToOutput = '~/results/report/',fileType = 'pdf',genomeObject = BSgenome.Mmusculus.UCSC.mm10,outputPlots = T,switchTestMethod = 'DEXSeq')

mySwitchList <- isoformSwitchTestDEXSeq(mySwitchList,showProgress = T,correctForConfoundingFactors = T)


#mySwitchList <- isoformSwitchTestDRIMSeq(switchAnalyzeRlist = mySwitchList,showProgress = T)
mySwitchList <- analyzeORF(switchAnalyzeRlist = mySwitchList,genomeObject = BSgenome.Mmusculus.UCSC.mm10)
mySwitchList <- analyzeAlternativeSplicing(switchAnalyzeRlist = mySwitchList)

mySwitchList <- analyzeCPAT(switchAnalyzeRlist = mySwitchList,pathToCPATresultFile = '/media/alexander/Data/stringtie.results/supplementary/CPAT.txt',codingCutoff = 0.723,removeNoncodinORFs = F)

save.image('~/results/all.Rdata')
load('~/results/all.Rdata')

###transcripts annotation?
#rt <- data.frame(rtracklayer::import('~/transcriptomes/splcing.detection.comparsion/pipeline/merged/merged.annotated.gtf'))
#mySwitchList$isoformSwitchAnalysis$isoform_id <- rt[match(mySwitchList$isoformSwitchAnalysis$isoform_id, rt$transcript_id),]$cmp_ref
#mySwitchList$AlternativeSplicingAnalysis$isoform_id <- rt[match(mySwitchList$AlternativeSplicingAnalysis$isoform_id, rt$transcript_id),]$cmp_ref


###ALTSPLICE
sw <- data.frame(mySwitchList$isoformSwitchAnalysis)
sw$test <- paste(sw$condition_2,'_', sw$condition_1,sep = '')

spl <- data.frame(mySwitchList$AlternativeSplicingAnalysis)
spl.ids <- spl$isoform_id
spl <- spl[,nchar(colnames(spl)) < 5]
rownames(spl) <- spl.ids

sw$test <- paste(sw$condition_2,sw$condition_1, sep = '_')
ggplot(data=sw) + geom_boxplot(aes(x = test, y = dIF)) 
ggplot(data=sw[sw$isoform_id %in% rownames(spl[spl$IR > 0,]),]) + geom_boxplot(aes(x = test, y = dIF))
ggplot(data=sw[sw$isoform_id %in% rownames(spl[spl$ES > 0,]),]) + geom_boxplot(aes(x = test, y = dIF))

extractSplicingEnrichment(mySwitchList)

spl.summary <- data.frame()
for(f in colnames(spl)){
  print(f)
  tmp <- sw[sw$isoform_id %in% rownames(spl[spl[[f]] > 0,]),]
  ret <- data.frame(100*(table(t(tmp$test)) / sum(table(t(tmp$test)))))
  if(nrow(ret) == 0){
    next
  }
  ret$id <- f
  spl.summary <- rbind(ret,spl.summary)
}


spl <- data.frame(mySwitchList$AlternativeSplicingAnalysis)
spl <- spl[grep('ENS', spl$isoform_id),]

spl$chr <- unlist(lapply(spl$isoform_id, get.chr))
spl <- spl[spl$chr != 'no',]


event <- 'ES'


info.es <- get.gc.length.event('ES')
info.ri <- get.gc.length.event('IR')
info.atss <- get.gc.length.event('ATSS')
info.atts <- get.gc.length.event('ATTS')

gcs <- bind_rows(info.es, info.ri,info.atts, info.atss)

ggplot(data=gcs) + geom_boxplot(aes(x = event, y= gc))
ggplot(data=gcs) + geom_boxplot(aes(x = event, y= len))
ggplot(data=gcs) + geom_point(aes(x = gc, y = len, color = event)) + theme_bw()

###splicing gene ontology
###intron retention
spl$rs <- rowSums(spl)
spl.src <- data.frame(mySwitchList$AlternativeSplicingAnalysis)

###exon skipping
cat(grep('ENS', spl.src[spl.src$ES > 0,][,1:4]$isoform_id,value = T), sep = '\n')


###SWITCH ANALYSIS
ggplot(data=sw) + geom_point(aes(x = dIF, y = -log10(padj), color = test, size = 3)) + ylim(c(0,30)) + theme_bw()


