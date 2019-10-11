library(SAJR)
mazin.ref <- read.csv('~/transcriptomes/splcing.detection.comparsion/as.events.csv',stringsAsFactors = F)
setwd("~/transcriptomes/splcing.detection.comparsion/data/sajr/")

data = loadSAData(ann.gff='~/transcriptomes/splcing.detection.comparsion/soft/sajr/ann.gtf',c(1,2,3,4,5,6,7,8,9,10))
data = setSplSiteTypes(data,'~/transcriptomes/splcing.detection.comparsion/soft/sajr/ann.gtf')


data.f = data[data$seg$type %in% c('ALT','INT') & data$seg$position %in% c('LAST','INTERNAL','FIRST') & apply(data$i+data$e>=10,1,sum)==2 & apply(data$ir,1,sd,na.rm=TRUE) > 0,]
all.alts = makeAlts(data$seg,'~/transcriptomes/splcing.detection.comparsion/soft/sajr/ann.gtf',remove.exn.ext = F)

mod = list(f=factor(c(rep('Tg1', 5), rep('Tg3',5))))
data.f.glm = fitSAGLM(data.f,terms(x ~ f),mod)
data.f.pv = calcSAPvalue(data.f.glm)
data.f.pv[,2] = p.adjust(data.f.pv[,2],method='BH')


###DPSI CALC???

data.sign = data.f[data.f.pv[,2] <= 0.05,]
length(data.sign)

data.sign[order(abs(data.sign$ir[,1]-data.sign$ir[,2]),decreasing = TRUE),][1:18,]
data.sign$seg

events <- data.frame(code = c('ad', 'aa', 'dd', 'da'), 
                     val = c('CE', 'AA', 'AD', 'RI'))

events[match(data.sign$seg$sites, events$code),]$val


intersect(mazin.ref$gene_id, unlist(lapply(strsplit(rownames(data.sign$ir),'.', fixed = T),function(x)x[[1]])))
rowMeans(data.f$ir[,seq(6:10)]) - rowMeans(data.f$ir[,seq(1:5)])

tt <- as.character(rowMeans(data.f$ir[,seq(6:10)]) - rowMeans(data.f$ir[,seq(1:5)]))
tt <- tt[complete.cases(tt)]
hist(tt)

