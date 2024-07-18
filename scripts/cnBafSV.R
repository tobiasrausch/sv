library(ggplot2)
library(scales)
library(gridExtra)
library(DNAcopy)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=T)
x$chr = factor(x$chr, levels=c("chr1", "chr5"))

# Segmentation
seg=segments.summary(segment(smooth.CNA(CNA(log(x[,6]), x$chr, x$start, data.type="logratio", sampleid="tumor")), undo.splits="sdundo", undo.SD=1.5))
seg$chr = factor(seg$chr, levels=c("chr1", "chr5"))

# Read-depth
p = ggplot(data=x, aes(x=start, y=x[,6]))
p = p + geom_point(pch=21, size=0.5)
p = p + xlab("Chromosome")
p = p + ylab("Copy-number")
p = p + scale_x_continuous(labels=comma)
p = p + scale_y_continuous(labels=comma, breaks=0:8, limits=c(0,8))
p = p + geom_segment(data=seg, aes(x=loc.start, y=exp(seg.median), xend=loc.end, yend=exp(seg.median)), colour="darkorange")
p = p + theme(axis.text.x = element_text(angle=45, hjust=1))
p = p + facet_wrap(~chr, scales="free_x")
ggsave("cov.png", p, width=10, height=7)
if (length(args) == 1) { quit(); }

# Read-depth + SVs
sv = read.table(args[2], header=F)
colnames(sv) = c("chr","start","end","type","id", "chr2", "pos2")
p2 = p + geom_curve(data=sv, aes(x=start, xend=end, col=type), y=6, yend=6, curvature=-1)
sv = sv[sv$type=="BND",]
sv$pos2 = as.numeric(sv$pos2)
sv$chr = factor(sv$chr, levels=c("chr1", "chr5"))
sv$chr2 = factor(sv$chr2, levels=c("chr1", "chr5"))
p2 = p2 + geom_segment(data=sv, aes(x=start, xend=start, y=6, yend=8, col=type))
sv$chr = sv$chr2
sv$start = sv$pos2
p2 = p2 + geom_segment(data=sv, aes(x=start, xend=start, y=6, yend=8, col=type))
p2 = p2 + labs(colour="SV type")
p2 = p2 + theme(legend.position="bottom")
ggsave("cov.png", p2, width=10, height=7)
if (length(args) == 2) { quit(); }

# Read-depth + SVs + BAF
baf = read.table(args[3], header=F)
colnames(baf) = c("chr", "pos", "baf")
baf$chr = factor(baf$chr, levels=c("chr1", "chr5"))
q = ggplot(data=baf)
q = q + geom_jitter(aes(x=pos, y=baf), fill="red", colour="red", size=0.1, shape=21, width=0, height=0.05)
q = q + geom_jitter(aes(x=pos, y=1-baf), fill="blue", colour="blue", size=0.1, shape=21, width=0, height=0.05)
q = q + ylab("B-Allele Frequency") + xlab("Chromosome")
q = q + scale_x_continuous(labels=comma)
q = q + scale_y_continuous(labels=comma, breaks=c(0, 0.5, 1), limits=c(0, 1))
q = q + theme(axis.text.x = element_text(angle=45, hjust=1))
q = q + facet_wrap(~chr, scales="free_x")
plot_grid(q, p2, align="v", nrow=2, rel_heights=c(1/2, 1/2))
ggsave(file="cov.png", width=10, height=7)
print(warnings())
