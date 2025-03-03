x <- read.table(file="diffusion_events_annotated.txt", header=TRUE, sep="\t", quote="")
b4 <- x$biome4_names
b98 <- x$biome98_names
b98 <- gsub(" $", "", b98)
b4_u <- sort(unique(b4))
b98_u <- sort(unique(b98))
L1 <- length(b4_u)
L2 <- length(b98_u)
m <- matrix(0, nrow=L1, ncol=L2, dimnames=list(b4_u,b98_u))
for (i in 1:length(b4)) {
  p <- match(b4[i], rownames(m))
  q <- match(b98[i], colnames(m))
  m[p,q] <- m[p,q] + 1
}
abbr_row <- gsub("ropical", "rop.", rownames(m))
abbr_col <- gsub("ropical", "rop.", colnames(m))
rownames(m) <- abbr_row
colnames(m) <- abbr_col

heatmap(m, Rowv=NA, Colv=NA, margins=c(16,14), ylab="past", xlab="present", cexRow = .7, cexCol = .7)

