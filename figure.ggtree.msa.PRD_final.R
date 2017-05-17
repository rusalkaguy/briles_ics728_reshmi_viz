#
# ggtree
#
# https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeImport.html

library("ggplot2")
library(ggtree) # source("https://bioconductor.org/biocLite.R"); biocLite("ggtree")
library(phylobase) # source("https://bioconductor.org/biocLite.R");  biocLite("phylobase")
library(Biostrings)

# -----------------------------------------------------------------------------------
# load sample -> GROUP annotation
#

brilesGroups = read.delim("../../v2/CODING_MOTIFS_ICS728.161121.groups.txt", sep="\t", header=TRUE, stringsAsFactors=F)
#bg = brilesGroups[,c(-1,-2)]
#row.names(bg) = brilesGroups$name
bg = data.frame(row.names=brilesGroups$name, grp=brilesGroups$ColorCode)
write.table(rownames(bg), "sample_names.CODING_MOTIFS_ICS728.161121.groups.txt")

# ================================================================= 
# load Reshmi's NJ tree (PRD Final.nex)
# had to convert from .nex to .fa with Paup 
# had to convert from .nex to .nwk with figtree
#msa = read.tree("../PRD alignnment_ final.nex")
msa = readAAStringSet("../converted/PRD alignnment_ final.paup.fa",format="fasta", use.names=T)
msa_names = names(msa)
write.table(msa_names, "sample_names.PRD alignnment_ final.paup.fa.txt")

# ---------------------------------------------------------------------------------------
# ggtree.getCols (AA color pallette)
#
# copied from getCols from ggtree https://github.com/GuangchuangYu/ggtree/blob/master/R/utilities.R
# 
# used to build legend for AA<->color map
## who got it from ChIPseeker
ggtree.getCols <- function (n) {
  col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
           "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
           "#ccebc5", "#ffed6f")
  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
            "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
            "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
            "#ffff99", "#b15928")
  colorRampPalette(col3)(n)
}
# --- create a colorMap from an MSA
msaColorMap = function(msa) {
  aln  = BStringSet(msa)
  #width_fun <- get_fun_from_pkg("Biostrings", "width")
  window <- c(1, width(aln)[1])
  slice <- seq(window[1], window[2], by = 1)
  seqs <- lapply(1:length(aln), function(i) {
    x <- toString(aln[i])
    seq <- substring(x, slice, slice)
    seq[seq == "?"] <- "-"
    seq[seq == "*"] <- "-"
    seq[seq == " "] <- "-"
    return(seq)
  })
  names(seqs) <- names(aln)
  #if (is.null(color)) {
  alphabet <- unlist(seqs) %>% unique
  alphabet <- alphabet[alphabet != "-"]
  color <- ggtree.getCols(length(alphabet))
  names(color) <- alphabet
  color <- c(color, `-` = NA)
  #}
  
  return(color)
}
# -- plot a key for a colorMap
barplotColorMap = function(colorMap, title="AA color map", prefix="") {
  # make graphic of the color  
  bp = rep(1, length(colorMap))
  names(bp) = names(colorMap)
  midpoints= barplot(bp, main=paste0(prefix,title)
                     , col=colorMap
                     , ylim=c(0,10), space=0
                     , axes=F, xaxt='n', yaxt='n', ann=F
  ) #  legend=names(bp),
  text(midpoints, 0.5, labels=names(bp))
}
# test
barplotColorMapLegend(msaColorMap(msa), prefix="test: ")

# ======= LOAD NJ TREE ==============================================
#
# njTree = read.tree("../converted/PRD tree_final.figtree.nwk")
# read.tree breaks when there are colons in the node name.

# -- wrapper around ggtree::read.tree to handle colons in node names
read.tree.with_colons = function(file) {

  # pre-load raw text
  treeText = scan(file=file, what="", sep="\n", quiet=TRUE, skip=0, comment.char="#")
  #cat(treeText, "\n\n")
  
  # convert ":" to "___" (triple under)
  treeTextFixed = gsub("'([^':]+):([^':]+):([^':]+)'", "\\1___\\2___\\3", treeText)
  treeTextFixed = gsub("'([^':]+):([^':]+)'", "\\1___\\2", treeTextFixed)
  treeTextFixed = gsub("'([^':]+)'", "\\1", treeTextFixed)
  #cat(treeTextFixed,"\n\n")
  # parse the coded tree with ggtree
  tree = ggtree::read.tree(text=treeTextFixed)
  
  # restore the original node names
  tree$tip.label =gsub("___",":",tree$tip.label)
  
  return(tree)
}
# --- actually read & parse the file
njTree = read.tree.with_colons("../converted/PRD tree_final.figtree.nwk")
nj_tip_names= njTree$tip.label

# save samplenames out
write.table(nj_tip_names, "sample_names.converted.PRD tree_final.figtree.nwk.txt")

write.table(data.frame(groups=sort(row.names(bg)),nj=sort(nj_tip_names),msa=sort(msa_names)), "sample_names.all.txt", row.names = T)

# --- check node names

cat("expected node names: ", length(nj_tip_names), "\n")
cat("matching node names: ", sum(nj_tip_names %in% msa_names), "\n")

# --- build ggtree graphics
ggNJ = ggtree(njTree)
ggNJlabeled=ggNJ+geom_tiplab(align=F,  size=1.3) # aes(color=grp)
ggNJlabeledAligned=ggNJ + geom_tiplab(align=T, linesize=0.1, size=1.3)

ggNJ
# == PDF --- MSA + tree ------------------------------------
# 8 pages
#   * 2 color schemes
#     * legend page
#     * unlabelled, labelled, labelled and right-aligned.

aaColorNA = msaColorMap(msa)
comment(aaColorNA) = c(gap="bar")
aaColorWhite = aaColorNA; aaColorWhite["-"] = "white"
comment(aaColorWhite) = c(gap="white")

pdf("PRD_tree_alignment_final.nex.ggtree.msaplot.pdf", width=10, height=7.5)

for( colMap in list(aaColorNA, aaColorWhite)) {
  figTitle = paste0("PRD_tree_alignment_final.nex.ggtree.msaplot (gap=",comment(colMap)[["gap"]],"):\n " )
  barplotColorMap(colMap,prefix=figTitle)
  print(msaplot(ggNJ,               BStringSet(msa), color = colMap) )
  # Warning: Ignoring unknown aesthetics: x, y
  print(msaplot(ggNJlabeled,        BStringSet(msa), color = colMap, offset=0.08) )
  print(msaplot(ggNJlabeledAligned, BStringSet(msa), color = colMap, offset=0.09) )
}
dev.off()

#save.image("figure.ggtree.msa.PRD_final.RData")
#load("figure.ggtree.msa.PRD_final.RData")

