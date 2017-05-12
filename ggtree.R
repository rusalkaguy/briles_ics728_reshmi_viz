#
# ggtree
#
# https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/treeImport.html

library("ggplot2")
# library("ggtree") not yet availble
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
library(ggtree)
# biocLite("phylobase")
library(phylobase)

#
# load annotation
#

brilesGroups = read.delim("../../v2/CODING_MOTIFS_ICS728.161121.groups.txt", sep="\t", header=TRUE, stringsAsFactors=F)
#bg = brilesGroups[,c(-1,-2)]
#row.names(bg) = brilesGroups$name
bg = data.frame(row.names=brilesGroups$name, grp=brilesGroups$ColorCode)


# load Reshmi's NJ tree
#msa = read.tree("../PRD alignnment_ final.nex")
msa = readAAStringSet("../converted/PRD alignnment_ final.paup.fa",format="fasta", use.names=T)
njTree = read.tree("../converted/PRD tree_final.figtree.nwk")
ggNJ = ggtree(njTree)
ggNJ
# ---------------------------------------
msaplot(ggNJ, BStringSet(msa)) #### ERROR - nj_tip_names have single-quotes sometimes!!!
#----------------------------------------
msa_names = names(msa)
nj_tip_names= njTree$tip.label
sum(nj_tip_names %in% msa_names)
#[1] 106
length(nj_tip_names)
#[1] 136
length(msa_names)
#[1] 136
nj_tip_names[!(nj_tip_names %in% msa_names)]
# [1] "'A66.1_PRD'"    "'126.127'"      "'160.161"       "'170.171"       "'129.130"       "'147.148"       "'155.156"      
# [8] "'175.176"       "'101.102.144"   "'135.136'"      "'141.142"       "'084.95'"       "'111.113.207'"  "'116.117.208'" 
# [15] "'112.114'"      "'104.105"       "'051.52"        "'048.86'"       "'022.23"        "'026.27'"       "'039.96'"      
# [22] "'074.75'"       "'001.2"         "'006.7"         "'056.58.77.83'" "'016.61.78'"    "'018.19"        "'029.30"       
# [29] "'033.93.119'"   "'069.70"       
msa_names[!(msa_names %in% nj_tip_names)]
# [1] "A66.1_PRD"            "001.2:5.73"           "006.7:15.59:60.76.92" "016.61.78"            "018.19:21.62.82"     
# [6] "022.23:25.44.66"      "026.27"               "029.30:32"            "033.93.119"           "039.96"              
# [11] "048.86"               "051.52:55.103"        "056.58.77.83"         "069.70:72"            "074.75"              
# [16] "084.95"               "101.102.144:145"      "104.105:109"          "111.113.207"          "112.114"             
# [21] "116.117.208"          "126.127"              "129.130:134"          "135.136"              "141.142:143"         
# [26] "147.148:154"          "155.156:159"          "160.161:162"          "170.171:174"          "175.176:180"         

# ================================================================= 
# iterate over K
# k=33 #DEBUG
for( k in seq(9,45,by=3) ) {
# load file
  file<-paste0("../../v2/aaf/all_uncoded.161121.k",k,".tre")
  tree = read.tree(file)
  comment(tree) =c(filename=file, k=k)

# graph
# ggtree(tree) + geom_tiplab() + geom_point(aes(color=isTip), size=5, alpha=.5)
  #+ geom_text(aes(label=branch.length, x=branch), vjust=-.5)


#
# annotate tree
#
# annoTree = tree %<+% brilesGroups
# Error in fix.by(by.x, x) : 'by' must specify a uniquely valid column

# need rownamesAsLabels in order not to interpret numeric rownames as rowNUMBERS!
dbat = phylo4d(tree,bg, rownamesAsLabels=T)

# display 
#ggtree(dbat) + geom_tiplab() + geom_point(aes(color=grp),size=5, alpha=0.5)


pdffile<-paste0("../../v2/aaf/all_uncoded.161121.k",k,".tre.ggtree.pdf")
pdf(pdffile,width=8,height=10)
# 
# with Group
#
p =ggtree(dbat)
pp = (p +geom_treescale(x=.02) + geom_tiplab(align=T, aes(color=grp), size=2))
ppp=gheatmap(pp, bg, offset=0.02, width =0.2)
print(ppp)

# display with AA
seqs = BStringSet(brilesGroups$seq, start=1)
names(seqs)=brilesGroups$name
print(msaplot(pp, seqs, offset=0.03))

dev.off()

}
#
# convert TREE to DATA.FRAME
#
#tree_data = fortify(tree)
#head(tree_data)

