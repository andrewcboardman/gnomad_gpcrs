convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

mouse.genes <- read.table('',stringsAsFactors = F)$V1
human.genes <- convertMouseGeneList(mouse.genes)
write.csv(human.genes,'human_gene_homologues.txt')
human.gpcr.genelist <- read.table('human_gpcr_genes.txt',stringsAsFactors = F)$V1
human.gpcr.genelist <- sapply(human.gpcr.genelist,function(x) {strsplit(x,",")[1][1]})




