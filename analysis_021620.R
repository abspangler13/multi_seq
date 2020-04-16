x <- unique(stMwSeqESC$cell_type)

chosen.cell.types<-grep("MEF",stMwSeq$cell_type)
chosen.cell.types <- c(chosen.cell.types,grep("Stem",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("B cell",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("T cell",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("Epithel",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("Skele",stMwSeq$newAnn))
chosen.cell.types <- c(chosen.cell.types,grep("Brain",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("Endoth",stMwSeq$cell_type))
chosen.cell.types <- c(chosen.cell.types,grep("Macrophage",stMwSeq$cell_type))

stMwSeq[grep("Macrophage",stMwSeq$cell_type),]

stMwSeq.small<-stMwSeq[chosen.cell.types,]
