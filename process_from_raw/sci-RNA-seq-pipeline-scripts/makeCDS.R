suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(stringr)
    library(monocle3)
})

args = commandArgs(trailingOnly = T)

mat.path = args[1]
gene.annotation.path = args[2]
cell.annotation.path = args[3]
output.dir = args[4]

load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
    df = read.table(
        mat.path,
        col.names = c("gene.idx", "cell.idx", "count"),
        colClasses = c("integer", "integer", "integer"))
    
    gene.annotations = read.table(
        gene.annotation.path,
        col.names = c("id", "gene_short_name"),
        colClasses = c("character", "character"))
    
    cell.annotations = read.table(
        cell.annotation.path,
        col.names = c("Cell", "sample"),
        colClasses = c("character", "factor"))
    
    rownames(gene.annotations) = gene.annotations$id
    rownames(cell.annotations) = cell.annotations$Cell
    
    # add a dummy cell to ensure that all genes are included in the matrix
    # even if a gene isn't expressed in any cell
    df = rbind(df, data.frame(
        gene.idx = c(1, nrow(gene.annotations)),
        cell.idx = rep(nrow(cell.annotations)+1, 2),
        count = c(1, 1)))
    
    
    mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
    mat = mat[, 1:(ncol(mat)-1)]
        

        

    print(dim(mat))
    rownames(mat) = gene.annotations$id
    colnames(mat) = cell.annotations$Cell
    
    pd = new("AnnotatedDataFrame", data = cell.annotations)
    fd = new("AnnotatedDataFrame", data = gene.annotations)

    cds = newCellDataSet(mat, phenoData = pd, featureData = fd)
    pData(cds)$n.umi = Matrix::colSums(exprs(cds))

    return(cds)
}

cds = load.cds(args[1],args[2],args[3])

saveRDS(object = cds,
    file = paste(output.dir, "/", "cds.RDS", sep = ""))
