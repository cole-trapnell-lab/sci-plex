suppressPackageStartupMessages({
    library(methods)
    library(dplyr)
    library(ggplot2)
})

args = commandArgs(trailingOnly = T)

if (length(args) < 2) {
    cat("Usage: Rscript knee-plot.R UMIs.per.cell.barcode output-dir/ [UMI-per-cell cutoff]", file = stderr())
    exit(1)
}

input.file = args[1]
output.dir = args[2]
if (length(args) >= 3)
    cutoff = as.integer(args[3])

df = read.table(
    args[1],
    col.names = c("sample", "barcode", "n.umi"),
    colClasses = c("factor", "character", "integer"))

df = df %>% group_by(1) %>% mutate(n.umi.rank = min_rank(-n.umi)) %>% ungroup()

plot = ggplot(
    df %>% arrange(-n.umi) %>% select(n.umi, n.umi.rank) %>% distinct(),
    aes(x = n.umi.rank, y = n.umi)) +
    geom_line(size = 0.8) +
    scale_x_log10(limits = c(10, NA),
                  breaks = c(10, 20, 40, 100, 200, 400, 1000, 2000, 4000, 8000, 16000)) +
    scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000)) +
    xlab("# of barcodes") +
    ylab("UMI count threshold") +
    theme_bw()

if (length(args) >= 3) {
    plot = plot +
        geom_hline(yintercept = cutoff, size = 1.2, color = "firebrick2")
}

ggsave(paste(args[2], "/all.samples.pdf", sep = ""),
    plot = plot, units = "in", width = 3.5*1.618, height = 3.5)

