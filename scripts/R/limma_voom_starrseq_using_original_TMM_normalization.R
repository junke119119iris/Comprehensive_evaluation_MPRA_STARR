# ========================
# Setup and Argument Parsing
# ========================
suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
    library(argparse)
    library(DESeq2)
    library(data.table)
})

# Setup argument parser
parser <- ArgumentParser()

# Define input arguments
parser$add_argument("-i", "--input", nargs="+", help="Input UMI count")
parser$add_argument("-o", "--outDir", nargs="+", help="Output directory")
parser$add_argument("--num_rep_DNA", type="integer", default=3, help="Number of DNA replicates")
parser$add_argument("--num_rep_RNA", type="integer", default=3, help="Number of RNA replicates")
parser$add_argument("--corresponding", required=F, action="store_true", default=FALSE, help="If DNA/RNA are corresponding pairs")
parser$add_argument("-s", required=F, action="store_true", default=FALSE, help="Separate orientation files")
parser$add_argument("--prefiltered", required=F, action="store_true", default=FALSE, help="Set if input count is already filtered")
parser$add_argument("--default_filter", required=F, action="store_true", default=FALSE)
parser$add_argument("--raw_filter", required=F, action="store_true", default=FALSE)
parser$add_argument("--cpm_filter", required=F, action="store_true", default=FALSE)
parser$add_argument("--count_cutoff", type="integer", default=10, help="Raw count cutoff value")
parser$add_argument("--num_lib_to_filter", type="integer", default=0, help="Min number of libraries to pass filter")
parser$add_argument("--DNA_filter", action="store_true", default=FALSE)
parser$add_argument("--average_lib_size", action="store_true", default=FALSE)
# parser$add_argument("--modified_TMM", required=F, action="store_true", default=FALSE, help="Use modified TMM normalization")

# Parse arguments
args <- parser$parse_args()

# ========================
# Load and prepare count matrix
# ========================
print(args$input)
count_mat = as.data.frame(fread(args$input, sep='\t', data.table=FALSE))
print(nrow(count_mat))

# Create DNA and RNA column labels
DNA_columns = paste0('DNA', c(1:args$num_rep_DNA))
RNA_columns = paste0('RNA', c(1:args$num_rep_RNA))

# Detect header format
if(ncol(count_mat) == (4+args$num_rep_DNA+args$num_rep_RNA)){
    col_name = c('seqnames', 'start', 'end', 'strand', DNA_columns, RNA_columns)
}else if(ncol(count_mat) == (4+1+args$num_rep_DNA+args$num_rep_RNA)){
    col_name = c('seqnames', 'start', 'end', 'name', 'strand', DNA_columns, RNA_columns)
}
colnames(count_mat) = col_name

# Create edgeR DGEList object
dge = DGEList(counts = count_mat[, c(DNA_columns, RNA_columns)], genes = count_mat[, c('seqnames', 'start', 'end', 'strand')])
print(dge$samples$lib.size)

# Add metadata
dge$samples$group = as.factor(c(rep('DNA', args$num_rep_DNA), rep('RNA', args$num_rep_RNA)))
dge$samples$replicate = as.factor(c(1:args$num_rep_DNA, 1:args$num_rep_RNA))

# ========================
# Filtering
# ========================
if(!args$prefiltered){
    print('Filter')
    if(args$cpm_filter){
        print('With CPM filter')
        cpm = cpm(dge, log=FALSE)

        # Choose how to define minimum library size
        if(args$DNA_filter){
            print('& Filter based on DNA libs')
            if(args$average_lib_size){
                print('& Use average lib size to calculate CPM threshold')
                lib_size = mean(dge$samples$lib.size[1:args$num_rep_DNA])
            }else{
                print('& Use minimal lib size to calculate CPM threshold')
                lib_size = min(dge$samples$lib.size[1:args$num_rep_DNA])
            }
        }else{
            print('& Filter based on all libs')
            if(args$average_lib_size){
                lib_size = mean(dge$samples$lib.size)
            }else{
                lib_size = min(dge$samples$lib.size)
            }
        }
        print('Calculate CPM threshold')
        cutoff = args$count_cutoff / (lib_size / 1e6)
        print('Cutoff')
        print(cutoff)

        if(args$num_lib_to_filter > 0){
            print('# of libs to filter')
            print(args$num_lib_to_filter)
            if(args$DNA_filter){
                if(args$num_rep_DNA == 1){
                    keep = cpm[, DNA_columns] >= cutoff
                }else{
                    keep = rowSums(cpm[, DNA_columns] >= cutoff) >= args$num_lib_to_filter
                }
            }else{
                keep = rowSums(cpm >= cutoff) >= args$num_lib_to_filter
            }
        }else{
            keep = rowSums(cpm >= cutoff) >= min(c(args$num_rep_DNA, args$num_rep_RNA))
        }
    }else{
        print("Filter with default edgeR's filterByExpr")
        # Default filtering method using edgeR’s filterByExpr
        keep = filterByExpr(dge, group = dge$samples$group)
    }

    print(table(keep))
    dge = dge[keep,, keep.lib.sizes=FALSE]
    # Ensure some RNA reads are retained
    check = rowSums(dge$counts[, RNA_columns] > 0) > 0
    print(table(check))
}

gc()  # Garbage collection

# Save filtered count matrix
if(!args$prefiltered){
    write.table(cbind(dge$genes, dge$counts), file=paste0(args$outDir, 'filtered_raw_count.txt'),
                sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    print('finished saving count matrix')
}

# ========================
# Normalization (TMM)
# ========================
dge = calcNormFactors(dge, method="TMM")
print(dge$samples$norm.factors)

# ========================
# Linear Modeling (limma-voom)
# ========================
# Create design matrix for DNA vs RNA
coldata = data.frame(
    RvD = c(rep(0, args$num_rep_DNA), rep(1, args$num_rep_RNA)),
    Rep = c(1:args$num_rep_DNA, 1:args$num_rep_RNA)
)
options(contrasts=c("contr.sum", "contr.poly"))
mdes = model.matrix(~RvD, data=coldata)

# If paired design (corresponding), apply duplicate correlation
if(args$corresponding){
    png(file=paste0(args$outDir, "mean_variance_trend_voom_plot.png"), width=5000, height=5000, res=300)
    vdata = voomWithQualityWeights(dge, design=mdes, plot=F)
    vcorrf = duplicateCorrelation(vdata, mdes, block=coldata$Rep)$consensus.correlation
    vdata = voomWithQualityWeights(dge, design=mdes, block=coldata$Rep, correlation=vcorrf, plot=T)
    vcorrf = duplicateCorrelation(vdata, mdes, block=coldata$Rep)$consensus.correlation
    dev.off()

    fit = eBayes(lmFit(vdata, mdes, block=coldata$Rep, correlation=vcorrf), robust=TRUE)
}else{
    png(file=paste0(args$outDir, "mean_variance_trend_voom_plot.png"), width=5000, height=5000, res=300)
    vdata = voomWithQualityWeights(dge, design=mdes, plot=T)
    dev.off()

    fit = eBayes(lmFit(vdata, mdes), robust=TRUE)
}

# Extract differential expression results
hits = topTable(fit, coef='RvD', number=Inf, sort.by="none")

# Write output
write.table(hits, file=paste0(args$outDir, 'limma_out.txt'), sep='\t', quote=FALSE)
print('finished')