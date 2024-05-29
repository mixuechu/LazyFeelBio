
library(VariantAnnotation)

source("r_packs/send_progress.R")

.collapseLists <- VariantAnnotation:::.collapseLists
.formatALT <- VariantAnnotation:::.formatALT
.formatInfo <- VariantAnnotation:::.formatInfo


.readVcf <- function(file, genome, param, row.names, progress_callback = NULL, ...)
{
    if (missing(genome))
        genome <- seqinfo(scanVcfHeader(file))
    if (!is(genome, "character") & !is(genome, "Seqinfo"))
        stop("'genome' must be a 'character(1)' or 'Seqinfo' object")
    if (is(genome, "Seqinfo")) {
        if (is(param, "ScanVcfParam"))
            chr <- names(vcfWhich(param))
        else
            chr <- seqlevels(param)
        ## confirm param seqlevels are in supplied Seqinfo
        if (any(!chr %in% seqnames(genome)))
            stop("'seqnames' in 'vcfWhich(param)' must be present in 'genome'")
    }
    # 更新进度：开始扫描VCF文件
    if (!is.null(progress_callback)) progress_callback("Scanning VCF file")
    vcf_data <- scanVcf(file, param=param, row.names=row.names, ...)
    # 更新进度：扫描完成，开始转换VCF
    if (!is.null(progress_callback)) progress_callback("Converting VCF data")
    .scanVcfToVCF(vcf_data, file, genome, param, progress_callback)
}

.scanVcfToVCF <- function(vcf, file, genome, param, progress_callback = NULL, ...)
{
    hdr <- scanVcfHeader(file)
    if (length(vcf[[1]]$GENO) > 0L)
        colnms <- colnames(vcf[[1]]$GENO[[1]])
    else
        colnms <- NULL

    # 更新进度：开始合并列表
    if (!is.null(progress_callback)) progress_callback("Collapsing lists")
    vcf <- .collapseLists(vcf, param)

    ## rowRanges
    rowRanges <- vcf$rowRanges
    if (length(rowRanges)) {
        if (is(genome, "character")) {
           if (length(seqinfo(hdr))) {
               merged <- merge(seqinfo(hdr), seqinfo(rowRanges))
               map <- match(names(merged), names(seqinfo(rowRanges)))
               seqinfo(rowRanges, map) <- merged
           }
           genome(rowRanges) <- genome
        } else if (is(genome, "Seqinfo")) {
            if (length(seqinfo(hdr)))
                reference <- merge(seqinfo(hdr), genome)
            else
                reference <- genome
            merged <- merge(reference, seqinfo(rowRanges))
            map <- match(names(merged), names(seqinfo(rowRanges)))
            seqinfo(rowRanges, map) <- merged
        }
    }
    values(rowRanges) <- DataFrame(vcf["paramRangeID"])

    ## fixed fields
    fx <- vcf[c("REF", "ALT", "QUAL", "FILTER")]
    fx$ALT <- .formatALT(fx$ALT)
    fixed <- DataFrame(fx[!sapply(fx, is.null)])

    # 更新进度：格式化 INFO 字段
    if (!is.null(progress_callback)) progress_callback("Formatting INFO fields")
    ## info
    info <- .formatInfo(vcf$INFO, info(hdr), length(rowRanges))

    ## colData
    colData <- DataFrame(Samples=seq_along(colnms), row.names=colnms)

    ## geno
    geno <- SimpleList(lapply(vcf$GENO, `dimnames<-`, NULL))

    vcf <- NULL
    VCF(rowRanges=rowRanges, colData=colData, exptData=list(header=hdr),
        fixed=fixed, info=info, geno=geno)
}

readVcfGz <- function(inputFile){
 data <- .readVcf(inputFile, param = ScanVcfParam(), row.names = TRUE, progress_callback = progress_callback)
 data
}