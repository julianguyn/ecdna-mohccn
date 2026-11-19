# helper functions for 1-CohortAnalysis/3-ecDNA_genomic_regions.R

#' parse AA Location strings ("['chr1:100-200', 'chr2:300-400']") into one
#' GRanges with one row per segment, tagged with sample, feature and filter flag
parse_ecDNA_segments <- function(df) {
    segs <- lapply(seq_len(nrow(df)), function(i) {
        loc <- regmatches(df$Location[i], gregexpr("chr[^:']+:[0-9]+-[0-9]+", df$Location[i]))[[1]]
        if (length(loc) == 0) return(NULL)
        p <- do.call(rbind, strsplit(loc, "[:-]"))
        data.frame(
            chr = p[, 1],
            start = as.numeric(p[, 2]),
            end = as.numeric(p[, 3]),
            sample = df$`Sample name`[i],
            feature = df$`Feature ID`[i],
            filter_flag = df$`Filter flag`[i]
        )
    })
    segs <- do.call(rbind, segs)
    makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
}

#' N-gap mask for hg38 (assembly gaps; also covers acrocentric p-arms)
get_N_mask <- function(chroms) {
    genome <- BSgenome.Hsapiens.UCSC.hg38
    mask <- lapply(chroms, function(chr) {
        ir <- as(masks(maskMotif(genome[[chr]], "N")), "NormalIRanges")
        data.frame(chr = rep(chr, length(ir)), start = start(ir), end = end(ir))
    })
    makeGRangesFromDataFrame(do.call(rbind, mask))
}

#' hg38 centromeres (cytoBand "acen" bands) from UCSC
#' Centromeres modelled alpha-satellite sequence in hg38 (not N),
#' attract spurious AA calls
get_centromere_mask <- function(chroms,
                                url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz") {
    tmp <- tempfile(fileext = ".txt.gz")
    download.file(url, tmp, quiet = TRUE)
    cyto <- read.table(tmp, sep = "\t",
                       col.names = c("chr", "start", "end", "band", "stain"))
    cen <- cyto[cyto$stain == "acen" & cyto$chr %in% chroms, ]
    cen$start <- cen$start + 1  # UCSC is 0-based
    reduce(makeGRangesFromDataFrame(cen))
}

#' Randomly place ecDNA segments in callable genome
#' Segments from the same feature on the same chromosome are moved together as
#' a block (keeping their spacing), so the null preserves ecDNA size and structure.
shuffle_blocks <- function(segs, block_id, block_start, block_width, chrom_len, mask, max_tries = 100) {
    n_blocks <- length(block_start)
    cum_len <- c(0, cumsum(as.numeric(chrom_len)))
    total_len <- cum_len[length(cum_len)]

    new_chr <- character(n_blocks)
    new_start <- numeric(n_blocks)
    todo <- seq_len(n_blocks)

    for (try in seq_len(max_tries)) {
        if (length(todo) == 0) break
        # uniform position in concatenated genome -> chromosome + position
        pos <- runif(length(todo), 0, total_len)
        chr_idx <- findInterval(pos, cum_len, rightmost.closed = TRUE)
        start <- floor(pos - cum_len[chr_idx]) + 1
        new_chr[todo] <- names(chrom_len)[chr_idx]
        new_start[todo] <- start

        # reject blocks running off the chromosome end
        off_end <- start + block_width[todo] - 1 > chrom_len[chr_idx]

        # reject blocks with any segment centred in a masked region
        # (midpoint rather than any overlap, so large segments can still span small gaps)
        shift <- new_start[block_id] - block_start[block_id]
        mid <- floor((start(segs) + end(segs)) / 2) + shift
        hits_mask <- overlapsAny(GRanges(new_chr[block_id], IRanges(mid, width = 1)), mask)
        bad_blocks <- unique(block_id[hits_mask])

        todo <- todo[off_end | todo %in% bad_blocks]
    }

    # anything still failing after max_tries is kept as-is (rare: very large blocks)
    shift <- new_start[block_id] - block_start[block_id]
    GRanges(new_chr[block_id], IRanges(start(segs) + shift, end(segs) + shift),
            sample = segs$sample)
}

#' drop LowCN ecDNA and/or short segments
filter_segments <- function(gr, drop_lowCN = TRUE, drop_small_segs = TRUE, min_seg_width = 1e3) {
    if (drop_lowCN) gr <- gr[gr$filter_flag != "LowCN"]
    if (drop_small_segs) gr <- gr[width(gr) >= min_seg_width]
    gr
}

#' count, per bin, the number of distinct samples with an ecDNA segment in it
count_samples_per_bin <- function(bins, gr) {
    hits <- findOverlaps(bins, gr)
    pairs <- unique(data.frame(bin = queryHits(hits), sample = gr$sample[subjectHits(hits)]))
    tabulate(pairs$bin, nbins = length(bins))
}

#' Bin enrichment function
#' 
bin_enrichment <- function(ecDNA_segments, label,
                           bin_size = 1e6, n_perm = 1000,
                           drop_lowCN = TRUE,
                           drop_small_segs = TRUE, min_seg_width = 1e3,
                           chroms = c(paste0("chr", 1:22), "chrX"),
                           mask_centromeres = TRUE,
                           max_masked_frac = 0.5,
                           p_method = c("poisson", "empirical"),
                           seed = 123) {

    set.seed(seed)
    gr <- ecDNA_segments

    # ---- filtering
    gr <- filter_segments(gr, drop_lowCN, drop_small_segs, min_seg_width)
    gr <- keepSeqlevels(gr, intersect(chroms, seqlevels(gr)), pruning.mode = "coarse")
    seqlevels(gr) <- chroms

    chrom_len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chroms]
    seqlengths(gr) <- chrom_len
    message(label, ": ", length(gr), " segments from ", length(unique(gr$feature)),
            " ecDNA in ", length(unique(gr$sample)), " samples")

    # ---- bins, dropping those that are mostly masked (N gaps +/- centromeres)
    mask <- get_N_mask(chroms)
    if (mask_centromeres) mask <- reduce(c(mask, get_centromere_mask(chroms)))
    bins <- tileGenome(chrom_len, tilewidth = bin_size, cut.last.tile.in.chrom = TRUE)
    ov <- findOverlaps(bins, mask)
    masked_bp <- tapply(width(pintersect(bins[queryHits(ov)], mask[subjectHits(ov)])),
                        queryHits(ov), sum)
    masked_frac <- numeric(length(bins))
    masked_frac[as.integer(names(masked_bp))] <- masked_bp / width(bins)[as.integer(names(masked_bp))]
    bins <- bins[masked_frac <= max_masked_frac]

    # ---- observed: number of samples with ecDNA overlapping each bin
    obs_counts <- count_samples_per_bin(bins, gr)

    # ---- permutation null: shuffle each ecDNA block across the callable genome
    block_key <- paste(gr$feature, seqnames(gr), sep = "|")
    block_id <- match(block_key, unique(block_key))
    block_start <- tapply(start(gr), block_id, min)
    block_end <- tapply(end(gr), block_id, max)
    block_width <- block_end - block_start + 1

    n_ge <- numeric(length(bins))    # perms with null count >= observed
    null_sum <- numeric(length(bins))  # for expected count per bin

    for (i in seq_len(n_perm)) {
        rand_gr <- shuffle_blocks(gr, block_id, block_start, block_width, chrom_len, mask)
        null_counts <- count_samples_per_bin(bins, rand_gr)
        n_ge <- n_ge + (null_counts >= obs_counts)
        null_sum <- null_sum + null_counts
        if (i %% 100 == 0) message("Perm ", i, " of ", n_perm)
    }

    # empirical p = (perms with count >= observed + 1) / (n_perm + 1)
    # floored at 1 / (n_perm + 1), so kept for reference only
    pval_perm <- (n_ge + 1) / (n_perm + 1)

    # expected samples per bin under the null; pseudocount avoids expected = 0
    # in bins no permutation ever hit
    expected <- (null_sum + 1) / (n_perm + 1)

    # Poisson upper-tail p: P(X >= observed | lambda = expected).
    # Per-bin count is a sum of rare per-sample events, so Poisson is a good
    # approximation, and it isn't limited by the number of permutations.
    pval_pois <- ppois(obs_counts - 1, lambda = expected, lower.tail = FALSE)

    p_method <- match.arg(p_method)
    pvals <- if (p_method == "poisson") pval_pois else pval_perm
    padj <- p.adjust(pvals, method = "BH")

    # ---- compile results
    results <- data.frame(
        chr = as.character(seqnames(bins)),
        start = start(bins),
        end = end(bins),
        count = obs_counts,          # samples with ecDNA in bin
        expected = expected,         # mean samples per bin under null
        fold_enrichment = obs_counts / expected,
        pval_perm = pval_perm,
        pval_pois = pval_pois,
        pval = pvals,                # p used for padj (per p_method)
        padj = padj,
        label = label
    )
    results$enriched <- ifelse(results$padj < 0.05 & results$count > 0, "Enriched", "Not Enriched")
    results <- results[order(results$padj, -results$count), ]
    return(results)
}


#' get bin frequency
#' % of samples with ecDNA in each bin; significance carried over from bin_enrichment
bin_frequency <- function(enrich_res, n_samples, sig_threshold = 0.05,
                          chroms = c(paste0("chr", 1:22), "chrX")) {
    freq <- enrich_res[enrich_res$chr %in% chroms, ]
    freq$n_samples <- n_samples
    freq$freq <- 100 * freq$count / n_samples
    freq$hotspot <- freq$padj < sig_threshold & freq$count > 0
    freq[order(match(freq$chr, chroms), freq$start), ]
}

#' get hg38 gene coordinates with HGNC symbols to label hotspots
get_gene_coords <- function() {
    g <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
    g$symbol <- suppressMessages(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = g$gene_id, column = "SYMBOL", keytype = "ENTREZID")
    )
    g[!is.na(g$symbol)]
}

#' get all hotpost regions
#' merge adjacent hotspot bins into regions
#' annotate each with no. samples with ecDNA in the region, peak bin frequency, and oncogenes (on ecDNA) that physically lie in the region,
#' ranked by number of samples
get_hotspot_regions <- function(freq, ecDNA_segments, ecDNA_df, n_samples,
                                drop_lowCN = TRUE,
                                drop_small_segs = TRUE, min_seg_width = 1e3,
                                n_genes = 2) {
    hs <- freq[freq$hotspot, ]
    if (nrow(hs) == 0) return(NULL)
    hs_gr <- makeGRangesFromDataFrame(hs, keep.extra.columns = TRUE)
    regions <- reduce(hs_gr)  # abutting 1 Mb bins are merged

    # same segment filters as bin_enrichment
    gr <- filter_segments(ecDNA_segments, drop_lowCN, drop_small_segs, min_seg_width)

    genes_gr <- get_gene_coords()
    feat_onco <- setNames(
        strsplit(gsub("\\[|\\]|'", "", ecDNA_df$Oncogenes), ",\\s*"),
        ecDNA_df$`Feature ID`
    )

    out <- lapply(seq_along(regions), function(i) {
        reg <- regions[i]
        bins_in <- hs_gr[overlapsAny(hs_gr, reg)]
        peak <- bins_in[which.max(bins_in$count)]

        # ecDNA (samples / features) overlapping the region
        segs_in <- gr[overlapsAny(gr, reg)]
        samples_in <- unique(segs_in$sample)

        # oncogenes carried by those ecDNA that lie within the region
        feats <- unique(segs_in$feature)
        onco <- unique(data.frame(
            sample = rep(segs_in$sample[match(feats, segs_in$feature)], lengths(feat_onco[feats])),
            gene = unlist(feat_onco[feats], use.names = FALSE)
        ))
        onco <- onco[onco$gene != "" & onco$gene %in% genes_gr$symbol[overlapsAny(genes_gr, reg)], ]
        gene_counts <- sort(table(onco$gene), decreasing = TRUE)
        top_genes <- names(gene_counts)[seq_len(min(n_genes, length(gene_counts)))]

        data.frame(
            chr = as.character(seqnames(reg)),
            start = start(reg),
            end = end(reg),
            n_bins = length(bins_in),
            n_samples_ecDNA = length(samples_in),
            freq = 100 * length(samples_in) / n_samples,
            peak_start = start(peak),
            peak_end = end(peak),
            peak_count = peak$count,
            peak_freq = 100 * peak$count / n_samples,
            min_padj = min(bins_in$padj),
            oncogenes = paste(names(gene_counts), collapse = ","),
            label = if (length(top_genes) > 0) paste(top_genes, collapse = "/")
                    else sprintf("%s:%g-%gMb", seqnames(reg), round((start(reg) - 1) / 1e6),
                                 round(end(reg) / 1e6))
        )
    })
    out <- do.call(rbind, out)
    out[order(-out$n_samples_ecDNA), ]
}

#' plot mahattan plot
#' frequency bin but colour by enrichment
plot_frequency <- function(freq, regions = NULL, sig_threshold = 0.05, n_labels = 15,
                           chroms = c(paste0("chr", 1:22), "chrX")) {

    # cumulative genomic positions for x axis (full chromosome lengths)
    chrom_len <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[chroms]
    offset <- setNames(c(0, cumsum(as.numeric(chrom_len)))[seq_along(chroms)], chroms)

    toPlot <- freq[freq$chr %in% chroms, ]
    toPlot$xmin <- toPlot$start + offset[toPlot$chr]
    toPlot$xmax <- toPlot$end + offset[toPlot$chr]
    toPlot$fill <- ifelse(toPlot$hotspot, "Hotspot",
                          ifelse(match(toPlot$chr, chroms) %% 2 == 0, "even", "odd"))

    axis_set <- data.frame(label = sub("chr", "", chroms), center = offset + chrom_len / 2)

    p <- ggplot(toPlot) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = freq, fill = fill)) +
      scale_fill_manual(
        values = c(odd = "#898989", even = "#AEAEAE", Hotspot = "#B13434"),
        breaks = "Hotspot",
        labels = paste0("Enriched (FDR < ", sig_threshold, ")"),
        name = NULL
      ) +
      scale_x_continuous(breaks = axis_set$center, labels = axis_set$label,
                         expand = c(0.005, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Chromosome", y = "Samples with ecDNA (%)") +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 9),
        legend.position = c(0.99, 0.98),
        legend.justification = c(1, 1),
        legend.key.size = unit(3, "mm")
      )

    # label the top hotspot regions at their peak bin
    if (!is.null(regions) && nrow(regions) > 0) {
        lab <- head(regions[order(-regions$peak_freq), ], n_labels)
        lab$x <- (lab$peak_start + lab$peak_end) / 2 + offset[lab$chr]
        p <- p + ggrepel::geom_text_repel(
          data = lab, aes(x = x, y = peak_freq, label = label),
          size = 2.5, fontface = "italic", nudge_y = 0.3,
          segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf
        )
    }
    return(p)
}