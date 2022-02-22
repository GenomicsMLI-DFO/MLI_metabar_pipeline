

ggpcrplate.modif <- function (metabarlist, legend_title = "well_values", FUN = function(metabarlist) {
  rowSums(metabarlist$reads)
}) 
{
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    function_values <- FUN(metabarlist)
    if (length(function_values) != nrow(metabarlist$pcrs)) {
      stop("provided information should have the length of pcrs")
    }
    if (!is.numeric(function_values)) {
      stop("provided information should be numeric")
    }
    cols_plate_design <- c("plate_no", "plate_col", "plate_row")
    if (!all(cols_plate_design %in% colnames(metabarlist$pcrs))) {
      stop("PCR plate design not properly provided: ", 
           paste(cols_plate_design[!cols_plate_design %in% 
                                     colnames(metabarlist$pcrs)], sep = ", "), 
           " missing !\n")
    }
    plate_design <- metabarlist$pcrs[, c("plate_no", "plate_col", 
                                         "plate_row", "control_type")]
    plate_design_levels <- c("extraction", "pcr", "sequencing", 
                             "positive")
    plate_design$control_type <- factor(plate_design$control_type, 
                                        levels = plate_design_levels)
    plate_design$well_values <- function_values
    plate_design$well_values[plate_design$well_values == 
                               0] <- NA
    ggplot(plate_design, aes(y = match(.data$plate_row, 
                                       LETTERS[1:8]), x = .data$plate_col)) + 
      geom_tile(aes(fill = .data$control_type), height = 1, width = 1, na.rm = TRUE) + 
      facet_wrap(~.data$plate_no, scales = "free") + theme_bw() + 
      scale_y_reverse(breaks = 1:8, labels = LETTERS[1:8]) + 
      scale_x_continuous(breaks = 1:12) + scale_fill_manual(values = c("brown", 
                                                                       "red", "pink", "cyan4"), na.translate = FALSE) + 
      geom_point(na.rm = TRUE, aes(size = .data$well_values)) + labs(x = NULL, y = NULL, 
                                      size = legend_title)
  }
}


# ggplate for contaminant

ggpcrplate.cont <- function(metabarlist, N = Inf){
  
  motus.common_contam <- metabarlist$motus %>% filter(not_a_all_conta == F | not_a_max_conta == F)%>% arrange(desc(count)) %>% head(N)
  reads.common_contam <- metabarlist$reads[,rownames(motus.common_contam)]
  
  
  metabarlist_contam <- metabarlist_generator(reads.common_contam, motus.common_contam, metabarlist$pcrs, metabarlist$samples)
  
  title.int = ifelse(N > length(rownames(motus.common_contam)), "#reads of the ALL contaminants", paste("#reads of the", N, "most abundant contaminants") )
  
  ggpcrplate.modif(metabarlist_contam, legend_title = title.int,
                   FUN = function(m){rowSums(m$reads)})
  
}

