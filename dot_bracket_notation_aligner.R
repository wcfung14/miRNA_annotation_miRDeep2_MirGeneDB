#' @title MicroRNAs hairpin dot-bracket notation strand aligner
#' @description Given a microRNA hairpin dot-bracket notation, separate it into top and bottom strand by its hairpin loop, make a miRNA structure by aligning the top and bottom strands. Optional: map miRNA sequence to dot-bracket notation.
# To view miRNA dot-bracket notation structure, open the csv, set font to monospaced typefaces (e.g. Consolas)
#' @param hairpin_db string; miRNA hairpin dot-bracket notation. Obtained after mapping miRNA dot-bracket notation to miRNA sequences ("dot_bracket_notation_seq_mapper.R")
#' @param align_from_loop logical; If TRUE, the alignment start from the loop (robust to different number of "( and ")" in top and bottom strands). If FALSE, alignment starts from first character of strands (alignment may be wrong if number of "(" in the top strand is not the same as ")" in the bottom strand)
#' @param mirna_precursor_seq string; **optional** miRNA precursor sequence (consensus.precursor.sequence)
#' @param mirna_mature_seq string; **optional** miRNA mature sequence from results of mirdeep2 (consensus.mature.sequence)
#' @param mirna_star_seq string; **optional** miRNA star sequence from results of mirdeep2 (consensus.starsequence)
#' @return dataframe of 8 columns of miRNA hairpin alignment information: ("top_strand_align", "strand_match", "bottom_strand_align", "hairpin_structure_db", "Remarks", "bottom_strand_align_seq", "hairpin_structure_seq")
#' @export

align_dot_bracket <- function(hairpin_db, align_from_loop = TRUE, mirna_precursor_seq = NULL, mirna_mature_seq = NULL, mirna_star_seq = NULL) {
  
  # check if number of loop found in hairpin dot-bracket notation is larger than 1 
  # this method may not work if number of loop is larger than 1 (i.e. both ends of the miRNA structure are loops)
  if (stringr::str_count(hairpin_db, "\\(\\.+\\)") > 1) {
    output <- data.frame("top_strand_align" = NA, "strand_match" = NA, 
                         "bottom_strand_align" = NA, "hairpin_structure_db" = NA,
                         "Remarks" = "More than one loop found in hairpin dot-bracket notation!", "top_strand_align_seq" = NA, 
                         "bottom_strand_align_seq" = NA, "hairpin_structure_seq" = NA, row.names=NULL)
    return(output)
  }

  # find the strand on top of miRNA hairpin structure by structure from start of hairpin dot-bracket notation to the first ")" encountered
  top_strand <- sub("(^.*\\(\\.+)\\).*$", "\\1", hairpin_db)
  top_strand_rev_inv <- gsub("\\(","\\)", stringi::stri_reverse(top_strand)) # reverse and change ")" to "(" in top strand (for dot-bracket alignment starting from the loop)

  # find the strand at the bottom of miRNA hairpin structure by structure from the first "("" to the end of hairpin dot-bracket notation
  bottom_strand <- sub("^.*\\(\\.+(\\).*$)", "\\1", hairpin_db) 
  bottom_strand_rev_inv <- gsub("\\)","\\(", stringi::stri_reverse(bottom_strand)) # reverse and change "(" to ")" in bottom strand (for dot-bracket alignment starting from first character of strand, NOT NEEDED)

  # define dot-bracket align function to align top and bottom strands
  align_top_bottom_db <- function(top, bottom) {
    i <- 1; j <- 1; d <- 1 # initialize index for top, bottom and combined strands
    top_align <- NULL; bottom_align <- NULL
    top_split <- unlist(stringr::str_split(top, ""))
    bottom_split <- unlist(stringr::str_split(bottom, ""))
    
    # dot-bracket notation alignment 
    while (!is.na(top_split[i]) | !is.na(bottom_split[j])) { # while either strand is not empty
      if (top_split[i] %in% bottom_split[j]) { # when top and bottom strands are identical (match or mismatch), +1 to all index
        # use operator "%in%" for checking of NA logic, cannot be check by "==")
        top_align[d] <- top_split[i]
        bottom_align[d] <- bottom_split[j]
        i <- i+1
        j <- j+1
        d <- d+1
      } else if ("." %in% top_split[i]) { # when top strand is a mismatch, insert gap "-" to bottom strand, +1 to top and combined index
        top_align[d] <- top_split[i]
        bottom_align[d] <- "-"
        i <- i+1
        j <- j
        d <- d+1
      } else if ("." %in% bottom_split[j]) { # when bottom strand is a mismatch, insert gap "-" to top strand, +1 to bottom and combined index
        top_align[d] <- "-"
        bottom_align[d] <- bottom_split[j]
        i <- i
        j <- j+1
        d <- d+1
      } else if (is.na(bottom_split[j]) & !is.na(top_split[i])) { # when bottom strand reaches its end but top strand hasn't finished, insert gap "-" to bottom strand, +1 to top and combined index
        top_align[d] <- top_split[i]
        bottom_align[d] <- "-"
        i <- i+1
        j <- j
        d <- d+1
      } else if (is.na(top_split[i]) & !is.na(bottom_split[j])) { # when top strand reaches its end but bottom strand hasn't finished, insert gap "-" to top strand, +1 to bottom and combined index
        top_align[d] <- "-"
        bottom_align[d] <- bottom_split[j]
        i <- i
        j <- j+1
        d <- d+1
      } else { # break the loop when both strands reach their end
        break
      }
    }
    
    top_align <- paste(top_align, collapse = "") # collapse top_align from list to string
    top_align <- sub("^(\\.+)(-+)", "\\2\\1", top_align) # rearrange location of the first gap and mismatch, so ".--" becomes "--.", making it more intuitive to read (needed for alignment starting from first character of strand only)
    strand_match <- gsub("\\(|\\)", "|", gsub("-|\\.", " ", paste(top_align, collapse = ""))) # make match strand (match: "|", mismatch: " ")
    bottom_align <- gsub("\\(", "\\)", paste(bottom_align, collapse = "")) # collapse bottom_align from list to string and change "(" to ")" 
    
    # reverse top_align, bottom_align and strand_match for readability (reading from left to right) if top_align contains ")" but not "(" (i.e. if top_strand_rev_inv was used and alignment starts from the loop)
    if (grepl("\\)", top_align)) { 
      top_align <- gsub("\\)","\\(", stringi::stri_reverse(top_align))
      strand_match <- stringi::stri_reverse(strand_match)
      bottom_align <- stringi::stri_reverse(bottom_align)
    }
    
    ## make miRNA structure with dot-bracket notation
    hairpin_structure_db <- paste(top_align, strand_match, bottom_align, sep="\n")
    # END for dot-bracket notation alignment

    hairpin_structure_df <- data.frame("top_strand_align" = top_align, "strand_match" = strand_match, 
                                       "bottom_strand_align" = bottom_align, "hairpin_structure_db" = hairpin_structure_db,
                                       "Remarks" = NA, "top_strand_align_seq" = NA, "bottom_strand_align_seq" = NA,
                                       "hairpin_structure_seq" = NA, row.names=NULL)
    return(hairpin_structure_df)
  }
 
  if (align_from_loop) {
    # align starting from the loop (robust to different number of "( and ")" )
    db_aligned <- align_top_bottom_db(top = top_strand_rev_inv, bottom = bottom_strand) 
  } else {
    # align starting from first character of strand (does not work if number of "(" in the top strand is not the same as ")" in the bottom strand)
    db_aligned <- align_top_bottom_db(top = top_strand, bottom = bottom_strand_rev_inv)
    # extra "(" in the bottom strand happens when the overhang is also a match but not part of the miRNA precursor predicted.
  }
  
  
  # Additional: map mirna_precursor_seq to hairpin_db (if mirna_precursor_seq, mirna_mature_seq and mirna_star_seq are supplied in the function)
  if (all(!is.null(mirna_precursor_seq), !is.null(mirna_mature_seq), !is.null(mirna_star_seq))) {

    ## convert mature and star seq in mirna_precursor_seq to UPPER case
    mirna_precursor_seq_UPPER <- sub(mirna_mature_seq, toupper(mirna_mature_seq), sub(mirna_star_seq, toupper(mirna_star_seq), mirna_precursor_seq))
    
    # map mirna_precursor_seq for top stand
    {
      j <- 1 # initialize index for mirna_precursor_seq_UPPER
      top_strand_align_seq = ""
      for (i in 1:nchar(db_aligned$top_strand_align)) {
        if (unlist(stringr::str_split(db_aligned$top_strand_align, pattern = ""))[i] == "-") { # if top strand is a gap, insert gap to top_strand_align_seq
          top_strand_align_seq[i] <- "-"
        } else { # if top strand is not a gap (i.e. is match or mismatch), insert mirna_precursor_seq_UPPER to top_strand_align_seq, +1 to mirna_precursor_seq_UPPER index
          top_strand_align_seq[i] <- unlist(stringr::str_split(mirna_precursor_seq_UPPER, ""))[j]
          j <- j + 1
        }
      }
      top_strand_align_seq <- paste(top_strand_align_seq, collapse = "") # collapse top_strand_align_seq from list to string
    }
    # map mirna_precursor_seq for bottom stand
    {
      j <- 1 # initialize index for mirna_precursor_seq_UPPER
      bottom_strand_align_seq = ""
      for (i in 1:nchar(db_aligned$bottom_strand_align)) {
        if (unlist(stringr::str_split(db_aligned$bottom_strand_align, pattern = ""))[i] == "-") { # if bottom strand is a gap, insert gap to bottom_strand_align_seq
          bottom_strand_align_seq[i] <- "-"
        } else { # if bottom strand is not a gap (i.e. is match or mismatch), insert the reverse of mirna_precursor_seq_UPPER to bottom_strand_align_seq, +1 to mirna_precursor_seq_UPPER index
          bottom_strand_align_seq[i] <- unlist(stringr::str_split(stringi::stri_reverse(mirna_precursor_seq_UPPER), ""))[j]
          j <- j + 1
        }
      }
      bottom_strand_align_seq<- paste(bottom_strand_align_seq, collapse = "") # collapse bottom_strand_align_seq from list to string
    }
    
    # make miRNA structure with nucleotide seq
    hairpin_structure_seq <- paste(top_strand_align_seq, db_aligned$strand_match, bottom_strand_align_seq, sep="\n")
    
    db_aligned[, "top_strand_align_seq"] = top_strand_align_seq
    db_aligned[, "bottom_strand_align_seq"] = bottom_strand_align_seq
    db_aligned[, "hairpin_structure_seq"] = hairpin_structure_seq
  }

  return(db_aligned)
} # function END  
#===================================================

# runs only when script is run by itself
if (sys.nframe() == 0) {
  # Scaffold_1043_1001
  align_dot_bracket(hairpin_db = "(.(((((((.(((((.((((.(((((......).))))))))))))).))))))).)))")
  align_dot_bracket(hairpin_db = "(.(((((((.(((((.((((.(((((......).))))))))))))).))))))).)))", 
                    align_from_loop = TRUE,
                    mirna_precursor_seq = "aggccugaaucuuugucucaaccuucuguaaugaaaggugagcaaaguuucaggugugu",
                    mirna_mature_seq = "gugagcaaaguuucaggugugu", mirna_star_seq = "aggccugaaucuuugucuca")
  
  # Scaffold_1043_1003
  align_dot_bracket(hairpin_db = "(((((((((..(((.((((.((((.((...)).)))))))))))..))))))))))..")
  align_dot_bracket(hairpin_db = "(((((((((..(((.((((.((((.((...)).)))))))))))..))))))))))..", 
                    align_from_loop = TRUE,
                    mirna_precursor_seq = "cgccugaaagcuugacucaaccuuugcauagucaaggugagcaaaguuucaggugugu",
                    mirna_mature_seq = "gugagcaaaguuucaggugugu", mirna_star_seq = "cgccugaaagcuugacucaaccu")
  
  # Scaffold_6695_20322
  align_dot_bracket(hairpin_db = "(.((.((((...((((.((.((((.((((.((.((((...)))).)).)))).)))).)).))))...)))).)).)))")
  align_dot_bracket(hairpin_db = "(.((.((((...((((.((.((((.((((.((.((((...)))).)).)))).)))).)).))))...)))).)).)))", 
                    align_from_loop = TRUE,
                    mirna_precursor_seq = "cuugaccgggagcuggccacgaggaccucuugaccgggagcuggccacgaggaccucuugaccgggagcuggccacgag",
                    mirna_mature_seq = "ugaccgggagcuggccacgag", mirna_star_seq = "cuugaccgggagcuggccacg")
}
