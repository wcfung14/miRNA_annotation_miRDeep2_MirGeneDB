#' @title Dot-bracket notation to microRNAs sequences mapper
#' @description Given a microRNA sequence, mature sequence, star sequence, and dot-bracket notation, map dot-bracket notation to miRNA sequences.
#' @param mirna_seq string; miRNA sequence from mirdeep2 .pdf results. Obtained after scrapping miRNA dot-bracket notation to miRNA sequences ("dot_bracket_notation_scrapper.R") (Has the same length with mirna dot_bracket_notation, not mirdeep2_res_db_seq$consensus.precursor.sequence)
#' @param mirna_mature_seq string; miRNA mature sequence from mirdeep2 results (consensus.mature.sequence)
#' @param mirna_star_seq string; miRNA star sequence from mirdeep2 results (consensus.star.sequence)
#' @param dot_bracket_notation string; dot-bracket notation from mirdeep2 results. Obtained after scrapping miRNA dot-bracket notation to miRNA sequences ("dot_bracket_notation_scrapper.R")
#' @return dataframe of 14 columuns of miRNA hairpin structure information: ("skip_this", "mature_loc_start", "mature_loc_end", "star_loc_start", "star_loc_end", "loop_loc_start", "loop_loc_end", "mature_first", "hairpin_db", "mature_db", "star_db", "extension_db", "loop_db", "Remarks")
#' @export

map_dot_bracket_notation_seq <- function(mirna_seq, mirna_mature_seq, mirna_star_seq, dot_bracket_notation) {
 
  skip_this <- ""
  # Catching errors
  # 1. Mature and star sequences can be found in both strands (e.g. Tco_Scaffold_6695_20322)
  if (stringr::str_count(mirna_seq, mirna_mature_seq) > 1 | stringr::str_count(mirna_seq, mirna_star_seq) > 1) {
    
    # check if number of mature/star sequence found in miRNA sequence is greater than 1 
    # if the number is greater than 1, finding hairpin sequence and aligning will be difficult -> skip
    output <- data.frame("skip_this" = TRUE, "mature_loc_start" = NA, "mature_loc_end" = NA, "star_loc_start" = NA, 
                   "star_loc_end" = NA, "loop_loc_start" = NA, "loop_loc_end" = NA, "mature_first" = NA, 
                   "hairpin_db" = NA, "mature_db" = NA, "star_db" = NA, "extension_db" = NA, "loop_db"= NA, 
                   "Remarks" = "Mature and star sequences can be found in both strands!", row.names=NULL)
    return(output)
  }
  
  # find start and end locus of mature, star sequences in miRNA seq / dot-bracket notation from results of mirdeep2
  # beware that the locus is relative to the results of mirdeep2 (precursor sequence), NOT relative to "consensus.precursor.sequence" (hairpin sequence) in "miRDeep2_miRNAs.norm_count.5p3p-arm-dominate.info.xls" (excel result output of mirdeep2)
  mature_loc <- stringr::str_locate(mirna_seq, mirna_mature_seq)
  star_loc <- stringr::str_locate(mirna_seq, mirna_star_seq)
  
  # divide dot-bracket notation into four substrings: mature, star, extension (hairpin - mature - star) and loop dot-bracket strings
  mature_db <- stringr::str_sub(string = dot_bracket_notation, start = mature_loc[, "start"], end = mature_loc[, "end"])
  star_db <- stringr::str_sub(string = dot_bracket_notation, start = star_loc[, "start"], end = star_loc[, "end"])
  
  mature_first <- "" # Boolean for whether mature seq is positioned before star seq
  hairpin_db <- "" # for finding the true loop inside hairpin sequences, not a false loop outside hairpin
  extension_db <- ""
 
  # find the dot-bracket notation of hairpin and extension
  if (mature_loc[, "end"] < star_loc[, "start"]) {
    mature_first <- TRUE
    
    hairpin_db <- stringr::str_sub(string = dot_bracket_notation, mature_loc[, "start"], star_loc[, "end"])
    extension_db <- stringr::str_sub(string = dot_bracket_notation, mature_loc[, "end"]+1, star_loc[, "start"]-1)
  } else {
    mature_first <- FALSE
      
    hairpin_db <- stringr::str_sub(string = dot_bracket_notation, star_loc[, "start"], mature_loc[, "end"])
    extension_db <- stringr::str_sub(string = dot_bracket_notation, star_loc[, "end"]+1, mature_loc[, "start"]-1)
  }
  
  # Catching errors
  # 2. Number of loops in the precursor sequence may be larger than 1 (e.g. Tco_Scaffold_1043_1001), so we have to grep the loop inside miRNA hairpin
  # Still, number of loops in the miRNA hairpin would still be larger than 1 (e.g. Tco_Scaffold_6695_20322, Hho_SczTNLB_6657_30897) 
  # for poorly aligned sequence (e.g. Tco_Scaffold_6695_20322), number of loops shown is 0 but is actually >1 
  # In that case, I can't do anything but to let you check the pdf file...
  if (stringr::str_count(hairpin_db, "\\(\\.+\\)") > 1) {
    # check if number of loop found in hairpin dot-bracket notation is larger than 1 
    # if the number is greater than 1, finding hairpin sequence and aligning will be difficult -> skip
    output <- data.frame("skip_this" = TRUE, "mature_loc_start" = NA, "mature_loc_end" = NA, "star_loc_start" = NA, 
                   "star_loc_end" = NA, "loop_loc_start" = NA, "loop_loc_end" = NA, "mature_first" = NA, 
                   "hairpin_db" = NA, "mature_db" = NA, "star_db" = NA, "extension_db" = NA, "loop_db"= NA, 
                   "Remarks" = "More than one loop found in hairpin dot-bracket notation!", row.names=NULL)
    return(output)
  }
  
  loop_loc <- stringr::str_locate(hairpin_db, "\\(\\.+\\)") # return loop loc in hairpin sequence, grep by structure that looks like "(....)"
  loop_loc <- loop_loc + min(mature_loc[, "start"], star_loc[, "start"])-1 # return loop loc in precursor sequence
  
  loop_db <- stringr::str_sub(string = dot_bracket_notation, start = loop_loc[, "start"], end = loop_loc[, "end"])
  
  output <- data.frame("skip_this" = FALSE, "mature_loc_start" = mature_loc[, "start"], "mature_loc_end" = mature_loc[, "end"],
                 "star_loc_start" = star_loc[, "start"], "star_loc_end" = star_loc[, "end"], 
                 "loop_loc_start" = loop_loc[, "start"], "loop_loc_end" = loop_loc[, "end"], "mature_first" = mature_first, 
                 "hairpin_db" = hairpin_db, "mature_db" = mature_db, "star_db" = star_db, 
                 "extension_db" = extension_db, "loop_db"= loop_db, "Remarks" = "", row.names=NULL)

  return(output)
} # function END
#===================================================

# runs only when script is run by itself
if (sys.nframe() == 0) {
  # Scaffold_1043_1001
  map_dot_bracket_notation_seq(mirna_seq = "uaucuccaucauuuucaugcgaaacuuucuggcaggccugaaucuuugucucaaccuucuguaaugaaaggugagcaaaguuucagguguguuugggggacuucuuuu", 
                               mirna_mature_seq = "gugagcaaaguuucaggugugu",
                               mirna_star_seq = "aggccugaaucuuugucuca",
                               dot_bracket_notation = "..(((((.....((((....))))......((((.(((((((.(((((.((((.(((((......).))))))))))))).))))))).))))..)))))........")
  
  # Scaffold_1043_1003
  map_dot_bracket_notation_seq(mirna_seq = "uuguaggucugugcauuuauacagcucuacuuggcgccugaaagcuugacucaaccuuugcauagucaaggugagcaaaguuucaggugugucagugggguuaaaaug", 
                               mirna_mature_seq = "gugagcaaaguuucaggugugu",
                               mirna_star_seq = "cgccugaaagcuugacucaaccu",
                               dot_bracket_notation = ".(((((((......))))))).(((((((((..((((((((((..(((.((((.((((.((...)).)))))))))))..))))))))))...)))))))))......")
  
  # Scaffold_6695_20322
  map_dot_bracket_notation_seq(mirna_seq = "ccacgaggaccucuugaccgggagcuggccacgaggaccucuugaccgggagcuggccacgaggaccucuugaccgggagcuggccacgaggaccucunnnnnnnnnnnnn", 
                               mirna_mature_seq = "ugaccgggagcuggccacgag",
                               mirna_star_seq = "cuugaccgggagcuggccacg",
                               dot_bracket_notation = "....((((.((((.((.((((...((((.((.((((.((((.((.((((...)))).)).)))).)))).)).))))...)))).)).)))).))))..............")
  
  
}
