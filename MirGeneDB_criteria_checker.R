#' @title MirGeneDB microRNAs criteria checker for miRNAs annotation
#' @description Check if the potential microRNAs (predicted by mirDeep2) fulfill the criteria 
#' of microRNAs of MirGeneDB (MirGeneDB 2.0, https://mirgenedb.org/information) 
#' by excel result output of mirdeep2 and miRNA's dot-bracket notation.
#' IMPORTANT: Rule 2 is NOT checked currently, and mirgenedb_rule_all_df does not check for rule 2
#' @param mirdeep2_res_db_seq_structure data.frame; excel result output of mirdeep2, after mapping microRNAs dot-bracket notation to miRNA sequences and return miRNA structure information (mature, star and loop sequences) ("dot_bracket_notation_seq_mapper.R")
#' @return dataframe of 9 columns ("mirna_id", "mirgenedb_rule_1", "mirgenedb_rule_2", mirgenedb_rule_3_AtLeast16bp", "mirgenedb_rule_3_ImperfectComplementarity", "mirgenedb_rule_4", "mirgenedb_rule_5", "mirgenedb_rule_6", "mirgenedb_rule_all")
#' @export

check_MirGeneDB_criteria <- function(mirdeep2_res_db_seq_structure) {
  mirna_id <- mirdeep2_res_db_seq_structure$provisional.id
  
  # check for MirGeneDB rule 1: both arms (mature and star) are expressed
  mirgenedb_rule_1_bothArmsExpressed <- mirdeep2_res_db_seq_structure$mature.read.count > 1 & mirdeep2_res_db_seq_structure$star.read.count > 1
  summary(mirgenedb_rule_1_bothArmsExpressed)
  
  
  # check for MirGeneDB rule 2: 5' ends of the reads are homogeneous
  # I have no idea what that means, so skip this first LOL
  mirgenedb_rule_2_5end_homogeneous <- NA
  
  
  # check for MirGeneDB rule 3: require imperfect complementarity, and base pairs in at least 16 of the ~22nt
  mirgenedb_rule_3_imperfectComplementarity <- stringr::str_count(mirdeep2_res_db_seq_structure$mature_db, pattern = "\\.") > 0 # check for imperfect complementarity by counting "." in mature sequence dot-bracket notation
  mirgenedb_rule_3_bpNumber <- stringr::str_count(mirdeep2_res_db_seq_structure$mature_db, pattern = "\\(|\\)") # check for number of base pairs by counting "(" or ")" in mature sequence dot-bracket notation
  mirgenedb_rule_3_atLeast16bp <- mirgenedb_rule_3_bpNumber >= 16 # check for base pairs in at least 16 nt
  
  mirgenedb_rule_3_imperfectComplementarity_atLeast16bp <- mirgenedb_rule_3_imperfectComplementarity & mirgenedb_rule_3_atLeast16bp
  
  summary(mirgenedb_rule_3_imperfectComplementarity)
  summary(mirgenedb_rule_3_bpNumber)
  summary(mirgenedb_rule_3_atLeast16bp)
  summary(mirgenedb_rule_3_imperfectComplementarity_atLeast16bp)
  
  
  # check for MirGeneDB rule 4: the 3p and 5p ends are offset by 2nt, or 3' end of the 3' arm is monouridilyated ("u"
  mirgenedb_rule_4_5end_offsetNumber <- ""
  mirgenedb_rule_4_3end_offsetNumber <- ""
  mirgenedb_rule_4_3end_nt <- ""
  
  # check for whether the 3p and 5p ends are offset by 2nt by checking position of max/min  mature and star seq (UPPER case) position
  for (i in 1:nrow(mirdeep2_res_db_seq_structure)) {
    # 5' end = minimum position of mature/star seq in bottom_strand_align_seq - minimum position of mature/star seq in top_strand_align_seq
    mirgenedb_rule_4_5end_offsetNumber[i] <- 
      min(unlist(gregexpr(pattern = "[A-Z]", mirdeep2_res_db_seq_structure$top_strand_align_seq[i]))) - min(unlist(gregexpr(pattern = "[A-Z]", mirdeep2_res_db_seq_structure$bottom_strand_align_seq[i])))
      
    # 3' end = maximum position of mature/star seq in bottom_strand_align_seq - maximum position of mature/star seq in top_strand_align_seq
    mirgenedb_rule_4_3end_offsetNumber[i] <-
      max(unlist(gregexpr(pattern = "[A-Z]", mirdeep2_res_db_seq_structure$top_strand_align_seq[i]))) - max(unlist(gregexpr(pattern = "[A-Z]", mirdeep2_res_db_seq_structure$bottom_strand_align_seq[i])))
    
    # check for whether the 3' end of the 3'arm is monouridilyated ("U")
    mirgenedb_rule_4_3end_nt[i] <-
      substr(mirdeep2_res_db_seq_structure$bottom_strand_align_seq[i], start = 1, stop = 1)
  }
  mirgenedb_rule_4_2nt_5end_2ntOffset <- abs(as.numeric(mirgenedb_rule_4_5end_offsetNumber)) >= 2
  mirgenedb_rule_4_2nt_3end_2ntOffset <- abs(as.numeric(mirgenedb_rule_4_3end_offsetNumber)) >= 2
  mirgenedb_rule_4_2nt_3end_isU <- tolower(mirgenedb_rule_4_3end_nt) %in% "u"
  
  mirgenedb_rule_4_2ntOffset_3endisU <- (mirgenedb_rule_4_2nt_5end_2ntOffset & mirgenedb_rule_4_2nt_3end_2ntOffset) | mirgenedb_rule_4_2nt_3end_isU
  
  summary(mirgenedb_rule_4_2nt_5end_2ntOffset)
  summary(mirgenedb_rule_4_2nt_3end_2ntOffset)
  summary(mirgenedb_rule_4_2nt_3end_isU)
  summary(mirgenedb_rule_4_2ntOffset_3endisU)
  
  
  # check for MirGeneDB rule 5: loop at least 8 nt long
  mirgenedb_rule_5_loopLength <- mirdeep2_res_db_seq_structure$loop_loc_end - mirdeep2_res_db_seq_structure$loop_loc_start - 1
  mirgenedb_rule_5_loopAtLeast8nt <- mirgenedb_rule_5_loopLength >= 8
  summary(mirgenedb_rule_5_loopAtLeast8nt)
  
  
  # check for MirGeneDB rule 6: mature RNA usually starts with "A" or "U"
  mirgenedb_rule_6_mature_nt <- ""
  for (i in 1:nrow(mirdeep2_res_db_seq_structure)) {
    mirgenedb_rule_6_mature_nt[i] <-
      substr(mirdeep2_res_db_seq_structure$consensus.mature.sequence[i], start = 1, stop = 1)
  }
  mirgenedb_rule_6_matureStartAorU <- grepl("a|u", x = mirgenedb_rule_6_mature_nt, ignore.case = TRUE)
  summary(mirgenedb_rule_6_matureStartAorU)
  
  # check for all MirGeneDB rules
  mirgenedb_rule_all <- mirgenedb_rule_1_bothArmsExpressed & 
    # mirgenedb_rule_2_5end_homogeneous & 
    mirgenedb_rule_3_imperfectComplementarity_atLeast16bp & mirgenedb_rule_4_2ntOffset_3endisU & mirgenedb_rule_5_loopAtLeast8nt & mirgenedb_rule_6_matureStartAorU
  summary(mirgenedb_rule_all)
  
  # SKIP mirgenedb_rule_2_5end_homogeneous first
  mirgenedb_rule_all_df <- data.frame(mirna_id, mirgenedb_rule_1_bothArmsExpressed, 
                                      mirgenedb_rule_2_5end_homogeneous,
                                      mirgenedb_rule_3_imperfectComplementarity, mirgenedb_rule_3_bpNumber,
                                      mirgenedb_rule_3_atLeast16bp, mirgenedb_rule_3_imperfectComplementarity_atLeast16bp,
                                      mirgenedb_rule_4_2nt_5end_2ntOffset, mirgenedb_rule_4_2nt_3end_2ntOffset,
                                      mirgenedb_rule_4_2nt_3end_isU, mirgenedb_rule_4_2ntOffset_3endisU,
                                      mirgenedb_rule_5_loopLength, mirgenedb_rule_5_loopAtLeast8nt,
                                      mirgenedb_rule_6_mature_nt, mirgenedb_rule_6_matureStartAorU, mirgenedb_rule_all)

  return(mirgenedb_rule_all_df)
} # function END
#===================================================

# runs only when script is run by itself
if (sys.nframe() == 0) {
  PARENT_DIR <- file.path(paste0("C:/Users/", Sys.info()[["user"]], "/Desktop/R/MirGeneDB_criteria_checker"))
  MIRDEEP2_RESULT_PATH = file.path(PARENT_DIR, "test/miRDeep2_miRNAs.norm_count.5p3p-arm-dominate.info.xls")
  
  mirdeep2_res_db_seq_structure <- read.csv(paste0(MIRDEEP2_RESULT_PATH, "_db_seq_structure.csv"))
  check_MirGeneDB_criteria(mirdeep2_res_db_seq_structure)
}
