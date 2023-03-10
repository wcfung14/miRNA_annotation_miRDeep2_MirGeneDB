#' @title microRNAs annotation by mirdeep2 and MirGeneDB criteria 
#' @description Analyze microRNAs prediction by miRDeep2 and check the predictions by MirGeneDB criteria for miRNAs annotation. 
#' Workflow: parsing mirdeep2 results, taking miRNA dot-bracket notations, mapping to miRNA sequences and obtain miRNA structure and hairpin information, make miRNA hairpins by aligning top and bottom strand, and checking for MirGeneDB "Unique structural features of microRNAs".
#' Require packages: pdftools, stringr, stringi
#' Known limitations:
#' 1. Alignment does not work if there are more than one loop found in miRNA sequences (will skip if "skip_this" == TRUE when running "dot_bracket_notation_seq_mapper.R") (e.g. Tco_Scaffold_6695_20322, Hho_SczTNLB_6657_30897) 
#' 2. Does not check for Rule 2 of MirGeneDB "Unique structural features of microRNAs"
#' 3. Some projected miRNA hairpin structures in mirdeep2 pdf files are not the same as that of dot-bracket notation (e.g. Tco_Scaffold_9089_33841)
#' 4. The "obs" miRNA sequence in mirdeep2 results is used in this analysis. The "exp" sequence (with Cyan color, from mirdeep2 prediction, to fulfill criteria) and its additional nucleotide is not included as part of "consensus.mature.sequence" or "consensus.star.sequence", and will not be used for analysis by this package.

  # 1. Initialize environment parameters
{
  PARENT_DIR <- file.path(paste0("C:/Users/", Sys.info()[["user"]], "/Desktop/R/miRNA_annotation_miRDeep2_MirGeneDB"))
  setwd(PARENT_DIR)
  PDF_FILE_DIR = file.path(PARENT_DIR, "test/mirdeep2_pdf")
  TXT_FILE_DIR = file.path(PARENT_DIR, "test/mirdeep2_txt")
  MIRDEEP2_RESULT_PATH = file.path(PARENT_DIR, "test/miRDeep2_miRNAs.norm_count.5p3p-arm-dominate.info.xls") 
}

# 2. Convert mirdeep2 results (.pdf) into .txt files
source("mirdeep2_pdf_to_txt.R")
pdf_file_list <- list.files(PDF_FILE_DIR, full.names = TRUE, pattern = '.pdf$')
for (i in pdf_file_list) {
  mirdeep2_pdf_to_txt(pdf_file = i, output_dir = TXT_FILE_DIR)
}

# 3. Scrap dot-bracket notation from .txt files
# return a vector of three elements ("mirna_id", "dot_bracket_notation" and "mirna_seq"), and merging them into a dataframe
source("dot_bracket_notation_scrapper.R")
dot_bracket_df <- data.frame()
txt_file_list <- list.files(TXT_FILE_DIR, full.names = TRUE, pattern = '.txt$')
for (i in txt_file_list) {
  res <- scrap_dot_bracket_notation(txt_file = i)
  dot_bracket_df <- rbind(dot_bracket_df, res)
}
colnames(dot_bracket_df) <- c("mirna_id", "dot_bracket_notation", "mirna_seq")

# 4. Merge mirdeep2 result (.xlsx) and dot-bracket notation into one dataframe
mirdeep2_res <- read.delim(MIRDEEP2_RESULT_PATH, header = TRUE, fill = TRUE, na.strings = "") 

mirdeep2_res_db <- merge(mirdeep2_res, dot_bracket_df, by.x="provisional.id", by.y="mirna_id", all.x = TRUE)
mirdeep2_res_db <- subset(mirdeep2_res_db, select = -X.miRNA)
summary(!is.na(mirdeep2_res_db$dot_bracket_notation)) # return number of miRNAs identified by mirDeep2 with dot-bracket notation

write.csv(mirdeep2_res_db, paste0(MIRDEEP2_RESULT_PATH, "_db.csv"), row.names = FALSE)

# 5. Map microRNAs dot-bracket notation to miRNA sequences and return miRNA structure information (mature, star and loop sequences)
# return a dataframe of 14 columns of miRNA structure information: ("skip_this", "mature_loc_start", "mature_loc_end", "star_loc_start", "star_loc_end", "loop_loc_start", "loop_loc_end", "mature_first", "hairpin_db", "mature_db", "star_db", "extension_db", "loop_db", "Remarks")
# Alignment does not work if there are more than one loop found in miRNA sequences. Will skip if "skip_this" == TRUE
source("dot_bracket_notation_seq_mapper.R")
mirna_map_df <- data.frame()

for (i in 1:nrow(mirdeep2_res_db)) {
  if (is.na(mirdeep2_res_db$mirna_seq[i])) {next}
  res <- map_dot_bracket_notation_seq(mirna_seq = mirdeep2_res_db$mirna_seq[i], # miRNA sequence from mirdeep2 .pdf results. (Has the same length with mirna dot_bracket_notation, not mirdeep2_res_db_seq$consensus.precursor.sequence)
                                      mirna_mature_seq = mirdeep2_res_db$consensus.mature.sequence[i], # required
                                      mirna_star_seq = mirdeep2_res_db$consensus.star.sequence[i], # required
                                      dot_bracket_notation = mirdeep2_res_db$dot_bracket_notation[i])
  mirna_id = mirdeep2_res_db$provisional.id[i]
  res <- cbind(mirna_id, res) # add mirna_id for merging
  
  mirna_map_df <- rbind(mirna_map_df, res)
}
mirdeep2_res_db_seq <- merge(mirdeep2_res_db, mirna_map_df, by.x="provisional.id", by.y="mirna_id", all.x = TRUE)
write.csv(mirdeep2_res_db_seq, paste0(MIRDEEP2_RESULT_PATH, "_db_seq.csv"), row.names = FALSE)

# 6. Align microRNAs hairpin dot-bracket notation. Optional but needed for MirGeneDB check: map miRNA sequence to dot-bracket notation.
# return a dataframe of 8 columns of miRNA hairpin structure information: ("top_strand_align", "strand_match", "bottom_strand_align", "hairpin_structure_db", "Remarks", "bottom_strand_align_seq", "hairpin_structure_seq")
source("dot_bracket_notation_aligner.R")
hairpin_structure_df <- data.frame()

for (i in 1:nrow(mirdeep2_res_db_seq)) {
  if (is.na(mirdeep2_res_db_seq$hairpin_db[i])) {next}
  res <- align_dot_bracket(hairpin_db = mirdeep2_res_db_seq$hairpin_db[i], 
                           align_from_loop = TRUE,
                           mirna_precursor_seq = mirdeep2_res_db_seq$consensus.precursor.sequence[i], 
                           mirna_mature_seq = mirdeep2_res_db_seq$consensus.mature.sequence[i],
                           mirna_star_seq = mirdeep2_res_db_seq$consensus.star.sequence[i])
  mirna_id = mirdeep2_res_db_seq$provisional.id[i]
  res <- cbind(mirna_id, res) # add mirna_id for merging
  
  hairpin_structure_df <- rbind(hairpin_structure_df, res)
}
mirdeep2_res_db_seq_structure <- merge(mirdeep2_res_db_seq, hairpin_structure_df, by.x="provisional.id", by.y="mirna_id", all.x = TRUE)
write.csv(mirdeep2_res_db_seq_structure, paste0(MIRDEEP2_RESULT_PATH, "_db_seq_structure.csv"), row.names = FALSE)
# To view miRNA dot-bracket notation structure, open the csv, set font to monospaced typefaces (e.g. Consolas)

# 7. Check for MirGeneDB "Unique structural features of microRNAs"
# return a dataframe of 16 columns of MirGeneDB criteria check("mirna_id", "mirgenedb_rule_1_bothArmsExpressed", "mirgenedb_rule_2_5end_homogeneous", "mirgenedb_rule_3_imperfectComplementarity", "mirgenedb_rule_3_bpNumber", "mirgenedb_rule_3_atLeast16bp", "mirgenedb_rule_3_imperfectComplementarity_atLeast16bp", "mirgenedb_rule_4_2nt_5end_2ntOffset", "mirgenedb_rule_4_2nt_3end_2ntOffset", "mirgenedb_rule_4_2nt_3end_isU", "mirgenedb_rule_4_2ntOffset_3endisU", "mirgenedb_rule_5_loopLength", "mirgenedb_rule_5_loopAtLeast8nt", "mirgenedb_rule_6_mature_nt", "mirgenedb_rule_6_matureStartAorU", "mirgenedb_rule_all")
#' IMPORTANT: Rule 2 is NOT checked currently, and mirgenedb_rule_all does not check for rule 2
source("MirGeneDB_criteria_checker.R")
mirgenedb_rule_all_df <- check_MirGeneDB_criteria(mirdeep2_res_db_seq_structure)
mirdeep2_res_db_seq_structure_checked <- merge(mirdeep2_res_db_seq_structure, mirgenedb_rule_all_df, by.x="provisional.id", by.y="mirna_id", all.x = TRUE)
summary(mirgenedb_rule_all_df)

write.csv(mirdeep2_res_db_seq_structure_checked, paste0(MIRDEEP2_RESULT_PATH, "_db_seq_structure_checked.csv"), row.names = FALSE)
