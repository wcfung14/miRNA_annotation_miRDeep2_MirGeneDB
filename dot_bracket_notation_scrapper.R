#' @title Dot-bracket notation from mirdeep2 results .txt files scrapper
#' @description Read and scrap microRNAs dot-bracket notation 
#' from mirdeep2 results .txt files (converted from .pdf using "mirdeep2_pdf_to_txt.R").
#' @param txt_file string; .txt file to be scrapped. Obtained after converting  mirdeep2 results (.pdf) into .txt files ("mirdeep2_pdf_to_txt.R")
#' @return vector of three elements: miRNA ID ("mirna_id"), its dot-bracket notation ("dot_bracket_notation_trimmed"), and its sequence ("mirna_seq_trimmed")
#' @export

scrap_dot_bracket_notation <- function(txt_file) {
  mirna_id <- sub(".txt$", "", basename(txt_file))
  
  txt_file_r <- readLines(txt_file)
  
  dot_bracket_notation <- txt_file_r[grep("\\(\\.+\\)", txt_file_r)] # grep dot-bracket notation by structure that looks like "(....)"
  dot_bracket_notation_trimmed <- gsub(" .*$", "", stringr::str_trim(stringr::str_squish(dot_bracket_notation))) # removing everything after the first space (squished tab)
  
  mirna_seq <- txt_file_r[grep("5'-", txt_file_r)] #grep miRNA sequence by "5' -"
  mirna_seq_trimmed <- stringr::str_trim(stringr::str_squish(sub("5'-(.*)-3'.*", "\\1", mirna_seq)))
  
  output <- c(mirna_id, dot_bracket_notation_trimmed, mirna_seq_trimmed)
  return(output)
} # function END
#===================================================

# runs only when script is run by itself
if (sys.nframe() == 0) {
  path = paste0("C:/Users/",Sys.info()[["user"]],"/Desktop/R/MirGeneDB_criteria_checker/test/mirdeep2_txt")
  file <- file.path(path, "Scaffold_1043_1001.txt")
  
  res <- dot_bracket_notation_scrapper(txt_file = file)
  print(res)
}
