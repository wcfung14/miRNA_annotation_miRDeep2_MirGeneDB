#' @title mirdeep2 results (.pdf) into .txt files converter
#' @description Given a pdf file path, convert it into text files for scrapping
#' @param pdf_file string; path to the pdf files to be converted
#' @param output_dir string; path to the folder for .txt output (default to be the same as *pdf_file*)
#' @export

mirdeep2_pdf_to_txt <- function(pdf_file, output_dir = dirname(pdf_file)) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  txt_file_path <- sub(".pdf$", ".txt", pdf_file)
  txt_file_path <- file.path(output_dir, basename(txt_file_path))

  cat(pdftools::pdf_text(pdf_file), file = txt_file_path)

  message(paste0("Finished converting from .pdf to .txt: ", basename(pdf_file)))
  message(paste0("output path: ", file.path(normalizePath(output_dir))))
} #function END
#===================================================

# runs only when script is run by itself
if (sys.nframe() == 0) {
  test_folder = paste0("C:/Users/",Sys.info()[["user"]],"/Desktop/R/MirGeneDB_criteria_checker/test")

  pdf_file_list = list.files(file.path(test_folder, "mirdeep2_pdf"), PATTERN = '.pdf$')
  
  for (i in pdf_file_list) {
    mirdeep2_pdf_to_txt(pdf_file = file.path(test_folder, "mirdeep2_pdf", i),
                    output_dir = file.path(test_folder, "mirdeep2_txt"))
    }
  }
 
