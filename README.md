# miRNA_annotation_miRDeep2_MirGeneDB
Analyze microRNAs prediction by miRDeep2 and check the predictions by MirGeneDB criteria for miRNAs annotation. 

Workflow: parsing mirdeep2 results, taking miRNA dot-bracket notations, mapping to miRNA sequences and obtain miRNA structure and hairpin information, make miRNA hairpins by aligning top and bottom strand, and checking for MirGeneDB "Unique structural features of microRNAs".
Require packages: `pdftools`, `stringr`, `stringi`

Workflow:
1.	Initialize environment parameters
2.	Convert mirdeep2 results (pdf files) into txt files
a.	`"mirdeep2_pdf_to_txt.R”`
3.	Scrap dot-bracket notation from txt files
a.	`"dot_bracket_notation_scrapper.R"`
4.	Merge mirdeep2 result (.xlsx) and dot-bracket notation into one dataframe
5.	Map microRNAs dot-bracket notation to miRNA sequences and return miRNA structure information (mature, star and loop sequences)
a.	`"dot_bracket_notation_seq_mapper.R"`
6.	Align microRNAs hairpin dot-bracket notation
a.	`“dot_bracket_notation_aligner.R"`
b.	To view miRNA dot-bracket notation structure, open the csv, set font to monospaced typefaces (e.g. Consolas)
7.	Check for MirGeneDB "Unique structural features of microRNAs"
a.	`"MirGeneDB_criteria_checker.R"`

Known limitations:
1.	Alignment does not work if there are more than one loop found in miRNA sequences (will skip intentionally) (e.g. Tco_Scaffold_6695_20322, Hho_SczTNLB_6657_30897) 
2.	Does not check for Rule 2 of MirGeneDB "Unique structural features of microRNAs"
3.	miRNA hairpin structure in mirdeep2 pdf files is not the same as that of dot-bracket notation (e.g. Tco_Scaffold_9089_33841)
4.	4. The "obs" miRNA sequence in mirdeep2 results is used in this analysis. The "exp" sequence (with Cyan color, from mirdeep2 prediction, to fulfill criteria) and its additional nucleotide is not included as part of "consensus.mature.sequence" or "consensus.star.sequence", and will not be used for analysis by this package.


## script
- `"mirdeep2_pdf_to_txt.R”` - Convert mirdeep2 results (pdf files) into txt files
- `"dot_bracket_notation_scrapper.R"` - Scrap dot-bracket notation from txt files
- `"dot_bracket_notation_seq_mapper.R"` - Map microRNAs dot-bracket notation to miRNA sequences and return miRNA structure information (mature, star and loop sequences)
- `“dot_bracket_notation_aligner.R"` - Align microRNAs hairpin dot-bracket notation
- `"MirGeneDB_criteria_checker.R"` - Check for MirGeneDB "Unique structural features of microRNAs"
