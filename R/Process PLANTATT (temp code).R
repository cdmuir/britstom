library(stringr)

text <- readLines(paste0(pathRawData, "/PLANTATT/test1.txt"))

# Remove " I "

text <- gsub(" I ", "", text)

# Space between columns

text <- gsub("(?<!,) ", "\t", text, perl = TRUE)

# Export

writeLines(text, con = paste0(pathProcData, "/PLANTATT/test1.txt"))
gregexpr("(?<!,) ", text, perl = TRUE)

str_match_all(text[1], "(?<!,) ")
word(text[1], start = 1, end = -1)

gregexpr("a", "abc abc")
