source("header.R")

writeLines(capture.output(session_info()), str_c("R/session-info-", user(),".txt"))
