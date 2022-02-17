# Function to extract stuffs from options


get.value <- function(NAME, file = "Options.txt"){
  OPTIONS <- readr::read_lines(file)
  NAME <- paste0(NAME, ":")
  res <-stringr::str_subset(OPTIONS, NAME) 
  
  res  <- stringr::str_remove(res, NAME)  
  res <- stringr::str_remove_all(res, " ")
  
  return(res)  
}



