# Function to extract stuffs from options


get.value <- function(NAME){
  OPTIONS <- readr::read_lines("Options.txt")
  NAME <- paste0(NAME, ":")
  res <-stringr::str_subset(OPTIONS, NAME) 
  
  res  <- stringr::str_remove(res, NAME)  
  res <- stringr::str_remove_all(res, " ")
  
  return(res)  
}


