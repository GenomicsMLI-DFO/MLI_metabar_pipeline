

# Function to extract important info from WORMS
# With the help of this demo
# http://www.marinespecies.org/aphia/webservice/Aphia_webservice_R_elaborate.txt

#install the required packages (comment if needed)
#install.packages("jsonlite", repos="http://cran.r-project.org")
#install.packages("httr")

#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)


#namesToMatch <- c("Anonyx nugax", "Salmo salar")

validateWORMS <- function(namesToMatch){
  
  #Convert the namesToMatch to a valid REST-url
  urlNamesPart <- ""
  for (index in 1:length(namesToMatch)) {
    urlNamesPart <- sprintf("%s&scientificnames[]=%s", urlNamesPart, namesToMatch[index]);
  }

  #The url can contain special characters that need to be converted
  urlNamesPart <- URLencode(urlNamesPart)

  #The dyanmic build of the URL causes an obsolete '&' at the beginning of the string, so remove it
  urlNamesPart <- substring(urlNamesPart, 2)
  # Add this part to allow non-marine species
  urlSearchPart <- paste0(urlNamesPart,"&marine_only=false")
  #Build the final REST-url
  url <- sprintf("http://www.marinespecies.org/rest/AphiaRecordsByMatchNames?%s", urlSearchPart)

  #Get the actual data from the URL
  matches <- fromJSON(url)

  res <- data.frame(Species = character(),
                    scientificname = character(),
                    rank = character(),
                    AphiaID = numeric(),
                    authority = character(),
                    status = character(),
                    valid_name = character(),

                    # add taxonomic info
                    kingdom = character(),
                    phylum = character(),
                    class = character(),
                    order = character(),
                    family = character(),
                    genus = character(),
                    n_results = numeric(),
                    
                    stringsAsFactors = F)

  # Function to handle more than one entry within the results
  
  fct.dup <- function(vec){
      res <- vec %>% unique() %>% paste(collapse = " | ")
      return(res)
  }
  
  
  #Handle the data (each requested name has an list of results)
for (matchesindex in 1:length(namesToMatch)) {
  #Get the results for the current index
  currentResultList <- matches[[matchesindex]]
  
  #Get the number of list entries for the first column
  numberOfResults <- nrow(currentResultList)
  
  #Handle empty data due to no matches found
  #if (is.na(currentResultList[[1]][[1]])) {
  #  numberOfResults <- 0
  #}
  
  if (numberOfResults > 0) {

    res[matchesindex,] <- c(namesToMatch[matchesindex], 
                            currentResultList$scientificname %>% fct.dup(),
                            currentResultList$rank %>% fct.dup(),
                            currentResultList$AphiaID %>% fct.dup(),
                            currentResultList$authority %>% fct.dup(),
                            currentResultList$status %>% fct.dup(),
                            currentResultList$valid_name %>% fct.dup(),
                            

                            # check if these are ok BEFORE running
                            currentResultList$kingdom %>% fct.dup(),
                            currentResultList$phylum %>% fct.dup(),
                            currentResultList$class %>% fct.dup(),
                            currentResultList$order %>% fct.dup(),
                            currentResultList$family %>% fct.dup(),
                            currentResultList$genus %>% fct.dup(),
                            numberOfResults
    )   
    
      
       
  } else  {
    
    res[matchesindex,] <- c(namesToMatch[matchesindex], 
                            rep(NA, ncol(res)-2),
                            numberOfResults)
      
    
    }
  }

  return(res)

  }






