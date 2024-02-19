# Function to remove special characters from names
clean_text <- function(character_vector, remove_special=T, remove_space=T, lower=T) {
  
  if(remove_space == T) {
    # remove leading and trailing spaces
    cleaned_names <- trimws(character_vector, which = "both")
    writeLines("Removed leading and trailing spaces.\n")
  }
  
  if(remove_special == T) {
    # remove special characters
    cleaned_names <- gsub("[^A-Za-z0-9 ]", "", cleaned_names)
    writeLines("Removed special characters.\n")
  }
  
  if(lower==T) {
    cleaned_names <- tolower(cleaned_names)
    writeLines("Characters changed to lowercase.\n")
  }
  
  return(cleaned_names)
}