# function extracts the number of rows of the model file that contain a parameter of interest
grep.fun <- function(parameter.name, model.file.text){
  ##'@param parameter.name - the name of a parameter of interest
  ##'@param model.file.text - the text of a model file parsed as a character vector
  
  # The rows with the parameter names included 
  parameter.test <- paste("\\b",parameter.name,"\\b",sep="")
  positions <- grep(parameter.test, model.file.text, value = FALSE)
  if(length(positions)==0){  positions <- grep(paste("\\b",parameter.name,sep=""),
                                             model.file.text, value = FALSE)}
  
  return(length(positions))
}