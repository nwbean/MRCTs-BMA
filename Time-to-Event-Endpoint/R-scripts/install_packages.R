#################################################################
# Install R Packages Using Bash Script
#################################################################

pkg_list <- c("simsurv")

for(pkg in pkg_list){
  install.packages(pkg, repos = "http://cran.us.r-project.org", 
                   lib = "../R-packages") 
  if(!require(pkg, character.only = TRUE)){
    quit(save = "no", status = 1, runLast = FALSE)  
  }
}
