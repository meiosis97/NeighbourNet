# Install Package:           'Ctrl + Shift + B'
# Check Package:             'Ctrl + Shift + E'
# Test Package:              'Ctrl + Shift + T'

# usethis::use_test("name")

# Install chocolate using windows powershell
# Rcmd choco install git.install

# usethis::use_git_config(user.name = "Yidi Deng", user.email = "meiosis97@gmail.com")
# usethis::create_github_token()
# git config --global --list
# git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git (Generate a colon to the local machine)

# cd myrepo (Treat the github repository as a directory)
# echo "A line I wrote on my local computer  " >> README.md (An example modification of the repository)
# git status (Check what modifications has been made)
# git add README.md (Add modifications that hopes to be committed on the local machine.)
# git commit -m "A commit from my local computer" (Commit the added changes to the local machine)
# git push (Synchronize the local modification to the github repository)

# Export data
# 1. load data
# 2. usethis::use_data(my_data)
# 3. write documentation https://r-pkgs.org/data.html

# Gitignore does not ignore existing file in the repository.

#' @import Matrix
NULL

#' @import Seurat
NULL

#' @importFrom magrittr `%>%`
NULL

