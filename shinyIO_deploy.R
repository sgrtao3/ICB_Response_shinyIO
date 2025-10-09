#install.packages("rsconnect")
options(pkgType = "binary", timeout = 600)
install.packages("rsconnect", repos = "https://cloud.r-project.org", type = "win.binary")



setwd("G:/Ranran/ICB_Response/shinyapp")   # <- change this
list.files()
list.files("data")
rsconnect::bundleApp()  

# run once locally
expr <- readr::read_tsv("G:/Ranran/ICB_Response/combined_expr.tsv")
saveRDS(expr, "G:/Ranran/ICB_Response/shinyapp/data/combined_expr.rds")
x <- readRDS("G:/Ranran/ICB_Response/shinyapp/data/combined_expr.rds")
x
stopifnot("gene" %in% names(x))
mat <- as.matrix(x[ , setdiff(names(x), "gene"), drop = FALSE])
rownames(mat) <- as.character(x$gene)
saveRDS(mat, "G:/Ranran/ICB_Response/shinyapp/data/combined_expr.rds") # <- new compact file

library(rsconnect)
rsconnect::setAccountInfo(name='ranrantao', token='624BFCF0F1E5BB6F5F2BC24F05733CAD', secret='kFXney96Anq5vJyHNDcnrDxY5tax1FZIfPSqf41r')
library(rsconnect)
# optional: give your app a friendly name (used in the URL)
rsconnect::deployApp("G:/Ranran/ICB_Response/shinyapp",appName = "ICB_Response")
