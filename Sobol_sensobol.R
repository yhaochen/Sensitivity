rm(list = ls())
graphics.off()
library(sensobol)
#Problem definition: return yes if the product of two numbers >=0.75
Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
}
para<-c("X1","X2")
N<-500000
R<-100
order <- "first"
type <-"LHS"#"QRN" "R"
matrices=c("A", "B","BA")
mat <- sobol_matrices(N = N, params = para,order = order,type="LHS",matrices = matrices)
y<-Reliability(mat)
ind <- sobol_indices(Y = y, N = N, params = para, order=order,
                     boot = TRUE, R = R, type = "norm", conf = 0.95)
