rm(list=ls(all=TRUE))

# Supondo que o arquivo de entrada não possui espaços e nem quebra de linhas.
print('file name:')
file_name<-scan(, what="")

# A entrada depth deve ser um inteiro
print('depth:')
depth<-as.numeric(scan(,what=""))
sample<-scan(file_name,what=character())

# O setup_BIC.R deve estar no diretório descrito abaixo:
source("~/Documents/setup_BIC.R")
tree <- setup_BIC(sample,depth)