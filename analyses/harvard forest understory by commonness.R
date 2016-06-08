d <- read.csv("/Users/danflynn/Desktop/hf106-03-shrub.csv")


common <- sort(apply(d[4:ncol(d)], 2, sum)/sum(apply(d[4:ncol(d)], 2, sum))*100, T)

sp <- read.csv("/Users/danflynn/Desktop/hf106-01-species-codes.csv")





sp[match(toupper(names(common)[1:25]), sp$code),]
