# HF tree traits summary
library(gdata)
d <- read.xls("~/Dropbox/Work/Harvard/Farm and Forest/HF Tree Traits 2015-06-09.xlsx")
summary(d)

d <- d[d$Summer.2015.Route != "",]

dt <- table(d$Species)

sort(dt, T)
dt[dt<4]

SpeciesList <- unique(d$Species)

length(SpeciesList)

# grant planning

