library(biodiversity)
data(species)
suppressWarnings({
cat(logLik(biodiv(species, conf=1)),"\n", file="res.out")
cat(logLik(biodiv(species, conf=2)),"\n", file="res.out", append=TRUE)
cat(logLik(biodiv(species, conf=3)),"\n", file="res.out", append=TRUE)
cat(logLik(biodiv(species, conf=4)),"\n", file="res.out", append=TRUE)
})
