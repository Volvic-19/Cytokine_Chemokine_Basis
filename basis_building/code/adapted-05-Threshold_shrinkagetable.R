## This script is to inspect data in shrinkage.e5.DT.tsv.gz, and further filtering shrinkage table with threshold: e4, e3,e2
# 2022-07-11

library(data.table)
setDTthreads(10)

# Load e5 data
e5 <- fread("../manifest/shrinkage.e5.DT.tsv.gz", tmpdir = "tmp") #25643160 = 641079*40
 


# Filter shrinkage value by e4, e3, e2

e4 <- e5[shrinkage > 1e-4]
e3 <- e5[shrinkage > 1e-3]
e2 <- e5[shrinkage > 1e-2]


length(unique(e5$pid))
[1] 641079
length(unique(e4$pid))
[1] 59290
length(unique(e3$pid))                                                        
[1] 6434
length(unique(e2$pid))
[1] 430

fwrite(e4, "../manifest/shrinkage.e4.DT.tsv.gz", sep="\t")
fwrite(e3, "../manifest/shrinkage.e3.DT.tsv.gz", sep="\t")
fwrite(e2, "../manifest/shrinkage.e2.DT.tsv.gz", sep="\t")


