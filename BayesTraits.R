# Run Bayestrait for all families for which input files are available
# it was verified that all datfiles have at least 3 lines (without
# having to delete any of them)

dir.create("BT_outfiles")
d <- dir("datfiles")
cat("13\n2\nRun\n", file="BT_instructions.txt")
system.command <- "BayesTraitsV4.exe tmp.nex tmp.dat < BT_instructions.txt"
cat(system.command, file="BTrun.bat")
get.family.name <- function(x) strsplit(x, ".dat")[[1]][1]
fams <- unlist(lapply(dir("datfiles"), get.family.name))
for (i in seq_along(fams)) {
  cat("doing", fams[i], "which is no.", i, "in the list\n")
  nexus.file <- paste0("nexfiles/", fams[i], ".nex")
  data.file <- paste0("datfiles/", fams[i], ".dat")
  file.copy(nexus.file, "tmp.nex", overwrite=TRUE)
  file.copy(data.file, "tmp.dat", overwrite=TRUE)
  system("BTrun")
  file.rename("tmp.dat.AncStates.txt", paste0("BT_outfiles/", fams[i], "_AncStates.txt"))
  file.remove("tmp.dat.Log.txt", "tmp.dat.Schedule.txt", "tmp.dat", "tmp.nex")
}
file.remove("BT_instructions.txt")
file.remove("BTrun.bat")

# for experimenting
nexus.file <- paste0("nexfiles/", "new_tree", ".nex")
data.file <- paste0("datfiles/", fams[i], ".dat")
file.copy(nexus.file, "tmp.nex", overwrite=TRUE)
file.copy(data.file, "tmp.dat", overwrite=TRUE)
system("BTrun")

old_tree <- read.nexus(paste0("nexfiles/", fams[i], ".nex"))
new_tree <- di2multi(old_tree, tol=.1)
write.nexus(new_tree, file="new_tree.nex")

new_tree <- remove.zero.brlen(old_tree)
library(dispRity)



