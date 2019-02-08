fasta_from_vcf <- function () {
  cat ("You will need a folder full of only .vcf files.\n")
  cat ("These will be all saved as one .fasta file, so keep populations\n")
  cat ("seperated and run fasta_from_vcf once per folder/population\n")
  cat ("\n")
  cat ("Enter the file path to your folder with a \"/\" at the end \n")
  cat ("e.g. \"/users/me/documents/r files/Asymptomatic VCF folder/\"\n")
  path <- readline (" My folder path <- ")
  path <- gsub("\"", "", path)

  cat ("What is the start position of the gene?\n")
  cat ("i.e. under \"POS\" in your vcf file, what number is the first base of your gene?\n")
  start.pos <- readline (" My gene start postion <- ")
  start.pos <- as.numeric(start.pos)

  cat (" \n")
  cat (" Enter the file path to your reference file\n")
  cat (" e.g. \"/users/me/documents/r files/reference.fasta\"\n")
  ref <- readline (" My ref .fasta file path <- ")
  ref <- gsub("\"", "", ref)

  cat ("Low quality SNPs under a set threshold can be excluded\n")
  cat ("Quality is determined by column \"QUAL\" in your vcf files\n")
  threshold <- readline (" My quality threshold <- ")
  threshold <- as.numeric(threshold)

  vcf.locs <- list.files(path, full.names = TRUE)
  vcf.names <- gsub("\\.[aA-zZ]*", "", list.files(path))
  vcf.num <- length(vcf.locs)

  ref.string <- toupper(as.matrix(paste(readLines(ref)[-1], collapse = "")))
  seq.length <- sum(nchar(ref.string))
  ref.matrix <- matrix(nrow = 1, ncol = seq.length)
  for (i in 1:seq.length){ref.matrix[1,i] <- substr(ref.string, i, i)}

  cat("\n")
  my.vcf.2.fasta <- matrix()
  for (i in 1:vcf.num){
    try({
      vcf.as.table <- read.table(vcf.locs[i])
      vcf.as.table <- cbind(as.numeric(vcf.as.table[ ,2]),
                          vcf.as.table[ ,4],
                          vcf.as.table[ ,5],
                          as.numeric(vcf.as.table[ ,6]))

      vcf.as.table[ ,1] <- vcf.as.table[ ,1] - start.pos + 1


      vcf.as.table[ ,2:3] <- gsub(1, "A", vcf.as.table[ ,2:3])
      vcf.as.table[ ,2:3] <- gsub(2, "C", vcf.as.table[ ,2:3])
      vcf.as.table[ ,2:3] <- gsub(3, "G", vcf.as.table[ ,2:3])
      vcf.as.table[ ,2:3] <- gsub(4, "T", vcf.as.table[ ,2:3])

      colnames(vcf.as.table) <- c("pos", "alt", "alt", "qual")

      if (ref.matrix[as.numeric(vcf.as.table[,1])][1] == vcf.as.table[1,2])
        colnames(vcf.as.table)[2] <- c("ref")
      if (ref.matrix[as.numeric(vcf.as.table[,1])][1] != vcf.as.table[1,2])
        colnames(vcf.as.table)[3] <- c("ref")

      vcf.as.table <- vcf.as.table[as.numeric(vcf.as.table[ ,4]) > threshold, ]
      vcf.as.table <- as.data.frame(vcf.as.table)
      seq.matrix <- ref.matrix
      for (j in 1:nrow(vcf.as.table)){
        seq.matrix[1,as.numeric(as.character(vcf.as.table$pos[j]))] <- as.character(vcf.as.table$alt[j])
      }
      my.vcf.2.fasta <- rbind(my.vcf.2.fasta, paste(">", vcf.names[i]), paste(seq.matrix, collapse = ""))
      cat("Reading .vcf", i, "of", vcf.num, "\r")
    }, silent = TRUE)
  }

  my.vcf.2.fasta <- as.matrix(my.vcf.2.fasta[-1,])
  cat("\n")
  cat("Finished!\n")
  folder.name <- substr(path, 1, nchar(path)-1)
  folder.name <- gsub("..*[aA-zZ]*/", "", folder.name)

  cat("What would you like to do with the .fasta file \n")
  cat("  1 - Save as a .fasta file in your R working directory\n")
  cat("      WARNING - this will create a new file (your .vcf files will remain unchanged)\n")
  cat("  2 - Save the .fasta table as an object (my.vcf.2.fasta) in your R environment\n", sep = "")
  choice <- readline (" Enter 1 or 2 <- ")

  if (choice == 1)
    {write.table(my.vcf.2.fasta, file = paste(folder.name, ".fasta", sep = ""), quote = FALSE,
                 row.names = FALSE, col.names = FALSE)
    cat("Saved as ", folder.name, ".fasta in your working directory", sep = "")}
  if (choice == 2)
    {my.vcf.2.fasta <<- my.vcf.2.fasta
    cat("Saved as my.vcf.2.fasta in your global environment")}

}
