vaxpack_input <- function () {
  cat ("The following analyses are intended for haploid organisms such as Plasmodium")
  cat ("\n")
  cat ("For this function to work you will need:\n")
  cat ("1 - A folder containing all files to be analysed\n")
  cat ("     These need to be aligned, the same length, and with no gaps! (e.g. with MEGA)\n")
  cat ("     Only accepts bases 'A/a', 'T/t', 'C/c', 'G/g' \n")
  cat ("     Gaps, or bases that are not AaTtCcGg, are called as reference, so accuracy is lost\n")
  cat ("     These replacements will be counted in the results table as \"invalid sites\"\n")
  cat ("\n")
  cat ("     Make a different file for different populations to compare them\n")
  cat ("     e.g. 'Asymptomatic.fasta, Symptomatic.fasta', or 'Brazil.fasta, Peru.fasta'\n")
  cat ("\n")
  cat ("2 - A reference file for the gene you are analysing\n")
  cat ("     (This cannot be in the same folder as the other files!)\n")
  cat ("\n")
  cat ("Accepted extensions are \".fasta\", \".fas\", \".fa\", and \".seq\"\n")
  cat ("\n")
  cat ("Enter the file path to your folder with a \"/\" at the end \n")
  cat ("e.g. \"/users/me/documents/r files/my fasta folder/\"\n")
  path <- readline (" My folder path <- ")
  path <- gsub("\"", "", path)

  if (substr(path, nchar(path), nchar(path)) != "/")
    stop("There needs to be a backslash at the end!")

  cat (" \n")
  cat (" Enter the file path to your reference file\n")
  cat (" e.g. \"/users/me/documents/r files/reference.fasta\"\n")
  ref <- readline (" My ref .fasta file path <- ")
  ref <- gsub("\"", "", ref)

  if (grepl(".fasta", ref, ignore.case = TRUE) == FALSE)
    if (grepl(".fas", ref, ignore.case = TRUE) == FALSE)
      if (grepl(".fa", ref, ignore.case = TRUE) == FALSE)
        if (grepl(".seq", ref, ignore.case = TRUE) == FALSE)
          stop("The reference file extension must be \".fasta\", \".fas\", \".fa\", or \".seq\"")

  cat (" \n")
  genename <- readline (" What is the name of the gene you're analysing? <- ")
  cat (" \n")




  #vaxpack
  #---------------------------------------------------------------------------------------------------
  st <- Sys.time()
  # processing of fasta files ---------------------------------------------------------------------------------------------------
  fasta.locs <- list.files(path, full.names = TRUE)
  fasta.names <- gsub("\\.[aA-zZ]*", "", list.files(path))
  fasta.num <- length(fasta.locs) #number of files inside folder

  cat (" Completing Step: 1 of", (fasta.num * 7) + 16, "\r")
  # building matrix for result table : option-1 --------------------
  ndrt <- matrix(nrow = fasta.num + 1, ncol = 11)
  rownames(ndrt) <- c(fasta.names, "Total")
  colnames(ndrt) <- c("n", "Seq. Length", "Invalid Sites", "Seg. sites (S)", "SNPS", "\u03a0 x 10^-3", "TD",
                      "NS", "SP", "nuc. h", "nuc. Hd")

  # prepare for ref to input comparison (consider everything as string)---------------------------------------------------------------------------------------------------
  cat (" Completing Step: 2 of", (fasta.num * 7) + 16, "\r")
  ref.string <- toupper(as.matrix(paste(readLines(ref)[-1], collapse = "")))
  seq.length <- sum(nchar(ref.string)) #length of ref
  ref.matrix <- matrix(nrow = 1, ncol = seq.length) #set ref seq to compare with input fasta as string, each base become a column
  totalpop <- matrix(ncol = seq.length) #make a matrix table
  for (i in 1:seq.length){ref.matrix[1,i] <- substr(ref.string, i, i)} #substr — Return part of a string
  # build the table for AA convsersion--------------------------------------------------------------------------------------------------
  cat (" Completing Step: 3 of", (fasta.num * 7) + 16, "\r")
  aalookup <-matrix(c("ATT", "I", "ATC", "I", "ATA", "I",
                      "CTT", "L", "CTC", "L", "CTA", "L", "CTG", "L", "TTA", "L", "TTG","L",
                      "GTT", "V", "GTC", "V", "GTA", "V", "GTG", "V",
                      "TTT", "F", "TTC", "F",
                      "ATG", "M",
                      "TGT", "C", "TGC", "C",
                      "GCT", "A", "GCC", "A", "GCA", "A", "GCG", "A",
                      "GGT", "G", "GGC", "G", "GGA", "G", "GGG", "G",
                      "CCT", "P", "CCC", "P", "CCA", "P", "CCG", "P",
                      "ACT", "T", "ACC", "T", "ACA", "T", "ACG", "T",
                      "TCT", "S", "TCC", "S", "TCA", "S", "TCG", "S", "AGT", "S", "AGC", "S",
                      "TAT", "Y", "TAC", "Y",
                      "TGG", "W",
                      "CAA", "Q", "CAG", "Q",
                      "AAT", "N", "AAC", "N",
                      "CAT", "H", "CAC", "H",
                      "GAA", "E", "GAG", "E",
                      "GAT", "D", "GAC", "D",
                      "AAA", "K", "AAG", "K",
                      "CGT", "R", "CGC", "R", "CGA", "R", "CGG", "R", "AGA", "R", "AGG", "R",
                      "TAA", "STOP", "TAG", "STOP", "TGA", "STOP"), nrow = 64, ncol = 2, byrow = TRUE)
  #---------------------------------------------------------------------------------------------------
  loop.count <- 1
  invalid.site.counter.total <- 0
  for (i in 1:fasta.num){
    cat (" Completing Step:", loop.count + 3, "of", (fasta.num * 7) + 16, "\r") #keep track
    onepop <- toupper(as.matrix(readLines(fasta.locs[i])[c(FALSE, TRUE)])) #reading input fasta
    rownames(onepop) <- readLines(fasta.locs[i])[c(TRUE, FALSE)] #sample ID
    seq.count <- nrow(onepop) #number of seq inside aligned fasta
    if (length(unique(nchar(onepop))) != 1)
      stop("All sequences must be the same length!")
    if ((unique(nchar(onepop))) != seq.length)
      stop("Reference is not the same length as your sequences!") #quality control
    #---------------------------------------------------------------------------------------------------
    cat (" Completing Step:", loop.count + 4, "of", (fasta.num * 7) + 16, "\r")
    invalid.site.counter <- 0
    if (seq.count > 1){
      wholegene <- matrix(nrow = seq.count, ncol = seq.length) #set matrix
      n.change.matrix <- wholegene
      for (j in 1:seq.count){
        splitgene <- strsplit(onepop[j,1], split = "")[[1]] #split base from fasta for each sample
        for (k in 1:seq.length){
          if (splitgene[k] %in% c("A", "T", "C", "G", "a", "t", "c", "g") == TRUE)
            wholegene[j,k] <- splitgene[k]
          if (splitgene[k] %in% c("A", "T", "C", "G", "a", "t", "c", "g") == FALSE)
            {wholegene[j,k] <- ref.matrix[1,k]
            invalid.site.counter <- invalid.site.counter + 1}
        }
      }#differentiate between valid and ambigious sites(count, and set as ref)
    #---------------------------------------------------------------------------------------------------
      cat (" Completing Step:", loop.count + 5, "of", (fasta.num * 7) + 16, "\r")
      for (j in 1:(seq.count-1)){
        n.change.matrix[j,] <- colSums(sweep(wholegene[(j):(seq.count),], 2, wholegene[(j),], FUN = `!=`,
                                             check.margin = FALSE) + 0)}
      n.change.matrix[seq.count,] <- ( wholegene[seq.count,] != wholegene[seq.count-1,] ) + 0
      n.change.matrix <- n.change.matrix*2
      snps <- colSums(n.change.matrix)
      seg.sites <- matrix()
      for (j in 1:seq.length){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}
      }
      seg.sites <- as.numeric(seg.sites[,-1])
      seg.sites.num <- length(seg.sites) #count number of seg sites i.e sites different from ref including SP, NS
      #set SPs, and NPs changes and collect the positions using AA table above--------------------------------------------------------------------------------------------------
      cat (" Completing Step:", loop.count + 6, "of", (fasta.num * 7) + 16, "\r")
      snp.count <- 0
      syn.count <- 0
      nsyn.count <- 0
      for (j in 1:seg.sites.num){
        snp.site <- seg.sites[j]
        if (seg.sites[j] %% 3 == 1){possible.AA.number <- length(
          unique(aalookup[match(unique(paste(wholegene[ ,snp.site],
                                             ref.matrix[ ,snp.site+1],
                                             ref.matrix[ ,snp.site+2], sep = "")), aalookup),2] ))} #remainder of seg sites divided by 3

        if (seg.sites[j] %% 3 == 2){possible.AA.number <- length(
          unique(aalookup[match(unique(paste(ref.matrix[ ,snp.site-1],
                                             wholegene[ ,snp.site],
                                             ref.matrix[ ,snp.site+1], sep = "")), aalookup),2] ))}

        if (seg.sites[j] %% 3 == 0){possible.AA.number <- length(
          unique(aalookup[match(unique(paste(ref.matrix[ ,snp.site-2],
                                             ref.matrix[ ,snp.site-1],
                                             wholegene[ ,snp.site], sep = "")), aalookup),2] ))}
        possible.nuc.number <- length(unique(wholegene[ ,seg.sites[j]]))
        snp.count  <- snp.count  + (possible.nuc.number - 1)
        nsyn.count <- nsyn.count + (possible.AA.number  - 1)
        syn.count  <- syn.count  + (possible.nuc.number - possible.AA.number)
      }
    #---------------------------------------------------------------------------------------------------
      cat (" Completing Step:", loop.count + 7, "of", (fasta.num * 7) + 16, "\r")
      hap.nucs <- unique(onepop[ ,1]) #unqiue string = hap
      hap.count <- length(hap.nucs) #number of hap
      hap.lookup <- cbind(hap.nucs, c(1:hap.count)) #give the number to each hap for hap ranking
      hap.table <- match(onepop[ ,1], hap.lookup) #gave each sample to respective assigned hap number
      hap.freq <- matrix(ncol = 2, nrow = hap.count) #built matrix according number of hap
      hap.freq[ ,1] <- c(1:hap.count)
      for (j in 1:hap.count){
        n.hap.counter <- 0
        for (k in 1:seq.count){
          if (hap.table[k] == j) {n.hap.counter <- n.hap.counter + 1}
        }
        hap.freq[j,2] <- n.hap.counter #put freq for each assigned hap
      }
     Hd <- (seq.count/(seq.count-1)) * (1-(sum(((hap.freq[,2])/seq.count)^2))) #hap diversity overall
    #---------------------------------------------------------------------------------------------------
      cat (" Completing Step:", loop.count + 8, "of", (fasta.num * 7) + 16, "\r")
      pi.per.nuc <- ((colSums(n.change.matrix) /  (seq.count)^2) / seq.length) * (seq.count/(seq.count-1)) #pi for each base
      tpi <- sum(colSums(n.change.matrix)) /  (seq.count) / seq.count  * (seq.count/(seq.count-1))
      pi <- sum(pi.per.nuc) #average pi (Nei method)
      #tajima's D equation_step by step calculation
      aone <- 0
      for (j in 1:(seq.count-1)){aone <- aone + 1/j}
      atwo <- 0
      for (j in 1:(seq.count-1)){atwo <- atwo + 1/(j^2)}
      bone <- (seq.count + 1) / (3*(seq.count - 1))
      btwo <- (2*(seq.count^2 + seq.count + 3)) / (9*seq.count*(seq.count-1))
      cone <- bone - (1/aone)
      ctwo <- btwo - ((seq.count+2)/(aone*seq.count)) + (atwo/(aone^2))
      eone <- cone / aone
      etwo <- ctwo / ((aone^2) + atwo)
      td <- ((tpi) -  (seg.sites.num / aone)) /
        (sqrt(  (eone*seg.sites.num)  + ( (etwo*seg.sites.num) * (seg.sites.num-1) ) ))
    # puting result to table---------------------------------------------------------------------------------------------------
      cat (" Completing Step:", loop.count + 9, "of", (fasta.num * 7) + 16, "\r")
      ndrt[i,1]  <- seq.count
      ndrt[i,2]  <- seq.length
      ndrt[i,3]  <- invalid.site.counter
      ndrt[i,4]  <- seg.sites.num
      ndrt[i,5]  <- snp.count
      ndrt[i,6]  <- round(pi*1000, 3)
      ndrt[i,7]  <- round(td, 3)
      ndrt[i,8]  <- nsyn.count
      ndrt[i,9]  <- syn.count
      ndrt[i,10]  <- hap.count
      ndrt[i,11] <- round(Hd, 3)}
    #---------------------------------------------------------------------------------------------------
    if (seq.count == 1){
      ndrt[i,1]  <- seq.count
      ndrt[i,2]  <- seq.length
      ndrt[i,3]  <- invalid.site.counter
      ndrt[i,4]  <- 0
      ndrt[i,5]  <- 0
      ndrt[i,6]  <- 0
      ndrt[i,7]  <- 0
      ndrt[i,8]  <- 0
      ndrt[i,9]  <- 0
      ndrt[i,10] <- 1
      ndrt[i,11] <- 0
    }
    totalpop <- rbind(totalpop, wholegene)
    loop.count <- loop.count + 7
    invalid.site.counter.total <- invalid.site.counter.total + invalid.site.counter
  }
  totalpop <- totalpop[-1,]
  seq.count <- dim(totalpop)[1]
  seq.length <- dim(totalpop)[2]
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 4, "of", (fasta.num * 7) + 16, "\r") #pi equation
  n.change.matrix <- matrix(nrow = seq.count, ncol = seq.length)
  for (j in 1:(seq.count-1)){
    n.change.matrix[j,] <- colSums(sweep(totalpop[(j):(seq.count),], 2, totalpop[(j),], FUN = `!=`,
                                         check.margin = FALSE) + 0)}
  n.change.matrix[seq.count,] <- ( totalpop[seq.count,] != totalpop[seq.count-1,] ) + 0
  n.change.matrix <- n.change.matrix*2
  snps <- colSums(n.change.matrix)
  seg.sites <- matrix()
  for (j in 1:seq.length){
    if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}
  } #getting position for segregation sites on nucleotide scales
  seg.sites <- as.numeric(seg.sites[,-1])
  seg.sites.num <- length(seg.sites)
  seg.codons <- unique(ceiling(seg.sites/3))
  seg.codons.num <- length(seg.codons)
  globalsnps <- snps
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 5, "of", (fasta.num * 7) + 16, "\r")
  snp.count <- 0
  syn.count <- 0
  nsyn.count <- 0
  for (j in 1:seg.sites.num){
    snp.site <- seg.sites[j]
    if (seg.sites[j] %% 3 == 1){possible.AA <- aalookup[match(paste(totalpop[ ,snp.site],
                                                                    ref.matrix[ ,snp.site+1],
                                                                    ref.matrix[ ,snp.site+2], sep = ""), aalookup),2]
    refAA <- aalookup[match(paste(ref.matrix[ ,snp.site],
                                  ref.matrix[ ,snp.site+1],
                                  ref.matrix[ ,snp.site+2], sep = ""), aalookup),2]

    possible.AA.number <- length(unique(possible.AA))
    if (refAA %in% possible.AA == FALSE) possible.AA.number <- possible.AA.number + 1}
    if (seg.sites[j] %% 3 == 2){possible.AA <- aalookup[match(paste(ref.matrix[ ,snp.site-1],
                                                                    totalpop[ ,snp.site],
                                                                    ref.matrix[ ,snp.site+1], sep = ""), aalookup),2]
    refAA <- aalookup[match(paste(ref.matrix[ ,snp.site-1],
                                  ref.matrix[ ,snp.site],
                                  ref.matrix[ ,snp.site+1], sep = ""), aalookup),2]
    possible.AA.number <- length(unique(possible.AA))
    if (refAA %in% possible.AA == FALSE) possible.AA.number <- possible.AA.number + 1}
    if (seg.sites[j] %% 3 == 0){possible.AA <- aalookup[match(paste(ref.matrix[ ,snp.site-2],
                                                                    ref.matrix[ ,snp.site-1],
                                                                    totalpop[ ,snp.site], sep = ""), aalookup),2]
    refAA <- aalookup[match(paste(ref.matrix[ ,snp.site-2],
                                  ref.matrix[ ,snp.site-1],
                                  ref.matrix[ ,snp.site], sep = ""), aalookup),2]
    possible.AA.number <- length(unique(possible.AA))
    if (refAA %in% possible.AA == FALSE) possible.AA.number <- possible.AA.number + 1}
    possible.nuc.number <- length(unique(totalpop[ ,seg.sites[j]]))
    snp.count  <- snp.count  + (possible.nuc.number - 1)
    nsyn.count <- nsyn.count + (possible.AA.number  - 1)
    syn.count  <- syn.count  + (possible.nuc.number - possible.AA.number)
  }
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 6, "of", (fasta.num * 7) + 16, "\r")
  tot.pop.string <- matrix(nrow = seq.count)
  for (i in 1:seq.count){tot.pop.string[i,1] <- paste(totalpop[i,], collapse = "")}
  hap.nucs <- unique(tot.pop.string[ ,1])
  hap.count <- length(hap.nucs) #hap based on nucleotide
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 7, "of", (fasta.num * 7) + 16, "\r")
  hap.lookup <- cbind(hap.nucs, c(1:hap.count))
  hap.table <- match(tot.pop.string[ ,1], hap.lookup)
  hap.freq <- matrix(ncol = 2, nrow = hap.count)
  hap.freq[ ,1] <- c(1:hap.count)
  for (j in 1:hap.count){
    n.hap.counter <- 0
    for (k in 1:seq.count){
      if (hap.table[k] == j) {n.hap.counter <- n.hap.counter + 1}}
    hap.freq[j,2] <- n.hap.counter}
  Hd <- (seq.count/(seq.count-1)) * (1-(sum(((hap.freq[,2])/seq.count)^2))) #total HD inside folder
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 8, "of", (fasta.num * 7) + 16, "\r")
  pi.per.nuc <- ((colSums(n.change.matrix) /  (seq.count)^2) / seq.length) * (seq.count/(seq.count-1))
  tpi <- sum(colSums(n.change.matrix)) /  (seq.count) / seq.count  * (seq.count/(seq.count-1))
  pi <- sum(pi.per.nuc)
  aone <- 0
  for (j in 1:(seq.count-1)){aone <- aone + 1/j}
  atwo <- 0
  for (j in 1:(seq.count-1)){atwo <- atwo + 1/(j^2)}
  bone <- (seq.count + 1) / (3*(seq.count - 1))
  btwo <- (2*(seq.count^2 + seq.count + 3)) / (9*seq.count*(seq.count-1))
  cone <- bone - (1/aone)
  ctwo <- btwo - ((seq.count+2)/(aone*seq.count)) + (atwo/(aone^2))
  eone <- cone / aone
  etwo <- ctwo / ((aone^2) + atwo)
  td <- ((tpi) -  (seg.sites.num / aone)) /
    (sqrt(  (eone*seg.sites.num)  + ( (etwo*seg.sites.num) * (seg.sites.num-1) ) )) #TD for total i.e all fasta files inside folder
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 9, "of", (fasta.num * 7) + 16, "\r") #filling up total row from result table option 1
  ndrt[fasta.num + 1, 1]  <- seq.count
  ndrt[fasta.num + 1, 2]  <- seq.length
  ndrt[fasta.num + 1, 3]  <- invalid.site.counter.total
  ndrt[fasta.num + 1, 4]  <- seg.sites.num
  ndrt[fasta.num + 1, 5]  <- snp.count
  ndrt[fasta.num + 1, 6]  <- round(pi*1000, 3) #muptiply by 1000 so we can show the value of 10^-3
  ndrt[fasta.num + 1, 7]  <- round(td, 3)
  ndrt[fasta.num + 1, 8]  <- nsyn.count
  ndrt[fasta.num + 1, 9]  <- syn.count
  ndrt[fasta.num + 1, 10] <- hap.count
  ndrt[fasta.num + 1, 11] <- round(Hd, 3)
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 10, "of", (fasta.num * 7) + 16, "\r")
  AAs.each.variable.codon <- matrix(nrow = seq.count, ncol = seg.codons.num)
  colnames(AAs.each.variable.codon) <- seg.codons #seg.condons = sites of segregation in AA
  ref.each.variable.codon <- matrix(nrow = 1, ncol = seg.codons.num)
  colnames(ref.each.variable.codon) <- seg.codons #for reference seq
  max.AA.variants <- 0
  for (i in 1:seg.codons.num){
    codon.num <- seg.codons[i]
    refAA <- aalookup[match(paste(ref.matrix[ ,(codon.num*3)-2],
                                  ref.matrix[ ,(codon.num*3)-1],
                                  ref.matrix[ ,(codon.num*3)], sep = ""), aalookup),2]
    ref.each.variable.codon[ ,i] <- refAA #look at reference sequence, and check
    possible.AA <- aalookup[match(paste(totalpop[ ,(codon.num*3)-2],
                                        totalpop[ ,(codon.num*3)-1],
                                        totalpop[ ,(codon.num*3)], sep = ""), aalookup),2]
    AAs.each.variable.codon[ ,i] <- possible.AA #extract segregation sites from input fasta and display
    possible.AA.number <- length(unique(c(refAA, possible.AA)))
    if (possible.AA.number > max.AA.variants) max.AA.variants <- possible.AA.number}
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 11, "of", (fasta.num * 7) + 16, "\r")
  did.AA.change <- colSums(sweep(AAs.each.variable.codon, 2, ref.each.variable.codon, FUN = `!=`) + 0) / seq.count #non-ref allele freq
  aa.variant.table <- matrix(nrow = max.AA.variants + 1, ncol = seg.codons.num)
  colnames(aa.variant.table) <- seg.codons
  aa.variant.table[1, ] <- did.AA.change #change of amino acid
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 12, "of", (fasta.num * 7) + 16, "\r")
  maf.AA.table <- matrix(nrow = 22, ncol = dim(aa.variant.table)[2]) #making matrix for AA - table option 3
  aalist <- c("I", "L", "V", "F", "M", "C", "A", "G", "P", "T", "S", "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "STOP")
  rownames(maf.AA.table) <- c(aalist, "REF")
  colnames(maf.AA.table) <- seg.codons #empty ready to put matrix
  #---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 13, "of", (fasta.num * 7) + 16, "\r")
  seg.codons.ref.variants.matrix <- rbind(ref.each.variable.codon, AAs.each.variable.codon) #puting ref and changes in fasta input into table
  for (i in 1:seg.codons.num){
    for (j in 1:max.AA.variants){
      possible.AA <- unique(seg.codons.ref.variants.matrix[ ,i])
      length(possible.AA) <- max.AA.variants
      aa.variant.table[j+1,i] <- sum(sweep(as.matrix(AAs.each.variable.codon[ ,i]), 2, possible.AA[j], FUN = `==`) + 0)
      maf.AA.table[match(possible.AA[j], aalist),i] <- round(((aa.variant.table[j+1,i]) * 100 / seq.count), 3)
    }
    maf.AA.table[22,i] <- possible.AA[1] #put respective AA changes into the table
  }
  maf.AA.table[is.na(maf.AA.table)] <- 0 #replace NA with zero, and the table is done
  # the following is to make the same as above but for nucleotide---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 14, "of", (fasta.num * 7) + 16, "\r")
  maf.nuc.table <- matrix(nrow = 5, ncol = seg.sites.num)
  possible.nucs <- c("A", "T", "G", "C")
  for (i in 1:seg.sites.num){
    for (j in 1:4){
      nuc.counter <- sum(sweep(
        as.matrix(totalpop[ ,seg.sites[i]]), 2, possible.nucs[j], FUN = `==`) + 0)
      maf.nuc.table[j,i] <- round(((nuc.counter * 100) / seq.count), 3)
    }
    maf.nuc.table[5,i] <- ref.matrix[seg.sites[i]]
  }
  maf.nuc.table[is.na(maf.nuc.table)] <- 0
  colnames(maf.nuc.table) <- seg.sites
  rownames(maf.nuc.table) <- c(possible.nucs, "REF")

#the following is for adding TD simulation values---------------------------------------------------------------------------------------------------
  cat (" Completing Step:", (fasta.num * 7) + 15, "of", (fasta.num * 7) + 16, "\r")
# building matrix in R
td.sim <- matrix(ncol = 9, nrow = 73)

td.sim[ ,1] <- c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
                 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                 49, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 175, 200, 250,
                 300, 350, 400, 450, 500, 600, 800, 1000)
td.sim[ ,2] <- c(-0.876, -1.255, -1.405, -1.498, -1.522, -1.553, -1.559, -1.572, -1.573, -1.58, -1.58,
                 -1.584, -1.583, -1.585, -1.584, -1.585, -1.584, -1.585, -1.584, -1.584, -1.583, -1.583,
                 -1.582, -1.582, -1.581, -1.581, -1.58, -1.58, -1.579, -1.579, -1.578, -1.578, -1.577,
                 -1.577, -1.576, -1.576, -1.575, -1.575, -1.574, -1.574, -1.573, -1.573, -1.572, -1.572,
                 -1.571, -1.571, -1.57, -1.568, -1.566, -1.565, -1.563, -1.561, -1.56, -1.559, -1.557,
                 -1.556, -1.555, -1.552, -1.55, -1.549, -1.547, -1.545, -1.542, -1.539, -1.534, -1.53,
                 -1.526, -1.523, -1.521, -1.519, -1.515, -1.51, -1.505)
td.sim[ ,3] <- c(2.081, 1.737, 1.786, 1.728, 1.736, 1.715, 1.719, 1.71, 1.713, 1.708, 1.71, 1.708, 1.709,
                 1.708, 1.709, 1.708, 1.71, 1.709, 1.711, 1.71, 1.712, 1.712, 1.712, 1.712, 1.713, 1.714,
                 1.714, 1.714, 1.715, 1.716, 1.716, 1.717, 1.717, 1.717, 1.718, 1.718, 1.719, 1.719, 1.72,
                 1.72, 1.721, 1.721, 1.721, 1.722, 1.722, 1.722, 1.723, 1.724, 1.726, 1.727, 1.729, 1.73,
                 1.731, 1.732, 1.733, 1.734, 1.735, 1.737, 1.739, 1.74, 1.741, 1.743, 1.746, 1.748, 1.752,
                 1.755, 1.757, 1.759, 1.761, 1.763, 1.765, 1.769, 1.772)
td.sim[ ,4] <- c(-0.876, -1.269, -1.478, -1.608, -1.663, -1.713, -1.733, -1.757, -1.765, -1.779, -1.783,
                 -1.791, -1.793, -1.798, -1.799, -1.802, -1.803, -1.805, -1.804, -1.806, -1.806, -1.807,
                 -1.807, -1.807, -1.807, -1.807, -1.807, -1.807, -1.806, -1.806, -1.806, -1.806, -1.805,
                 -1.805, -1.804, -1.804, -1.804, -1.803, -1.803, -1.803, -1.802, -1.802, -1.801, -1.801,
                 -1.8, -1.8, -1.8, -1.797, -1.795, -1.793, -1.791, -1.79, -1.788, -1.786, -1.784, -1.783,
                 -1.781, -1.779, -1.776, -1.774, -1.771, -1.769, -1.765, -1.76,-1.754, -1.748, -1.744,
                 -1.74, -1.737, -1.734, -1.728, -1.721, -1.715)
td.sim[ ,5] <- c(2.232, 1.834, 1.999, 1.932, 1.975, 1.954, 1.975, 1.966, 1.979, 1.976, 1.985, 1.984, 1.99,
                 1.99, 1.996, 1.996, 2.001, 2.001, 2.005, 2.006, 2.009, 2.01, 2.013, 2.014, 2.017, 2.018,
                 2.02, 2.021, 2.023, 2.024, 2.026, 2.027, 2.029, 2.03, 2.031, 2.032, 2.033, 2.034, 2.036,
                 2.037, 2.038, 2.039, 2.04, 2.041, 2.042, 2.042, 2.044, 2.048, 2.052, 2.055, 2.058, 2.061,
                 2.064, 2.066, 2.069, 2.071, 2.073, 2.077, 2.08, 2.084, 2.086, 2.089, 2.095, 2.1, 2.107,
                 2.114, 2.119, 2.123, 2.127, 2.13, 2.135, 2.143, 2.15)
td.sim[ ,6] <- c(-0.876, -1.275, -1.54, -1.721, -1.83, -1.916, -1.967, -2.014, -2.041, -2.069, -2.085, -
                   2.103, -2.113, -2.126, -2.132, -2.141, -2.146, -2.152, -2.153, -2.16, -2.162, -2.165,
                 -2.167, -2.17, -2.171, -2.173, -2.173, -2.175, -2.175, -2.177, -2.177, -2.178, -2.178,
                 -2.179, -2.178, -2.179, -2.179, -2.179, -2.179, -2.179, -2.179, -2.179, -2.179, -2.179,
                 -2.178, -2.178, -2.178, -2.177, -2.175, -2.173, -2.171, -2.17, -2.168, -2.166, -2.164,
                 -2.162, -2.16, -2.157, -2.153, -2.15, -2.147, -2.144, -2.138, -2.132, -2.122, -2.114,
                 -2.107, -2.101, -2.096, -2.092, -2.084, -2.072, -2.062)
td.sim[ ,7] <- c(2.324, 1.901, 2.255, 2.185, 2.313, 2.296, 2.362, 2.359, 2.401, 2.403, 2.432, 2.436,
                 2.457, 2.461, 2.478, 2.483, 2.496, 2.501, 2.512, 2.516, 2.526, 2.53, 2.538, 2.542, 2.549,
                 2.553, 2.559, 2.563, 2.569, 2.572, 2.577, 2.58, 2.585, 2.588, 2.592, 2.595, 2.599, 2.601,
                 2.605, 2.608, 2.611, 2.613, 2.617, 2.619, 2.622, 2.624, 2.627, 2.638, 2.649, 2.658, 2.666,
                 2.673, 2.681, 2.687, 2.693, 2.699, 2.704, 2.713, 2.722, 2.73, 2.736, 2.743, 2.757, 2.768,
                 2.787, 2.802, 2.814, 2.824, 2.833, 2.84, 2.853, 2.873, 2.887)
td.sim[ ,8] <- c(-0.876, -1.276, -1.556, -1.761, -1.909, -2.023, -2.105, -2.174, -2.223, -2.267, -2.299,
                 -2.329, -2.35, -2.372, -2.387, -2.403, -2.414, -2.426, -2.434, -2.443, -2.449, -2.457,
                 -2.461, -2.467, -2.471, -2.475, -2.478, -2.482, -2.484, -2.487, -2.489, -2.492, -2.493,
                 -2.495, -2.496, -2.498, -2.499, -2.5, -2.501, -2.502, -2.502, -2.503, -2.504, -2.504,
                 -2.505, -2.505, -2.505, -2.506, -2.506, -2.506, -2.505, -2.504, -2.502, -2.5, -2.499,
                 -2.497, -2.495, -2.492, -2.488, -2.484, -2.481, -2.477, -2.47, -2.462, -2.449, -2.439,
                 -2.43, -2.422, -2.415, -2.409, -2.398, -2.382, -2.369)
td.sim[ ,9] <- c(2.336, 1.913, 2.373, 2.311, 2.524, 2.519, 2.64, 2.649, 2.729, 2.741, 2.798, 2.811,
                 2.854, 2.866, 2.9, 2.911, 2.939, 2.95, 2.973, 2.983, 3.002, 3.011, 3.029, 3.037, 3.052,
                 3.06, 3.073, 3.08, 3.092, 3.099, 3.11, 3.116, 3.126, 3.132, 3.141, 3.147, 3.155, 3.16,
                 3.168, 3.173, 3.18, 3.185, 3.191, 3.196, 3.202, 3.207, 3.212, 3.235, 3.256, 3.274, 3.291,
                 3.306, 3.32, 3.333, 3.345, 3.355, 3.366, 3.385, 3.401, 3.416, 3.43, 3.443, 3.47, 3.492,
                 3.529, 3.558, 3.581, 3.6, 3.617, 3.632, 3.657, 3.694, 3.722)
colnames(td.sim) <- c("n", 0.1, 0.1, 0.05, 0.05, 0.01, 0.01, 0.001, 0.001)

# assiging to sensible name and unassinged variables will be removed later on---------------------------------------------------------------------------------------------------
cat (" Completing Step:", (fasta.num * 7) + 16, "of", (fasta.num * 7) + 16, "\r")
  vp.AA.VARIANT.TABLE <- aa.variant.table
  vp.AAs.EACH.VARIABLE.CODON <- AAs.each.variable.codon
  vp.REF.EACH.VARIABLE.CODON <- ref.each.variable.codon
  vp.BASIC.RESULTS <- ndrt
  vp.MAF.AA.TABLE <- maf.AA.table
  vp.MAF.NUC.TABLE <- maf.nuc.table
  vp.GLOBAL.SNPS <- globalsnps
  vp.SEQ.NUM <- seq.count #total samples in fasta files
  vp.SEQ.LENGTH <- seq.length
  vp.GENE.NAME <- genename
  vp.SEG.CODONS <- seg.codons
  vp.POP.NUM <- fasta.num
  vp.MAX.AA.VARIANTS <- max.AA.variants
  vp.A1 <- aone
  vp.E1 <- eone
  vp.E2 <- etwo
  vp.TD.SIM <- td.sim
  vp.ALL.SEQUENCES.MATRIX <- totalpop
  vp.REF.MATRIX <- ref.matrix
  vp.SEG.SITES <- seg.sites

  vp.data <- list()
  vp.data <<- list(vp.AA.VARIANT.TABLE, vp.AAs.EACH.VARIABLE.CODON, vp.REF.EACH.VARIABLE.CODON,
                   vp.BASIC.RESULTS, vp.MAF.NUC.TABLE, vp.MAF.AA.TABLE, vp.GLOBAL.SNPS, vp.SEQ.NUM, vp.SEQ.LENGTH,
                   vp.GENE.NAME, vp.SEG.CODONS, vp.POP.NUM, vp.MAX.AA.VARIANTS, vp.A1, vp.E1, vp.E2, vp.TD.SIM,
                   vp.ALL.SEQUENCES.MATRIX, vp.REF.MATRIX, vp.SEG.SITES)

  cat ("Completed! Step:", (fasta.num * 7) + 16, "of", (fasta.num * 7) + 16, "\r")
  cat ("\n")
  et <- Sys.time()
  cat("vaxpack_input() took", round(et - st, digits = 2), units.difftime(et - st), "\n")
  cat ("Now use vaxpack_output() to get your results! \n")
}
