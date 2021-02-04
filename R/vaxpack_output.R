#two different huge functions have written, and results are extracted by number which acts like mini functions
#vaxpact_input mainly does the job of accepting input and core calculation.
#vaxpack_output does extract the results, and carry out graphical output.
vaxpack_output <- function () {

  if (exists("vp.data") != TRUE)
    stop("Please run the input function 'vaxpack_input()' first")
#assign everything needed into vp.data
  vp.data <- vp.data
  vp.AA.VARIANT.TABLE <- vp.data[[1]]
  vp.AAs.EACH.VARIABLE.CODON <- vp.data[[2]]
  vp.REF.EACH.VARIABLE.CODON <- vp.data[[3]]
  vp.BASIC.RESULTS <- vp.data[[4]] #table 1
  vp.MAF.NUC.TABLE <- vp.data[[5]]
  vp.MAF.AA.TABLE <- vp.data[[6]]
  vp.GLOBAL.SNPS <- vp.data[[7]]
  vp.SEQ.NUM <- vp.data[[8]]
  vp.SEQ.LENGTH <- vp.data[[9]]
  vp.GENE.NAME <- vp.data[[10]]
  vp.SEG.CODONS <- vp.data[[11]]
  vp.POP.NUM <- vp.data[[12]]
  vp.MAX.AA.VARIANTS <- vp.data[[13]]
  vp.A1 <- vp.data[[14]]
  vp.E1 <- vp.data[[15]]
  vp.E2 <- vp.data[[16]]
  vp.TD.SIM <- vp.data[[17]]
  vp.ALL.SEQUENCES.MATRIX <- vp.data[[18]]
  vp.REF.MATRIX <- vp.data[[19]]
  vp.SEG.SITES <- vp.data[[20]]
  vp.total.row.num <- nrow(vp.BASIC.RESULTS)

  cat("\014")
  cat ("Choose your output:\n")
  cat ("TABLES\n")
  cat ("1  - Results Table \n")
  cat ("2  - Haplotype Table\n")
  cat ("3  - Minor Allele Frequency Table (Nucleotides)\n")
  cat ("4  - Minor Allele Frequency Table (Amino Acids)\n")
  cat ("\n")
  cat ("GRAPHS\n")
  cat ("5  - Haplotype Population Pie Chart\n")
  cat ("6  - AA Variant Percentage Column Graph\n")
  cat ("7  - Phylogenetic Tree\n")
  cat ("8  - Haplotype Accumulation Plot\n")
  cat ("\n")
  cat ("SLIDING SCALE\n")
  cat ("9  - Polymorphism\n")
  cat ("10 - Tajima's D\n")
  cat ("11 - Nucleotide Diversity\n")
  cat ("12 - All Sliding Scale Graphs Overlapped \n")
  cat ("\n")
  cat ("13 - Sliding Scale data table \n")
  cat ("\n")

  choice <- as.numeric( readline ("Enter a number from the selection above -  ")) #readline does read the user defined number
  if (choice %in% c(1:13) == FALSE) stop("There is no output option for that number")
  cat ("\n")

  #RESULTS TABLE-------------------------------------------------------------------------
  if (choice == 1) {
    cat("Minimal haplotypes will be calculated using only the segregation sites where polymorphism\n")
    cat("is found in at least 'x' percent of the population at that site")
    threshold <- as.numeric(readline ("'x' <-  "))

    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,]<=threshold/100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    colnames(vp.AAs.EACH.VARIABLE.CODON) <- vp.SEG.CODONS #segregation sites AA
    AAs.each.variable.codon.over.threshold <- vp.AAs.EACH.VARIABLE.CODON
    AAs.each.variable.codon.over.threshold <- AAs.each.variable.codon.over.threshold[ ,
      colnames(AAs.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)] #this return amino acid sites that have MAF over user specified values
    min.hap.strings <- matrix(nrow = vp.SEQ.NUM)
    #AA_based Hap calculation
    for (i in 1:vp.SEQ.NUM){
      min.hap.strings[i,1] <- paste(AAs.each.variable.codon.over.threshold[i,], collapse = "")}
    min.hap.count.col <- matrix(nrow = vp.POP.NUM + 1)
    min.hap.div.col <- matrix(nrow = vp.POP.NUM + 1)
    counter <- 0
    for (i in 1:vp.POP.NUM){
      min.hap.count.col[i,1] <- length(
        unique(min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), 1]))
      min.hap.AAs.col <- matrix(nrow = min.hap.count.col[i,1])
      #the following code will give you unique hap above specified threshold
      min.hap.AAs.col[ ,1] <- unique(
        min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), ])
      min.hap.table <- match(min.hap.strings[(1+counter):((vp.BASIC.RESULTS[i,1]) + counter), ],
                             min.hap.AAs.col)
      min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count.col[i,1]) #Hap based on AA and its frequency
      min.hap.freq[ ,1] <- c(1:min.hap.count.col[i,1])
      for (j in 1:min.hap.count.col[i,1]){
        AA.hap.counter <- 0
        for (k in 1:vp.BASIC.RESULTS[i,1]){
          if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
        min.hap.freq[j,2] <- AA.hap.counter}
      min.hap.div.col[i,1] <-
        (vp.BASIC.RESULTS[i,1]/
           (vp.BASIC.RESULTS[i,1]-1)) * (1-(sum(((min.hap.freq[,2])/vp.BASIC.RESULTS[i,1])^2))) #hap diversity statistic
      counter <- counter + (vp.BASIC.RESULTS[i,1])}
#for total as above for result table option 1
      min.hap.count.col[vp.total.row.num,1] <- length(unique(min.hap.strings[,1]))
      min.hap.AAs.col <- matrix(nrow = min.hap.count.col[vp.total.row.num,1])
      min.hap.AAs.col[ ,1] <- unique(min.hap.strings[ ,1])
      min.hap.table <- match(min.hap.strings[ ,1], min.hap.AAs.col)
      min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count.col[vp.total.row.num,1])
      min.hap.freq[ ,1] <- c(1:min.hap.count.col[vp.total.row.num,1])
      for (j in 1:min.hap.count.col[vp.total.row.num,1]){
        AA.hap.counter <- 0
        for (k in 1:vp.BASIC.RESULTS[vp.total.row.num,1]){
          if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
        min.hap.freq[j,2] <- AA.hap.counter}
      min.hap.div.col[vp.total.row.num,1] <-
        (vp.BASIC.RESULTS[vp.total.row.num,1]/
           (vp.BASIC.RESULTS[vp.total.row.num,1]-1)) * (1-(sum(((min.hap.freq[,2])/
                vp.BASIC.RESULTS[vp.total.row.num,1])^2)))

    vp.RESULTS.TABLE <<- cbind(vp.BASIC.RESULTS,
                               round(min.hap.count.col, 3), round(min.hap.div.col, 3))
    colnames(vp.RESULTS.TABLE)[12:13] <- c("Minimal AA h", "Minimal AA Hd")
    cat("\n")
    cat("Saved as \"vp.RESULTS.TABLE\", use write.csv() to save to excel \n")
    cat("e.g. write.csv(vp.RESULTS.TABLE, file = \"my.results.table.in.excel\")")
    View(vp.RESULTS.TABLE)
    }


  #HAPLOTYPE RESULTS TABLE-----------------------------------------------------------------------------
  if (choice == 2) {
    cat("Haplotypes will be calculated using only the segregation sites where polymorphism\n")
    cat("is found in at least 'x' percent of the population at that site")
    threshold <- as.numeric(readline ("'x' <-  "))
    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,]<=threshold/100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    AAs.each.variable.codon.over.threshold <- vp.AAs.EACH.VARIABLE.CODON
    AAs.each.variable.codon.over.threshold <- AAs.each.variable.codon.over.threshold[ ,
      colnames(AAs.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    ref.each.variable.codon.over.threshold <- vp.REF.EACH.VARIABLE.CODON
    ref.each.variable.codon.over.threshold <- ref.each.variable.codon.over.threshold[ ,
      colnames(ref.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    AA.changes.per.hap <- rowSums(sweep(unique(AAs.each.variable.codon.over.threshold), 2,
                                        ref.each.variable.codon.over.threshold, FUN = `!=`,
                                        check.margin = FALSE) +0)
    min.hap.strings <- matrix(nrow = vp.SEQ.NUM)
    for (i in 1:vp.SEQ.NUM){
      min.hap.strings[i,1] <- paste(AAs.each.variable.codon.over.threshold[i,], collapse = "")}
    min.hap.count <- length(AA.changes.per.hap)
    min.hap.AAs <- unique(min.hap.strings)
    min.hap.table <- match(min.hap.strings, min.hap.AAs)
    min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count)
    min.hap.freq[ ,1] <- c(1:min.hap.count)
    for (j in 1:min.hap.count){
      AA.hap.counter <- 0
      for (k in 1:vp.SEQ.NUM){
        if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
      min.hap.freq[j,2] <- AA.hap.counter}
    min.hap.freq <- cbind(min.hap.freq, AA.changes.per.hap)
    min.hap.freq <- as.data.frame(min.hap.freq)
    min.hap.freq <- min.hap.freq[order(min.hap.freq$V2, decreasing = TRUE),]
    min.hap.freq <- cbind(min.hap.freq[ ,2],
                          round((min.hap.freq[ ,2]/vp.SEQ.NUM)*100, 3), min.hap.freq[ ,3])
    rownames(min.hap.freq) <- c(1:min.hap.count)
    colnames(min.hap.freq) <- c("Frequency", "Proportion %", "AA Difference vs Reference")

    vp.HAPLOTYPE.TABLE <<- min.hap.freq
    cat("\n")
    cat("Saved as \"vp.HAPLOTYPE.TABLE\", use write.csv() to save to excel \n")
    cat("e.g. write.csv(vp.HAPLOTYPE.TABLE, file = \"my.haplotype.table.in.excel\")")
    suppressWarnings(View(vp.HAPLOTYPE.TABLE))
    }

  #MAF NUCs----------------------------------------------------------------------------------
  if (choice == 3) {
    vp.MAF.NUC.TABLE <<- vp.MAF.NUC.TABLE
    cat("\n")
    cat("Saved as \"vp.MAF.NUC.TABLE\", use write.csv() to save to excel \n")
    cat("e.g. write.csv(vp.MAF.NUC.TABLE, file = \"my.MAF.NUC.TABLE.in.excel\")")
    suppressWarnings(View(vp.MAF.NUC.TABLE))}

  #MAF AA---------------------------------------------------------------
  if (choice == 4) {
    vp.MAF.AA.TABLE <<- vp.MAF.AA.TABLE
    cat("\n")
    cat("Saved as \"vp.MAF.AA.TABLE\", use write.csv() to save to excel \n")
    cat("e.g. write.csv(vp.MAF.AA.TABLE, file = \"my.MAF.AA.table.in.excel\")")
    suppressWarnings(View(vp.MAF.AA.TABLE))}


  #HAPLOTYPE PIE CHART--------------------------------------------------------------------
  if (choice == 5) {

    cat("Haplotypes will be calculated using only the segregation sites where polymorphism\n")
    cat("is found in at least 'x' percent of the population at that site")
    threshold <- as.numeric(readline ("'x' <-  "))
    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,]<=threshold/100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    AAs.each.variable.codon.over.threshold <- vp.AAs.EACH.VARIABLE.CODON
    AAs.each.variable.codon.over.threshold <- AAs.each.variable.codon.over.threshold[ ,
      colnames(AAs.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    ref.each.variable.codon.over.threshold <- vp.REF.EACH.VARIABLE.CODON
    ref.each.variable.codon.over.threshold <- ref.each.variable.codon.over.threshold[ ,
      colnames(ref.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    AA.changes.per.hap <- rowSums(sweep(unique(AAs.each.variable.codon.over.threshold), 2,
                                        ref.each.variable.codon.over.threshold, FUN = `!=`,
                                        check.margin = FALSE) +0)
    min.hap.strings <- matrix(nrow = vp.SEQ.NUM)
    for (i in 1:vp.SEQ.NUM){
      min.hap.strings[i,1] <- paste(AAs.each.variable.codon.over.threshold[i,], collapse = "")}
    min.hap.count <- length(AA.changes.per.hap)
    min.hap.AAs <- unique(min.hap.strings)
    min.hap.table <- match(min.hap.strings, min.hap.AAs)
    min.hap.freq <- matrix(ncol = 2, nrow = min.hap.count)
    min.hap.freq[ ,1] <- c(1:min.hap.count)
    for (j in 1:min.hap.count){
      AA.hap.counter <- 0
      for (k in 1:vp.SEQ.NUM){
        if (min.hap.table[k] == j) {AA.hap.counter <- AA.hap.counter + 1}}
      min.hap.freq[j,2] <- AA.hap.counter}
    min.hap.freq <- cbind(min.hap.freq, AA.changes.per.hap)
    min.hap.freq <- as.data.frame(min.hap.freq)
    min.hap.freq <- min.hap.freq[order(min.hap.freq$V2, decreasing = TRUE),]
    min.hap.freq <- cbind(min.hap.freq[ ,2], round((min.hap.freq[ ,2]/vp.SEQ.NUM)*100, 3),
                          min.hap.freq[ ,3])
    rownames(min.hap.freq) <- c(1:min.hap.count)
    colnames(min.hap.freq) <- c("Frequency", "Proportion %", "AA Difference vs Reference")

    vp.HAPLOTYPE.TABLE <- as.data.frame(min.hap.freq)
    vp.HAPLOTYPE.TABLE <- cbind(vp.HAPLOTYPE.TABLE, c(1:min.hap.count))

    vp.Haplotype.Pie.Chart <<-plot_ly(vp.HAPLOTYPE.TABLE, labels = c(1:min.hap.count),
                      values = vp.HAPLOTYPE.TABLE$Frequency,
                      type = 'pie', showlegend = FALSE, textposition = FALSE,
                      hovertext = paste(vp.HAPLOTYPE.TABLE$`AA Difference vs Reference`,
                                        "AA Differences vs Reference"),
                      hoverinfo = 'text') %>%
      layout(title = paste('Haplotype Distribution of', vp.GENE.NAME),
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    cat("\n")
    cat("Saved as \"vp.Haplotype.Pie.Chart\"\n")
    cat("Hover over pie chart to see the number of AA differences from the reference\n")
    return(vp.Haplotype.Pie.Chart)
  }

  #AA VARIANT GRAPH---------------------------------------------------------------------------
  if (choice == 6) {
    cat("Only segregation sites which have a non-reference variation\n")
    cat("in at least 'x' percent of the population will be graphed\n")
    threshold <- as.numeric(readline ("'x' <-  "))
    vp.AA.VARIANT.TABLE <- as.data.frame(vp.AA.VARIANT.TABLE)
    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,] <= threshold / 100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    minimal.variant.table <- rbind(colnames(variant.codons.over.threshold),
                                   variant.codons.over.threshold)
    var.table.row.names <- c("SITE", "Non-ref Freq.", "REF", "ALT 1",
                             "ALT 2", "ALT 3", "ALT 4", "ALT 5", "ALT 6",
                             "ALT 7", "ALT 8", "ALT 9", "ALT 10", "ALT 11",
                             "ALT 12", "ALT 13", "ALT 14","ALT 15", "ALT 16",
                             "ALT 17", "ALT 18", "ALT 19", "ALT 20")
    length(var.table.row.names) <- vp.MAX.AA.VARIANTS + 2
    rownames(minimal.variant.table) <- var.table.row.names
    minimal.variant.table <- minimal.variant.table[-2, ]
    min.var.plot.data <- matrix(nrow = vp.MAX.AA.VARIANTS*dim(minimal.variant.table)[2], ncol = 3)
    min.var.plot.data[ ,1] <- as.numeric(minimal.variant.table[1, ])
    for (i in 1:vp.MAX.AA.VARIANTS){
      min.var.plot.data[((dim(minimal.variant.table)[2])*(i-1)+1):((dim(minimal.variant.table)[2])*(i)),
                        2] <- as.matrix(rownames(minimal.variant.table)[i+1])}
    value.column <- matrix(nrow = 1, ncol = 1)
    for (i in 1:vp.MAX.AA.VARIANTS){value.column <- as.matrix(cbind(value.column,
                                                                    minimal.variant.table[(i+1),]))}
    min.var.plot.data <- as.data.frame(min.var.plot.data)
    min.var.plot.data[ ,3] <- as.numeric(as.character(value.column[ ,-1]))
    min.var.plot.data[is.na(min.var.plot.data)] <- 0
    colnames(min.var.plot.data) <- c("codons", "variable", "value")
    min.var.plot.data$codons = factor(min.var.plot.data$codons, levels = unique(min.var.plot.data$codons))
    vp.AA.Variant.Graph <<- ggplot(data = min.var.plot.data,
                                      aes(x = min.var.plot.data$codons,
                                          y = min.var.plot.data$value)) +
      geom_bar(aes( fill = min.var.plot.data$variable), stat = 'identity')+
      scale_fill_brewer(palette="Set1",
                        name = "Amino Acid")+
      labs(subtitle = vp.GENE.NAME, y = "Percentage (%)", x = "Codon",
           title = "Amino Acid Variation")+
      scale_y_discrete(limits = c((vp.SEQ.NUM/10), (vp.SEQ.NUM/5), (3*vp.SEQ.NUM/10),
                                  (4*vp.SEQ.NUM/10), (vp.SEQ.NUM/2), (6*vp.SEQ.NUM/10),
                                  (7*vp.SEQ.NUM/10), (8*vp.SEQ.NUM/10), (9*vp.SEQ.NUM/10), vp.SEQ.NUM),
                       labels = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))#, axis.text.x = element_text(angle = 45))
    cat("\n")
    cat("Saved as \"vp.AA.Variant.Graph\"\n")
    return(vp.AA.Variant.Graph)}

  #PHYLOGENETIC TREE-----------------------------------------------------------------------------
  if (choice == 7)
  #Phylogenic tree
    {
    cat("A neighbour joining tree will be drawn, only considering segregation sites\n")
    cat("where a non-reference variation is found in at least 'x' percent of the\n")
    cat("population at that site\n")
    threshold <- as.numeric(readline ("'x' <-  "))

    snps.to.include <- matrix()
    for (i in 1:length(vp.SEG.SITES)){
      if(vp.MAF.NUC.TABLE[5,i] == "A")
        if((100 - as.numeric(vp.MAF.NUC.TABLE[1,i])) > threshold)
          snps.to.include <- cbind(snps.to.include, vp.SEG.SITES[i])
      if(vp.MAF.NUC.TABLE[5,i] == "T")
        if((100 - as.numeric(vp.MAF.NUC.TABLE[2,i])) > threshold)
          snps.to.include <- cbind(snps.to.include, vp.SEG.SITES[i])
      if(vp.MAF.NUC.TABLE[5,i] == "G")
        if((100 - as.numeric(vp.MAF.NUC.TABLE[3,i])) > threshold)
          snps.to.include <- cbind(snps.to.include, vp.SEG.SITES[i])
      if(vp.MAF.NUC.TABLE[5,i] == "C")
        if((100 - as.numeric(vp.MAF.NUC.TABLE[4,i])) > threshold)
          snps.to.include <- cbind(snps.to.include, vp.SEG.SITES[i])
    }
    snps.to.include <- snps.to.include[,-1]

    all.to.bin <- matrix(nrow = vp.SEQ.NUM, ncol = vp.SEQ.LENGTH)
    for (i in 1:vp.SEQ.NUM){
      all.to.bin[(i),] <- vp.REF.MATRIX
    }

    for (i in 1:vp.SEQ.NUM){
      for (j in 1:length(snps.to.include)){
        all.to.bin[(i),snps.to.include[j]] <- vp.ALL.SEQUENCES.MATRIX[i,snps.to.include[j]]
      }
    }

    genebin <- as.DNAbin(all.to.bin)

    gene.distance <- dist.gene(genebin,pairwise.deletion = T, method = "percentage")
    gene.distance <- as.matrix(gene.distance)
    tree.by.population <- nj(gene.distance)

    by.population.labels <- matrix(nrow = vp.SEQ.NUM)
    counter <- 0
    for (i in 1:vp.POP.NUM){
      for (j in 1:vp.BASIC.RESULTS[i,1]){
        by.population.labels[(j + counter),1] <- i
      }
      counter <- counter + vp.BASIC.RESULTS[i,1]
    }

    tree.colours <- c('firebrick1', 'dodgerblue', 'green', 'orange', 'gold',
                      'darkviolet', 'black', 'deeppink', 'lightsalmon', 'lightseagreen',
                      'midnightblue', 'magenta', 'rosybrown', 'plum', 'snow4', 'tomato',
                      'slateblue2', 'snow3', 'khaki1', 'darkgreen', 'cyan')
    other.colours <- colors()[colors() %in% tree.colours == FALSE]
    tree.colours <- c(tree.colours, other.colours)

    by.population.labels <- as.numeric(by.population.labels)
    by.population.labels <- tree.colours[by.population.labels]

    phylo_tree_labels <- (c(rep("", vp.SEQ.NUM)))

    vp.Phylogenetic.Tree <<- function(){
      plot.phylo(tree.by.population, "fan", lab4ut="axial", show.tip.label = F,
                 cex=0.5, pch=c(15), rotate.tree = 330, edge.width=1.5, edge.color="black")
      tiplabels(phylo_tree_labels, cex = 1.5, frame = "none", col = by.population.labels, pch = 20)
      legend(x = 'right', inset=c(-0.35,0),
             legend = row.names(vp.BASIC.RESULTS)[1:vp.POP.NUM], bty = "n",
             col = unique(by.population.labels), cex=0.8, pch=c(16), xpd = NA)
      title(main = "Populations of Haplotypes", xpd = TRUE, outer = FALSE)
      title(main = vp.GENE.NAME, xpd = TRUE, line = 0.5, outer = FALSE)
    }
    vp.Phylogenetic.Tree()
    cat("\n")
    cat("This plot has been saved as a function, \"vp.Phylogenetic.Tree()\"\n")
    cat("You may need to use Zoom, or increase the width of the plot while\n")
    cat("exporting in order to see the legend clearly\n")
  }

  #HAPLOTYPE ACCUMULATION PLOT-------------------------------------------------------------------
  if (choice == 8) {
    cat("Haplotypes will be calculated using only the segregation sites where polymorphism\n")
    cat("is found in at least 'x' percent of the population at that site")
    threshold <- as.numeric(readline ("'x' <-  "))
    variant.codons.over.threshold <- vp.AA.VARIANT.TABLE[,!vp.AA.VARIANT.TABLE[1,]<=threshold/100]
    if (threshold == 0) variant.codons.over.threshold <- vp.AA.VARIANT.TABLE
    AAs.each.variable.codon.over.threshold <- vp.AAs.EACH.VARIABLE.CODON
    AAs.each.variable.codon.over.threshold <- AAs.each.variable.codon.over.threshold[ ,
      colnames(AAs.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    ref.each.variable.codon.over.threshold <- vp.REF.EACH.VARIABLE.CODON
    ref.each.variable.codon.over.threshold <- ref.each.variable.codon.over.threshold[ ,
      colnames(ref.each.variable.codon.over.threshold) %in% colnames(variant.codons.over.threshold)]
    AA.changes.per.hap <- rowSums(sweep(unique(AAs.each.variable.codon.over.threshold), 2,
                                        ref.each.variable.codon.over.threshold, FUN = `!=`,
                                        check.margin = FALSE) +0)
    min.hap.strings <- matrix(nrow = vp.SEQ.NUM)
    for (i in 1:vp.SEQ.NUM){
      min.hap.strings[i,1] <- paste(AAs.each.variable.codon.over.threshold[i,], collapse = "")}
    min.hap.count <- length(AA.changes.per.hap)
    min.hap.AAs <- unique(min.hap.strings)
    min.hap.table <- match(min.hap.strings, min.hap.AAs)
    accumulation.data <- matrix(nrow = length(min.hap.table),
      ncol = max(min.hap.table))
    rownames(accumulation.data) <- 1:length(min.hap.table)
    colnames(accumulation.data) <- 1:max(min.hap.table)
    for (i in 1:length(min.hap.table)){
     for (j in 1:max(min.hap.table)){
       accumulation.data[i,j] <- ifelse(colnames(accumulation.data)[j] == min.hap.table[i],1,0)
     }
    }
    sp_overall <- specaccum(accumulation.data, method = "rarefaction", permutations = 100)
    vp.Haplotype.Accumulation.Plot <<- function(){
      plot(sp_overall, col = "blue", lwd = 2, ci.type = "poly",
         ci.lty = 0, ci.col = "lightblue", xlab = paste(vp.GENE.NAME, "sequences"),
         ylab = "Haplotypes", main = paste(vp.GENE.NAME, "Haplotype Accumulation Plot"))
    }
    vp.Haplotype.Accumulation.Plot()
    cat("\n")
    cat("This plot has been saved as a function, \"vp.Haplotype.Accumulation.Plot()\"\n")

  }
  #POLYMORHPISM GRAPH (S)-----------------------------------------------------------------------------
  if (choice == 9) {
    window.size <- as.numeric(readline ("Sliding window size <-  "))
    step.size <- as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.S <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
     xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      sliding.window.S[1,i] <- seg.sites.num
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.S <- as.data.frame(t(sliding.window.S))


    vp.S.Graph <<- ggplot(data = sliding.window.S,
                          aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "Gold")+
      geom_area(fill = "Yellow")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump+50)*6/label.scaler),
                         labels = xlabels)+
      labs(subtitle = vp.GENE.NAME, y = "S", x = "Nucleotide", title = "Polymorphic Sites")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    cat("\n")
    cat("\n")
    cat("Saved as \"vp.S.Graph\"\n")
    return(suppressMessages(vp.S.Graph))
  }

  #TAJIMA'S D GRAPH---------------------------------------------------------------------------
  if (choice == 10) {
    window.size <- as.numeric(readline ("Sliding window size <-  "))
    step.size <- as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.TD <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
    xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      tpi <- sum(snps) /  (vp.SEQ.NUM) / vp.SEQ.NUM  * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      td <- ((tpi) -  (seg.sites.num / vp.A1)) /
        (sqrt(  (vp.E1*seg.sites.num)  + ( (vp.E2*seg.sites.num) * (seg.sites.num-1) ) ))
      sliding.window.TD[1,i] <- td + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.TD <- as.data.frame(t(sliding.window.TD))
    sliding.window.TD[is.na(sliding.window.TD)] <- 0

    td.sim.for.count <- matrix(ncol = 9)
    for (i in 1:72){
      if (vp.TD.SIM[i,1] <= vp.SEQ.NUM)
        if (vp.TD.SIM[(i+1),1] > vp.SEQ.NUM)
          td.sim.for.count <- vp.TD.SIM[i, ]
    }
    if (vp.SEQ.NUM >= 1000) td.sim.for.count <- vp.TD.SIM[73, ]

    vp.TD.Graph <<- suppressMessages(ggplot(data = sliding.window.TD,
                           aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "DarkRed")+
      geom_area(fill = "Red")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump+50)*6/label.scaler),
                         labels = xlabels)+
      scale_y_continuous(sec.axis = sec_axis(~.+0, breaks = (td.sim.for.count[2:9]),
                                             labels = c("0.1", "0.1",
                                                        "0.05", "0.05",
                                                        "0.01", "0.01",
                                                        "0.001", "0.001"),
                                             name = "p value for Tajima's D"))+
      geom_hline(yintercept = (td.sim.for.count[2:9]),
                 colour = "DarkRed", linetype = "dashed", size = 0.4)+
      labs(subtitle = vp.GENE.NAME, y = "Tajima's D", x = "Nucleotide", title = "Tajima's D")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
    cat("\n")
    cat("Saved as \"vp.TD.Graph\"\n")
    return(suppressMessages(vp.TD.Graph))
  }

  #NUCLEOTIDE DIVERSITY GRAPH-----------------------------------------------------------------------
  if (choice == 11) {
    window.size <- as.numeric(readline ("Sliding window size <-  "))
    step.size <- as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.pi <- matrix(ncol = number.of.windows, nrow = 1)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
    xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      sliding.window.pi[1,i] <- pi + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}

    sliding.window.pi <- as.data.frame(t(sliding.window.pi))
    sliding.window.pi[is.na(sliding.window.pi)] <- 0

    vp.pi.Graph <<- ggplot(data = sliding.window.pi,
                           aes(y = V1, x = 1:number.of.windows), na.rm=TRUE)+
      geom_line(colour = "DarkBlue")+
      geom_area(fill = "Blue")+
      theme_classic()+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump +50) *6/label.scaler),
                         labels = xlabels)+
      labs(subtitle = vp.GENE.NAME, y = "\u03a0", x = "Nucleotide", title = "Nucleotide Diversity")+
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    cat("\n")
    cat("Saved as \"vp.pi.Graph\"\n")
    return(suppressMessages(vp.pi.Graph))
  }

  #ALL SLIDING SCALES-------------------------------------------------------------------------------------
  if(choice == 12) {
    window.size <- as.numeric(readline ("Sliding window size <-  "))
    step.size <- as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.stats <- matrix(ncol = number.of.windows, nrow = 3)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler
    xlabels.values <- seq(from = xlabels.values[1], to = xlabels.values[6], by = xlabels.values[1])

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      tpi <- sum(snps) /  (vp.SEQ.NUM) / vp.SEQ.NUM  * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      td <- ((tpi) -  (seg.sites.num / vp.A1)) /
        (sqrt(  (vp.E1*seg.sites.num)  + ( (vp.E2*seg.sites.num) * (seg.sites.num-1) ) ))
      sliding.window.stats[1,i] <- pi + 0
      sliding.window.stats[2,i] <- td + 0
      sliding.window.stats[3,i] <- seg.sites.num
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}
    windowstart <- 1
    windowsites <- paste(windowstart, "-", (windowstart + window.size - 1), sep = "")
    for (i in 1:number.of.windows){
      windowsites <- c(windowsites, paste(windowstart, "-", (windowstart + window.size - 1), sep = ""))
      windowstart <- windowstart + step.size}
    windowsites <- windowsites[-1]
    colnames(sliding.window.stats) <- windowsites
    rownames(sliding.window.stats) <- c("Pi", "TD", "S")

    sliding.window.stats <- as.data.frame(t(sliding.window.stats))
    sliding.window.stats[is.na(sliding.window.stats)] <- 0

    Pi.scaler <- (max(sliding.window.stats$S)) / (max(sliding.window.stats$Pi))
    TD.scaler <- (max(sliding.window.stats$S)) / (max(sliding.window.stats$TD))


    sliding.window.stats$Pi <- sliding.window.stats$Pi*Pi.scaler
    sliding.window.stats$TD <- sliding.window.stats$TD*TD.scaler


    stacked.sliding.window.stats <- matrix(nrow = 3*number.of.windows, ncol = 3)
    counter <- 0
    for (i in 1:3){
      stacked.sliding.window.stats[(1+counter):(number.of.windows+counter),3] <-
        1:number.of.windows
      stacked.sliding.window.stats[(1+counter):(number.of.windows+counter),2] <-
        as.numeric(sliding.window.stats[ ,i])
      stacked.sliding.window.stats[(1+counter):(number.of.windows+counter),1] <-
        colnames(sliding.window.stats)[i]
      counter <- counter + number.of.windows
    }

    stacked.sliding.window.stats <- as.data.frame(stacked.sliding.window.stats)
    colnames(stacked.sliding.window.stats) <- c("variable", "value", "window")

    stacked.sliding.window.stats$value <-
      as.numeric(levels(stacked.sliding.window.stats$value))[stacked.sliding.window.stats$value]
    stacked.sliding.window.stats$window <-
      as.numeric(levels(stacked.sliding.window.stats$window))[stacked.sliding.window.stats$window]

    td.sim.for.count <- matrix(ncol = 9)
    for (i in 1:72){
      if (vp.TD.SIM[i,1] <= vp.SEQ.NUM)
        if (vp.TD.SIM[(i+1),1] > vp.SEQ.NUM)
          td.sim.for.count <- vp.TD.SIM[i, ]
    }
    if (vp.SEQ.NUM >= 1000) td.sim.for.count <- vp.TD.SIM[73, ]

    vp.Combined.Pop.Gen.Stats.Graph <<- ggplot(stacked.sliding.window.stats,
                                               aes(x = stacked.sliding.window.stats$window,
                                                   y = stacked.sliding.window.stats$value,
                                                   col = variable,
                                                   fill = variable), na.rm=TRUE)+
      geom_line(size = 0.7)+
      geom_area(alpha = 0.6, position = 'identity')+
      theme_classic()+
      labs(y = "Scaled Values",
           x = "Nucleotide", title = paste("Statistics Across", vp.GENE.NAME, "Scaled To TD"))+
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks = xlabels.values,
                         limits = c(0, (label.jump+50)*6/label.scaler),
                         labels = xlabels)+
      scale_y_continuous(breaks = c(min(stacked.sliding.window.stats$value),
                                    -2*TD.scaler, 0, 2*TD.scaler,
                                    max(stacked.sliding.window.stats$value)),
                         labels = c(paste(round(min(stacked.sliding.window.stats$value)/TD.scaler, 3)),
                                    "-2","0", "2",
                                    paste(round(max(stacked.sliding.window.stats$value)/TD.scaler, 3))),
                         sec.axis = sec_axis(~.+0, breaks = (td.sim.for.count[2:9])*TD.scaler,
                                             labels = c("0.1", "0.1",
                                                        "0.05", "0.05",
                                                        "0.01", "0.01",
                                                        "0.001", "0.001"),
                                             name = "p value for Tajima's D"))+
      geom_hline(yintercept = 0, colour = "Black", size = 0.7)+
      geom_hline(yintercept = (td.sim.for.count[2:9])*TD.scaler,
                 colour = "DarkRed", linetype = "dashed", size = 0.4)+
      scale_fill_manual(values = c("Blue", "Yellow", "Red"), name = "Statistic")+
      scale_colour_manual(values = c("DarkBlue", "Gold", "DarkRed"), name = "Statistic")
    cat("\n")
    cat("\n")
    cat("Confidence limits of Tajima's D determined by comparison to original \n")
    cat("simulation of beta distribution by F. Tajima (1989)\n")
    cat("\n")
    cat("Saved as \"vp.Combined.Pop.Gen.Stats.Graph\"\n")
    return(suppressMessages(vp.Combined.Pop.Gen.Stats.Graph))
  }

  #SLIDING SCALE DATA TABLE---------------------------------------------------------------------
  if(choice == 13) {
    window.size <- as.numeric(readline ("Sliding window size <-  "))
    step.size <- as.numeric(readline ("Step size <-  "))
    number.of.windows <- (ceiling((vp.SEQ.LENGTH / step.size) - (window.size / step.size)))
    sliding.window.stats <- matrix(ncol = number.of.windows, nrow = 3)
    windowstart <- 1
    windowend <- window.size

    label.jump <- ceiling((vp.SEQ.LENGTH / 100) / 7) * 100
    label.scaler <- vp.SEQ.LENGTH / number.of.windows
    xlabels.values <- matrix(ncol = 6)
    for (i in 1:6){
      xlabels.values[1,i] <- label.jump + (i-1)*label.jump
    }
    xlabels <- paste(xlabels.values)
    xlabels.values <- xlabels.values / label.scaler

    for (i in 1:number.of.windows){
      snps <- vp.GLOBAL.SNPS[windowstart:windowend]
      seg.sites <- matrix(0)
      for (j in 1:window.size){
        if (snps[j] != 0) {seg.sites <- cbind(seg.sites, as.matrix(paste(j)))}}
      seg.sites.num <- length(seg.sites) - 1
      pi.per.nuc <- ((snps /  (vp.SEQ.NUM)^2) / window.size) * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      tpi <- sum(snps) /  (vp.SEQ.NUM) / vp.SEQ.NUM  * (vp.SEQ.NUM/(vp.SEQ.NUM-1))
      pi <- sum(pi.per.nuc)
      td <- ((tpi) -  (seg.sites.num / vp.A1)) /
        (sqrt(  (vp.E1*seg.sites.num)  + ( (vp.E2*seg.sites.num) * (seg.sites.num-1) ) ))
      sliding.window.stats[1,i] <- seg.sites.num
      sliding.window.stats[2,i] <- pi + 0
      sliding.window.stats[3,i] <- td + 0
      windowstart <- windowstart + step.size
      windowend <- windowend + step.size
      cat("Calculating window number", i, "of", number.of.windows, "\r")}
    windowstart <- 1
    windowsites <- paste(windowstart, "-", (windowstart + window.size - 1), sep = "")
    for (i in 1:number.of.windows){
      windowsites <- c(windowsites, paste(windowstart, "-", (windowstart + window.size - 1), sep = ""))
      windowstart <- windowstart + step.size}
    windowsites <- windowsites[-1]
    colnames(sliding.window.stats) <- windowsites
    rownames(sliding.window.stats) <- c("S", "Pi", "TD")

    sliding.window.stats <- as.data.frame(t(sliding.window.stats))
    sliding.window.stats[is.na(sliding.window.stats)] <- 0

    td.sim.for.count <- matrix(ncol = 9)
    for (i in 1:72){
      if (vp.TD.SIM[i,1] <= vp.SEQ.NUM)
        if (vp.TD.SIM[(i+1),1] > vp.SEQ.NUM)
          td.sim.for.count <- vp.TD.SIM[i, ]
    }
    if (vp.SEQ.NUM >= 1000) td.sim.for.count <- vp.TD.SIM[73, ]

    TD_p_values <- matrix(nrow = number.of.windows)
    for (i in 1:number.of.windows){
      if (sliding.window.stats[i,3] < td.sim.for.count[2])
        TD_p_values[i,1] <- 0.1
      if (sliding.window.stats[i,3] > td.sim.for.count[3])
        TD_p_values[i,1] <- 0.1
      if (sliding.window.stats[i,3] < td.sim.for.count[4])
        TD_p_values[i,1] <- 0.05
      if (sliding.window.stats[i,3] > td.sim.for.count[5])
        TD_p_values[i,1] <- 0.05
      if (sliding.window.stats[i,3] < td.sim.for.count[6])
        TD_p_values[i,1] <- 0.01
      if (sliding.window.stats[i,3] > td.sim.for.count[7])
        TD_p_values[i,1] <- 0.01
      if (sliding.window.stats[i,3] < td.sim.for.count[8])
        TD_p_values[i,1] <- 0.001
      if (sliding.window.stats[i,3] > td.sim.for.count[9])
        TD_p_values[i,1] <- 0.001
    }
    vp.Sliding.Window.Stats <<- cbind(sliding.window.stats, TD_p_values)
    colnames(vp.Sliding.Window.Stats)[4] <- ("TD p value")
    vp.Sliding.Window.Stats[is.na(vp.Sliding.Window.Stats)] <- ">0.1"

  View(vp.Sliding.Window.Stats)
  cat("\n")
  cat("\n")
  cat("Confidence limits of Tajima's D determined by comparison to original \n")
  cat("simulation of beta distribution by F. Tajima (1989)\n")
  cat("\n")
  cat("Saved as \"vp.Sliding.Window.Stats\", use write.csv() to save to excel\n")
  cat("e.g. write.csv(vp.Sliding.Window.Stats, file = \"my.sliding.window.stats.in.excel\")")
  }
}
