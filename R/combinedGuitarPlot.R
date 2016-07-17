combinedGuitarPlot <- function(ct, comLength = c(0.136,0.459, 0.405)){ 
  # stack
  adjust=1
  
  # plot the figure
  # extract information of mRNA and lncRNA
  ct$weight <- ct$count # as numeric
  ct1 <- ct[ct$category=="mRNA",] # mRNA
  ct2 <- ct[ct$category=="lncRNA",] # lncRNA
  
  # disable notes
  pos=Feature=weight=NULL
  
  # remove DNA
  id1 <- which(match(ct1$comp,c("Front","Back")) >0 )
  ct1 <- ct1[-id1,]
  id2 <- which(match(ct2$comp,c("Front","Back")) >0 )
  ct2 <- ct2[-id2,]  
  
  # normalize feature
  featureSet <- as.character(unique(ct$Feature))
  for (i in 1:length(featureSet)) {
    id <- (ct1$Feature==featureSet[i])
    ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
    
    id <- (ct2$Feature==featureSet[i])
    ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
  }
  
  p2 <- 
    ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) + 
    ggtitle("Distribution on lncRNA")  +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
    xlab("") + 
    ylab("Frequency") +
    geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
    annotate("text", x = 0.5, y = -0.2, label = "lncRNA")+
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    theme(legend.position="bottom")
  
	weight <- comLength/sum(comLength)
      names(weight) <- c("5'UTR","CDS","3'UTR")

  # density
      cds_id <- which(ct1$comp=="CDS")
      utr3_id <- which(ct1$comp=="UTR3")
      utr5_id <- which(ct1$comp=="UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]
      
      # re-normalization
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      
      # stratch
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
      ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]
      
      p1 <- 
        ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        xlab("") + 
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = x[1]/2, y = -0.2, label = "5'UTR") +
        annotate("text", x = x[1] + weight[2]/2, y = -0.2, label = "CDS") +
        annotate("text", x = x[2] + weight[3]/2, y = -0.2, label = "3'UTR") + 
        theme(legend.position="bottom") +
        geom_vline(xintercept= x[1:2], linetype="dotted") +
        annotate("rect", xmin = 0, xmax = x[1], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")
  
  .multiplot(p1, p2, cols=2)
  
}