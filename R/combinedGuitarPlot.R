combinedGuitarPlot <- function(ct){ 
  # stack
  
  
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
  
  pos_adjust <- match(ct1$comp,c("UTR5","CDS","UTR3"))-1
  ct1$pos <- ct1$pos + pos_adjust
  p1 <- 
    ggplot(ct1, aes(x=pos, group=Feature,colour=factor(Feature), weight=3*weight))  +
    ggtitle("Distribution on mRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
    xlab("") + 
    ylab("Frequency") +
    geom_density(adjust=adjust,aes(fill=factor(Feature)),alpha=0.2) +
    annotate("text", x = 0.5, y = -0.2, label = "5'UTR") +
    annotate("text", x = 1.5, y = -0.2, label = "CDS") +
    annotate("text", x = 2.5, y = -0.2, label = "3'UTR") + 
    geom_vline(xintercept=1:2, linetype="dotted") + 
    theme(legend.position="bottom") +
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    annotate("rect", xmin = 2, xmax = 3, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    annotate("rect", xmin = 1, xmax = 2, ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")
  
  
  .multiplot(p1, p2, cols=2)
  
}