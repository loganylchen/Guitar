\name{combinedGuitarPlot}
\alias{combinedGuitarPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
combinedGuitarPlot
}
\description{
combine multiple GuitarPlots, especially when they are from different species and cannot directly handled with GuitarPlot function together.
}
\usage{
combinedGuitarPlot(ct)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{ct}{
Count of overlapping features, which is can be the output from GuitarPlot function. See example. 
}
}
\details{
combine multiple GuitarPlots, especially when they are from different species.
}
\value{
A figure showing the transcriptomic distribution of the genomic features will be generated. Post-editing with Adobe Illustrator or other graphic software is recommended.
}

\examples{

# read genomic features
narrowPeak <- system.file(
  "extdata", "m6A_hg19_1000peaks_macs2.narrowPeak", 
  package="Guitar")

# bam imported as GAlignments
m6A_Bcell <- narrowPeaktoGRanges(narrowPeak) 

# generate a list of genomic features
m6A_Bcell_1 <- m6A_Bcell[1:300]
m6A_Bcell_2 <- m6A_Bcell[301:600]
m6A_Bcell_3 <- m6A_Bcell[601:900]
m6A_Bcell_4 <- m6A_Bcell[201:900]

# Make Guitar coordiantes
txdb_file <- system.file("extdata", "hg19_toy.sqlite", 
                         package="Guitar")
txdb <- loadDb(txdb_file)
gc_txdb <- makeGuitarCoordsFromTxDb(txdb,noBins =10)

# Guitar plot 1
feature1_hg19 <- list(m6A_Bcell_1, m6A_Bcell_2) 
names(feature1_hg19) <- c("m6A_1","m6A_2")
ct1 <- GuitarPlot(feature1_hg19, returnCount = TRUE,
           GuitarCoordsFromTxDb = gc_txdb)

# Guitar plot 2
feature2_hg19 <- list(m6A_Bcell_3, m6A_Bcell_4) 
names(feature2_hg19) <- c("m6A_3","m6A_4")
ct2 <- GuitarPlot(feature2_hg19, returnCount = TRUE,
           GuitarCoordsFromTxDb = gc_txdb)

# combine two Gutiar Plot
ct <- rbind(ct1,ct2)
combinedGuitarPlot(ct)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
