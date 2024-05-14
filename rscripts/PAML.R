# Set workinng directory
setwd("/Users/anjaligupta/Library/CloudStorage/OneDrive-UniversityofKansas/PhD/SR_drive/PAML")

# Import data
library(readxl)
Book2 <- read_excel("Book2.xlsx")

# Filter NAs
Book2 <- subset(Book2, is.na(omega_dNdS)==FALSE)
Book2$omega <- ifelse(Book2$omega_dNdS > 1, ">1",
                       ifelse(Book2$omega_dNdS < 1, "<1", 
                              ifelse(Book2$omega_dNdS == 1, "1", "NA")))


require(GenomicFeatures)
txdb <- makeTxDbFromGFF("Daff_ST_ChrX.gtf", format="gtf")
genes <- as.data.frame(genes(txdb))

exons <- as.data.frame(exonsBy(txdb), by = c("gene"))
transcripts <- as.data.frame(transcriptsBy(txdb, by = "gene"))
gene <- as.data.frame(genes, row.names = NULL)

transcripts$Gene_ID <- transcripts$tx_name

Book <- merge(Book2, transcripts, by="Gene_ID", all.x=TRUE)

# Visualization
library(ggplot2)

ggplot(Book, aes(x=as.numeric(start),
                 y=as.numeric(log10(omega_dNdS)),
                 color=as.factor(omega))) +
  geom_point(alpha=0.5) + 
  labs(x="Position on X chromosome",
       y="log10 Omega (dN/dS)",
       color="omega (dN/dS)") +
  theme_bw() +
  theme(legend.position = "bottom",
        title = element_text(size = 14, face = "bold"),
        text = element_text(size = 14))

ggplot(Book, aes(x=as.numeric(start),
                  y=as.numeric(log10(omega_dNdS)),
                  color=as.factor(omega),
                 label = group_name)) +
  geom_point(alpha=0.5, position = position_jitterdodge()) +
  geom_text(size = 2, vjust = -0.5, position = position_jitterdodge()) + 
  labs(x="Position on X",
       y="Omega (dN/dS)",
       color="omega (dN/dS)") +
  theme_bw() +
  theme(legend.position = "bottom")

