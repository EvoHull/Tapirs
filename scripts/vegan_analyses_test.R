#Set working directory to where the r script is
setwd("~/github/Tapirs/scripts")

#Read the converted tsv files
kraken2_table <- t(read.table(file = '~/Downloads/table.from_biom.tsv', 
                           sep = '\t', header = TRUE, row.names = 1))

blast_table <- t(read.table(file = '~/Downloads/tableblast.from_biom.tsv', 
                            sep = '\t', header = TRUE, row.names = 1))
#Now lets test this in vegan

#install and load packages (perhaps look for a way to check for package and install if missing)
install.packages("vegan")
install.packages("ape")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggfortify")
library("vegan")
library("ape")
library("dplyr")
library("ggplot2")
library("ggfortify")
#This would be an example of how to check for packages and install if missing
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

#testing that vegan works on this

#beta diversity
kraken_beta_diversity <- vegdist(kraken2_table, method="bray")
blast_beta_diversity <- vegdist(blast_table, method="bray")

#pca
kraken_pca <- prcomp(kraken_beta_diversity)
blast_pca <- prcomp(blast_beta_diversity)

#How do i combine kraken and blast output

autoplot(kraken_pca, blast_pca)
#Appears to NOT PLOT both - doesn't work, gives error with color. so I think I have to add the method to the sample name and then merge the tables

row.names(kraken2_table) <- paste0(row.names(kraken2_table), "_kraken")
row.names(blast_table) <- paste0(row.names(blast_table), "_blast")
#adds the method of the result to the sample names (to differentiate between blast and kraken outputs)

#in order to then bind the 2 tables the rbind.fill command from the pylr package to "fill" the empty columns - the columns they don't have in common
install.packages("plyr")
library("plyr")
all_methods_table <- rbind.fill(data.frame(kraken2_table), data.frame(blast_table))

#Repeating the previous analyses with this new table
all_methods_betadiv <- vegdist(all_methods_table, method="bray", na.rm = TRUE)
all_methods_pca <- prcomp(all_methods_betadiv)
plot(all_methods_pca$x[,1], all_methods_pca$x[,2], pch=20, main="PCA plot of Kraken2 and BLAST results", 
     col=rep(c("red","blue"), each=nrow(kraken2_table)), ylab="PCA2", xlab="PCA1")
legend("bottomright", legend=c("Kraken2", "BLAST"), title = "Analysis method", fill = c("red","blue"), cex=0.8)

##Saving the PCA plot
# Step 1: Call the pdf command to start the plot
png(file = "~/github/Tapirs/results/pca_plot_test.png",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches

# Step 2: Create the plot with R code
plot(all_methods_pca$x[,1], all_methods_pca$x[,2], pch=20, main="PCA plot of Kraken2 and BLAST results", 
     col=rep(c("red","blue"), each=nrow(kraken2_table)), ylab="PCA2", xlab="PCA1")
legend("bottomright", legend=c("Kraken2", "BLAST"), title = "Analysis method", fill = c("red","blue"), cex=0.8)

# Step 3: Run dev.off() to create the file!
dev.off()


