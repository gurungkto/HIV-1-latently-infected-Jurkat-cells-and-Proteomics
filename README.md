# HIV-1 latently infected Jurkat cells and Proteomics
Approximately 1 gram cell pellet per cell line was processed for HLA-peptide extraction. Peptides were extracted from the HLA molecules, fractionated by RP-HPLC, and ran on nano-LC MS/MS system (TTOF 5600). Raw data from the mass spectrometer was processed within PEAKS software at 1% FDR to generate peptide hits based on peptide spectrum matches.
A total of **595** peptides were identified in **uninfected Jurkat E6-1** cell line of which **272** were unique peptides.

```{r, echo=TRUE, warning=FALSE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
data.e61 <- read.csv("protein-peptides.csv", header = T)
table.e61 <- table(data.e61$Length)
table.e61 <- data.frame(table.e61)
library(ggplot2)
ggplot(table.e61, aes(table.e61$Var1, table.e61$Freq, group = 1)) + geom_point() + 
  geom_line(color = "red") + 
  labs(x = "Peptide length" , y = "Frequency") +
  theme(axis.line = element_line(size = 1),plot.title = element_text(hjust = 0.5)) +
  ggtitle("Frequency of #mer Jurkat E6-1 MHC I peptides")
```

A total of **636** peptides were identified in HIV-1 infected **Jurkat DF731** cell line of which **269** were unique peptides.

```{r, echo=TRUE, warning=FALSE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF731 W632 1g")
data.731 <- read.csv("protein-peptides.csv", header = T)
table.731 <- table(data.731$Length)
table.731 <- data.frame(table.731)
library(ggplot2)
ggplot(table.731, aes(table.731$Var1, table.731$Freq, group = 1)) + geom_point() + 
  geom_line(color = "red") + 
  labs(x = "Peptide length" , y = "Frequency") +
  theme(axis.line = element_line(size = 1),plot.title = element_text(hjust = 0.5)) +
  ggtitle("Frequency of #mer Jurkat DF-731 MHC I peptides")
```

A total of **427** peptides were identified in **HIV-1 infected Jurkat DFRC3** cell line of which **170** were unique peptides.

```{r, echo=TRUE, warning=FALSE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF-RC3 W632 1g")
data.RC3 <- read.csv("protein-peptides.csv", header = T)
table.RC3 <- table(data.RC3$Length)
table.RC3 <- data.frame(table.RC3)
library(ggplot2)
ggplot(table.RC3, aes(table.RC3$Var1, table.RC3$Freq, group = 1)) + geom_point() + 
  geom_line(color = "red") + 
  labs(x = "Peptide length" , y = "Frequency") +
  theme(axis.line = element_line(size = 1),plot.title = element_text(hjust = 0.5)) +
  ggtitle("Frequency of #mer Jurkat DF-RC3 MHC I peptides")
```


### Comparision between **uninfected Jurkat E6-1** and **HIV-1 infected Jukat DF731** Cell lines

```{r, echo=TRUE, warning=FALSE,  message=FALSE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
data.e61peptide <- read.csv("peptide.csv", header = T)

setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF731 W632 1g")
data.731peptide <- read.csv("peptide.csv", header = T)

matching.peptides <- data.731peptide$Peptide[data.731peptide$Peptide %in% data.e61peptide$Peptide]
length(matching.peptides)
```

Between uninfected Jurkat E6-1 and infected Jurkat DF731, **155** were overlapping peptides.

```{r, echo=TRUE, message=FALSE, warning=FALSE}
library(VennDiagram)
draw.pairwise.venn(272, 269, 155, category = c("Jurkat E6-1", "HIV-1 infected Jurkat DF731"),
                   fill = c("green", "red"), alpha = rep(0.5,2), cat.pos = c(0,0), 
                   cat.dist = rep(0.025,2))
```

That means 117 peptides were *unique* to Jurkat E6-1 and 114 peptides were *unique* to HIV-1 infected Jurkat DF731.
Performing database match using Anthony Purcell's Latent HIV-1 database (upregulated genes/proteins) http://hivlatency.erc.monash.edu/searchDBwithExpression?SearchDataset=Mo-HIVvsMockinLatency&SearchEffect=up&SubmitExpression=Submit, 45 peptides in Jurkat E6-1 matched the database whereas 57 peptides in Jurkat DF731 matched the database. Further analysis showed only 2 overlapping genes/proteins between the Jurkat E6-1 and Jurkat DF731 which means 43 peptides were *unique* to Jurkat E6-1 and 55 peptides were *unique* to HIV-1 infected Jurkat DF731.

The reason why the 2 overlapping genes/proteins were observed between the Jurkat E6-1 and Jurkat DF731 was because the two genes/proteins (TNKS2 and TRIP12) were associated with different peptide sequences within each cell line. In other words, different sequences were sampled from the same proteins between the two cell lines.

```{r, echo=TRUE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
data.vlookupUniqueE61 <- read.csv("vlookupUniqueE61vsDF731.csv", header = T)
data.vlookupUniqueE61 <- data.frame(lapply(data.vlookupUniqueE61, as.character), 
                                    stringsAsFactors = F)
dt <- data.frame(matrix(unlist(strsplit(data.vlookupUniqueE61$Proteins, split = " "))))
colnames(dt) <- "strings"
dt <- as.data.frame(dt[(grep("^GN=", dt$strings)),])
colnames(dt) <- "Genes"
dt <- data.frame(substr(dt$Genes, 4, length(dt$Genes)))
colnames(dt) <- "Genes"

setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF731 W632 1g")
data.vlookupUniqueDF731 <- read.csv("vlookupUniqueDF731vsE61.csv", header = T)
data.vlookupUniqueDF731 <- data.frame(lapply(data.vlookupUniqueDF731, as.character), 
                                      stringsAsFactors = F)
dt.731 <- data.frame(matrix(unlist(strsplit(data.vlookupUniqueDF731$PROTEINS, split = " "))))
colnames(dt.731) <- "strings"
dt.731 <- as.data.frame(dt.731[(grep("^GN=", dt.731$strings)),])
colnames(dt.731) <- "Genes"
dt.731 <- data.frame(substr(dt.731$Genes, 4, length(dt.731$Genes)))
colnames(dt.731) <- "Genes"

# Read hiv-1 latency associated genes
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan")
latency <- read.csv("HIV-1 latency upregulated peptides.csv",  header = T)

# compare to E61 genes
matching.E61.genes <- dt$Genes[dt$Genes %in% latency$Gene.Name]
length(matching.E61.genes) # Number of genes that matched the HIV-1 latency database in Jurkat E6-1 
matching.E61.genes <- data.frame(matching.E61.genes)
colnames(matching.E61.genes) <- "Genes"

# compare to DF731 genes
matching.DF731.genes <- dt.731$Genes[dt.731$Genes %in% latency$Gene.Name]
length(matching.DF731.genes)# Number of genes that matched the HIV-1 latency database in 
                            # HIV-1 infected Jurkat DF731
matching.DF731.genes <- data.frame(matching.DF731.genes)
colnames(matching.DF731.genes) <- "Genes"

# next question.... are those that matched in each group the same peptides...what are unique to df731?

matching.genes.in.both.samples <- matching.E61.genes$Genes[matching.E61.genes$Genes %in% matching.DF731.genes$Genes]
length(matching.genes.in.both.samples) # only 2 HIV-1 associated genes were overlapping 
                                       # between E6-1 and DF731
matching.genes.in.both.samples <- data.frame(matching.genes.in.both.samples)
colnames(matching.genes.in.both.samples) <- "Genes"

unique.matching.E61.genes <- matching.E61.genes[!(matching.E61.genes$Genes %in% matching.genes.in.both.samples$Genes),]
length(unique.matching.E61.genes) # 43 genes in E61 uniquely match the HIV latency database
unique.matching.E61.genes <- data.frame(unique.matching.E61.genes)
colnames(unique.matching.E61.genes) <- "Genes"

unique.matching.DF731.genes <- matching.DF731.genes[!(matching.DF731.genes$Genes %in% matching.genes.in.both.samples$Genes),]
length(unique.matching.DF731.genes) # 55 genes in DF731 sample uniquely match the HIV latency database
unique.matching.DF731.genes <- data.frame(unique.matching.DF731.genes)
colnames(unique.matching.DF731.genes) <- "Genes"

```

The latent HIV-1 associated host peptides that were *unique* to each group are presented in tables below for each cell line.

### Latent HIV-1 associated host peptides in Jurkat E6-1 as compared to Jurkat DF731

```{r, echo=TRUE}
library(DT)
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
unique.E61.peptides <- read.csv("Unique-E61-peptides-vs-DF731.csv", header = T)
unique.matching.E61.proteins <- unique.E61.peptides[which(unique.E61.peptides$Genes %in% unique.matching.E61.genes$Genes),]
datatable(unique.matching.E61.proteins[ , c(1,2,3,4,5,6,7,10,15,16)], filter = "top")

```


### Latent HIV-1 associated host peptides in HIV-1 infected Jurkat DF731 as compared to Jurkat E6-1

```{r, echo=TRUE}
library(DT)
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF731 W632 1g")
unique.DF731.peptides <- read.csv("Unique-DF731-peptides-vs-E61.csv", header = T)
unique.matching.DF731.proteins <- unique.DF731.peptides[which(unique.DF731.peptides$Genes %in% unique.matching.DF731.genes$Genes),]
datatable(unique.matching.DF731.proteins[ , c(1,2,3,4,5,6,7,10,15,16)], filter = "top")
```

## Similarly, comparing the  **uninfected Jurkat E6-1** and **HIV-1 infected Jukat DFRC3** Cell lines

```{r, echo=TRUE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
data.e61peptide <- read.csv("peptide.csv", header = T)

setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF-RC3 W632 1g")
data.RC3peptide <- read.csv("peptide.csv",  header = T)

matching.peptides <- data.RC3peptide$Peptide[data.RC3peptide$Peptide %in% data.e61peptide$Peptide]
```

Between uninfected Jurkat E6-1 and infected Jurkat DFRC3, **101** were overlapping peptides.

```{r, echo=TRUE, message=FALSE, warning=FALSE}
library(VennDiagram)
draw.pairwise.venn(272, 170, 101, category = c("Jurkat E6-1", "HIV-1 infected Jurkat DFRC3"),
                   fill = c("green", "red"), alpha = rep(0.5,2), cat.pos = c(0,0), 
                   cat.dist = rep(0.035,2))
```

That means 171 peptides were *unique* to Jurkat E6-1 and 69 peptides were *unique* to HIV-1 infected Jurkat DFRC3.
Performing database match using Anthony Purcell's Latent HIV-1 database (upregulated genes/proteins) http://hivlatency.erc.monash.edu/searchDBwithExpression?SearchDataset=Mo-HIVvsMockinLatency&SearchEffect=up&SubmitExpression=Submit, 66 peptides in Jurkat E6-1 matched the database whereas 28 peptides in Jurkat DF731 matched the database. Further analysis showed only 1 overlapping genes/proteins between the Jurkat E6-1 and Jurkat DF731 which means 65 peptides were *unique* to Jurkat E6-1 and 27 peptides were *unique* to HIV-1 infected Jurkat DFRC3.

Again, the reason why the 1 overlapping genes/protein was observed between the Jurkat E6-1 and Jurkat DFRC3 was because the gene/protein (TNKS2) was associated with different peptide sequences within each cell line. In other words, different sequences were sampled from the same protein between the two cell lines.

```{r, echo=TRUE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
data.vlookupUniqueE61 <- read.csv("vlookupUniqueE61vsRC3.csv", header = T)
data.vlookupUniqueE61 <- data.frame(lapply(data.vlookupUniqueE61, as.character), stringsAsFactors = F)

dt <- data.frame(matrix(unlist(strsplit(data.vlookupUniqueE61$Proteins, split = " "))))
colnames(dt) <- "strings"
dt <- as.data.frame(dt[(grep("^GN=", dt$strings)),])
colnames(dt) <- "Genes"
dt <- data.frame(substr(dt$Genes, 4, length(dt$Genes)))
colnames(dt) <- "Genes"

### 
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF-RC3 W632 1g")
data.vlookupUniqueRC3 <- read.csv("vlookupUniqueRC3vsE61.csv", header = T)
data.vlookupUniqueRC3 <- data.frame(lapply(data.vlookupUniqueRC3, as.character), stringsAsFactors = F)

dt.RC3 <- data.frame(matrix(unlist(strsplit(data.vlookupUniqueRC3$Proteins, split = " "))))
colnames(dt.RC3) <- "strings"
dt.RC3 <- as.data.frame(dt.RC3[(grep("^GN=", dt.RC3$strings)),])
colnames(dt.RC3) <- "Genes"

dt.RC3 <- data.frame(substr(dt.RC3$Genes, 4, length(dt.RC3$Genes)))
colnames(dt.RC3) <- "Genes"

# Read hiv-1 latency associated genes
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan")
latency <- read.csv("HIV-1 latency upregulated peptides.csv",  header = T)

# compare to E61 genes
matching.E61.genes <- dt$Genes[dt$Genes %in% latency$Gene.Name]
length(matching.E61.genes) # 66 peptides from E6-1 matched the HIV-1 latency database
matching.E61.genes <- data.frame(matching.E61.genes)
colnames(matching.E61.genes) <- "Genes"

# compare to DFRC3 genes
matching.DFRC3.genes <- dt.RC3$Genes[dt.RC3$Genes %in% latency$Gene.Name]
length(matching.DFRC3.genes) # 28 peptides from RC3 matched the HIV-1 latency database
matching.DFRC3.genes <- data.frame(matching.DFRC3.genes)
colnames(matching.DFRC3.genes) <- "Genes"

# next question.... are those that matched in each group the same peptides...what are unique to dfRC3?

matching.genes.in.both.samples <- matching.E61.genes$Genes[matching.E61.genes$Genes %in% matching.DFRC3.genes$Genes]
length(matching.genes.in.both.samples) # Only 1 HIV-1 latency associated gene/protein overlapped                                              # between E61 and RC3 cell lines
matching.genes.in.both.samples <- data.frame(matching.genes.in.both.samples)
colnames(matching.genes.in.both.samples) <- "Genes"

# only 1 gene (TNKS2) is common to both uninfected and infected jurkats that matched 
# the HIV latency associated genes

unique.matching.E61.genes <- matching.E61.genes[!(matching.E61.genes$Genes %in% matching.genes.in.both.samples$Genes),]
length(unique.matching.E61.genes) # 65
unique.matching.E61.genes <- data.frame(unique.matching.E61.genes)
colnames(unique.matching.E61.genes) <- "Genes"

unique.matching.DFRC3.genes <- matching.DFRC3.genes[!(matching.DFRC3.genes$Genes %in% matching.genes.in.both.samples$Genes),]
length(unique.matching.DFRC3.genes) #27
unique.matching.DFRC3.genes <- data.frame(unique.matching.DFRC3.genes)
colnames(unique.matching.DFRC3.genes) <- "Genes"
```

The latent HIV-1 associated host peptides that were *unique* to each group are presented in tables below for each cell line.

### Latent HIV-1 associated host peptides in Jurkat E6-1 as compared to Jurkat DFRC3

```{r, echo=TRUE}
library(DT)
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT E61 W632 1g")
unique.E61.peptides <- read.csv("Unique-E61-peptides-vs-RC3.csv", header = T)
unique.matching.E61.proteins <- unique.E61.peptides[which(unique.E61.peptides$Genes %in% unique.matching.E61.genes$Genes),]
datatable(unique.matching.E61.proteins[ , c(1,2,3,4,5,6,7,10,15,16)], filter = "top")

```


### Latent HIV-1 associated host peptides in HIV-1 infected Jurkat DFRC3 as compared to Jurkat E6-1

```{r, echo=TRUE}
library(DT)
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF-RC3 W632 1g")
unique.DFRC3.peptides <- read.csv("Unique-DFRC3-peptides-vs-E61.csv", header = T)
unique.matching.DFRC3.proteins <- unique.DFRC3.peptides[which(unique.DFRC3.peptides$Genes %in% unique.matching.DFRC3.genes$Genes),]
datatable(unique.matching.DFRC3.proteins[ , c(1,2,3,4,5,6,7,10,15,16)], filter = "top")
```

## How many genes overlap between the Jurkat DF731 and Jurkat DFRC3?

20 genes (18 peptides) overlapped the HIV-1 infected Jurkat DF731 and Jurkat DFRC3

```{r, echo=TRUE}
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF731 W632 1g")
data.jrkt731 <- read.csv("Unique-HIV1-related-DF731-peptides.csv", header = T)
setwd("C:/Users/HGURUNG1/Desktop/HG/Projects/West and Belshan/JRKT DF-RC3 W632 1g")
data.jrktRC3 <- read.csv("Unique-HIV1-related-DFRC3-peptides.csv", header = T)

# How many genes overlap ?

overlap.genes <- data.jrkt731$Genes[data.jrkt731$Genes %in% data.jrktRC3$Genes]
length(unique(overlap.genes))

overlap.peptides <- data.jrkt731$Peptide[data.jrkt731$Peptide %in% data.jrktRC3$Peptide]
length(unique(overlap.peptides))

# identify the two peptides that do not overlap

# overlapping.peptides.from.genes <- data.jrkt731[data.jrkt731$Genes %in% overlap.genes, ]
# two.nonoverlapping.peptides <- overlapping.peptides.from.genes[!overlapping.peptides.from.genes$Peptide %in% overlap.peptides,]
# two.nonoverlapping.peptides$Peptide

# The two peptides were from KIF18B gene that were not present in RC3

```

#### The table below contains the final list of 18 unique peptides/genes/proteins that overlapped HIV-1 infected Jurkats DF731 and Jurkat DFRC3

```{r, echo=TRUE}
library(DT)
datatable(data.jrkt731[data.jrkt731$Peptide %in% overlap.peptides, c(1,2,3,4,5,6,7,10,15,16)], 
          filter = "top")
```
