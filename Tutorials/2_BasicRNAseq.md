MTSchmitz

Basic RNAseq Analysis
=====================

See the end of this page for vocab and other hints \#\# Loading a File

The first step will be to download the example expected counts file GhostOfMT/Data/GSE90053\_H1\_EC.txt. You can do this by downloading the whole repository from github or by downloading the individual files. When you download it, the file will probably at the path ~/Downloads/GSE90053\_H1\_EC.txt. You can move it wherever you want and then right click and "Get Info" to find the path to the file (without the file name itself tacked on the end, which you will need to add after you copy and paste the path from there).

We will start by installing the packages we need:

``` r
source("http://bioconductor.org/biocLite.R")
#install packages if not already installed
biocLite("EBSeq")
biocLite("rafalib")
biocLite("RColorBrewer")
biocLite("ggplot2")
```

Once we have those, we can load up the RNAseq file, which you downloaded from the Data folder in this Github Repo

read.csv2 is a versatile function that you can use to read data files that are delineated with a bunch of stuff. It's not magic though, and if your files are irregular or tampered with, it won't save you.

The parameters of this function are very important. header is a T/F parameter stating whether or not the first line of the file is the column names. row.names is a numeric parameter asking which column should be rownames (enter 0 if none, or you have duplicated gene names). The sep value is what the separator is for this file type. "" is code for tab. You can do " " for space, or "," for comma-separated variable (csv) files. RSEM outputs tab-separated variable files (.tab or .tsv) and this line should work on RSEM output files out of the box.

type ?read.csv2 into the prompt in RStudio and it will give you more info!

You should change the path to the file (pathToFile) within the quotes to the path on your computer.

``` r
#load the packages you'll need:
#for normalizing
library(EBSeq)
#for changing view layout
library(rafalib)
#for the color palette
library(RColorBrewer)
#for fancy plots
library(gplots)
library(ggplot2)
```

``` r
pathToFile <- "~/code/GhostOfMT/Data/GSE90053_H1_EC.txt"
#NOTE THE PARAMETERS
counts <- read.csv2(file=pathToFile,header = T,row.names = 1,sep = "\t")
#counts <- counts[,-ncol(counts)] #uncomment this line if the last column is gene name description text
#It removes the last column of the data from the data set
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
head(counts)
```

    ##        H1Rosa7_d0 H1Rosa7_d1 H1Rosa7_d2 H1Rosa7_d3 H1Rosa7_d4 H1Rosa7_d5
    ## A1BG        12.00       7.00      13.00       4.00       1.00          9
    ## A1CF         0.00       0.00       2.87       0.00       0.00          0
    ## A2LD1        0.00       1.00       6.00       9.00       7.00          5
    ## A2M          0.00       0.00       0.00       0.00       0.00          0
    ## A2ML1       20.08       7.64       8.40       4.04       1.14          0
    ## A4GALT      19.00      19.00      34.00      41.00      21.00         19
    ##        H1Rosa7_d6 H1Rosa7_d7 H1Rosa7_d8 H1Rosa7_d10 H1Rosa7_d12
    ## A1BG        14.00       6.00       5.00           1       11.00
    ## A1CF         0.00       1.05       0.00           0        7.51
    ## A2LD1       13.00       3.00       2.00          10        6.00
    ## A2M          5.00       0.00       0.00           0        0.00
    ## A2ML1        0.90       0.00       1.89           0        7.82
    ## A4GALT       6.93      20.00       4.00           0        0.00
    ##        H1Rosa7_d14 H1Rosa7_d16 H1Rosa7_d18 H1Rosa7_d20 H1Rosa7_d21
    ## A1BG         12.00        8.00          15       27.00          12
    ## A1CF          0.58        0.07           0        0.00           0
    ## A2LD1         0.00        2.00           3        6.00           0
    ## A2M           9.00       29.00          53       83.00         160
    ## A2ML1         0.00        1.02           5        1.02           0
    ## A4GALT        0.00        0.00           1        0.00           1
    ##        H1Rosa7_d22 H1Rosa7_d24 H1Rosa7_d26 H1Rosa7_d30 H1Rosa7_d32
    ## A1BG             0           3        5.00          10        4.00
    ## A1CF             0           0        1.01           1        0.00
    ## A2LD1            3           0        1.00           2        2.00
    ## A2M            118          31       19.00          16       26.00
    ## A2ML1            0           0        0.00           0        2.21
    ## A4GALT           3           0        0.00           0        0.00
    ##        H1Rosa7_d34 H1Rosa7_d36 H1Rosa7_d38 H1Rosa7_d40 H1Rosa7_d42
    ## A1BG          5.00           8           7          15        8.00
    ## A1CF          0.00           0           0           0        3.65
    ## A2LD1         2.00           3           0          11        0.00
    ## A2M          62.00          67          96          80       70.00
    ## A2ML1         7.85           0           0           0        5.38
    ## A4GALT        0.00           0           3           0        0.00

Most of data you will be dealing with will probably be RSEM outputs. These outputs are usually structured as a directory with a name like "Sub\_0146\_Developmental\_Clock\_Human\_hg19\_\_890bf46a388bf6b4". Inside will be a mess of files, the most likely ones you will be working with are genes.no\_mt.**ec**.tab for counts data (NOT NORMALIZED) or TPMs genes.no\_mt.**tpm**.renorm.tab.

**QUICK NOTE:** If you have opened your files in Excel, always check the genes "MARCH1" and "SEPT1" to make sure that these gene names haven't been transformed into dates. If they have and you share/publish this data, people will instantly know you are a noob and discount any analysis you do!

You can get the number of mapped reads in each sample by getting the sum of the count numbers for each column (Genes are rows, Samples are columns)

``` r
colSums(counts)
```

    ##  H1Rosa7_d0  H1Rosa7_d1  H1Rosa7_d2  H1Rosa7_d3  H1Rosa7_d4  H1Rosa7_d5 
    ##   2535843.2   2395744.7   3123291.4   2277212.7   2140561.8   1821945.8 
    ##  H1Rosa7_d6  H1Rosa7_d7  H1Rosa7_d8 H1Rosa7_d10 H1Rosa7_d12 H1Rosa7_d14 
    ##   2485047.4   3322224.2   2669096.6   1711863.8   1426164.3   1807225.6 
    ## H1Rosa7_d16 H1Rosa7_d18 H1Rosa7_d20 H1Rosa7_d21 H1Rosa7_d22 H1Rosa7_d24 
    ##   1923557.4   1928532.0   1708037.4   1644661.4   1316266.1   1451445.5 
    ## H1Rosa7_d26 H1Rosa7_d30 H1Rosa7_d32 H1Rosa7_d34 H1Rosa7_d36 H1Rosa7_d38 
    ##   1311997.9   1520466.2   1509128.7   1550805.8   1522965.3   1956161.9 
    ## H1Rosa7_d40 H1Rosa7_d42 
    ##   1529980.2    908853.4

If everything looks good, now we have to normalize the counts data, because the data is Expected Counts. TPMs are already normalized, so you should skip this step if you're working with TPMs. **YOU HAVE TO NORMALIZE COUNTS OR ALL LIFE IS MEANINGLESS** (sample comparisons will be confounded by read depth)

``` r
normCounts <- GetNormalizedMat(counts, MedianNorm(counts))
#show plots in 1 row 2 columns
mypar(1,2)
plot(counts["PAX6",],main = "PAX6 Unnormalized Reads")
plot(normCounts["PAX6",],main = "PAX6 Normalized Reads",ylab="Normalized counts")
```

![](2_BasicRNAseq_files/figure-markdown_github/Normalizing-1.png)

It also might be helpful to parse the timepoints from the sample names. (Quite helpful.) The split parameter is a short bit of text that divides the different parts of longer strings of characters (sample names). This is why smart, standard naming of samples in the database is paramount!

``` r
splitList <- strsplit(x = colnames(normCounts),split = "_d")
head(splitList)
```

    ## [[1]]
    ## [1] "H1Rosa7" "0"      
    ## 
    ## [[2]]
    ## [1] "H1Rosa7" "1"      
    ## 
    ## [[3]]
    ## [1] "H1Rosa7" "2"      
    ## 
    ## [[4]]
    ## [1] "H1Rosa7" "3"      
    ## 
    ## [[5]]
    ## [1] "H1Rosa7" "4"      
    ## 
    ## [[6]]
    ## [1] "H1Rosa7" "5"

This gives us something tricky to use. First we use sapply to change the list to a matrix, where the first row is the cell type and the second row is the time point:

``` r
splitMatrix <- sapply(splitList, "[", 1:max(sapply(splitList,length)) )
print(dim(splitMatrix))
```

    ## [1]  2 26

``` r
print(splitMatrix)
```

    ##      [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]     
    ## [1,] "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7"
    ## [2,] "0"       "1"       "2"       "3"       "4"       "5"       "6"      
    ##      [,8]      [,9]      [,10]     [,11]     [,12]     [,13]     [,14]    
    ## [1,] "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7"
    ## [2,] "7"       "8"       "10"      "12"      "14"      "16"      "18"     
    ##      [,15]     [,16]     [,17]     [,18]     [,19]     [,20]     [,21]    
    ## [1,] "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7"
    ## [2,] "20"      "21"      "22"      "24"      "26"      "30"      "32"     
    ##      [,22]     [,23]     [,24]     [,25]     [,26]    
    ## [1,] "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7" "H1Rosa7"
    ## [2,] "34"      "36"      "38"      "40"      "42"

That was a bit complicated... In plain english, this is what happened. Assign the following to the variable "splitMatrix": loop through the objects stored in splitList, for each one, take elements number 1 through the minimum number of elements in the whole splitlist. It should work for any strsplit output, as long as your sample namings are consistent!

Lastly, we'll take the second row of this matrix, which should be the day at which each sample was taken, and convert it to a numeric type (it is currently characters ("9") instead of numbers (9) )

``` r
tps <- as.numeric(splitMatrix[2,])
colnames(normCounts) <- tps
```

We have the timepoints (in days) captured as "tps"" and also have renamed the column names of normCounts to the timepoints, so we're ready for some plotting!

Now for some exploratory data analysis!!!
-----------------------------------------

To make a heatmap of the spearman (rank) correlations of the samples with each other: **note:** Spearman rank correlation doesn't need to be normalized counts. Be careful though, pearson correlations do need to be normalized. Also, TPM values are scaled by the length of each specific gene, which changes the order. Therefore **you CAN'T compare TPMs to Counts**

``` r
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
spearmanCorrelations <- cor(normCounts,method = "spearman")
#col is the color scheme, trace is dumb so set it to none, cexRow/Col is the size of the text (<1 smaller, >1 larger)
heatmap.2( spearmanCorrelations ,trace = "none",col = cols,cexRow = .7, cexCol = .7)
```

![](2_BasicRNAseq_files/figure-markdown_github/Spearman-1.png)

Or the pearson correlations (with samples not clustered in the heatmap, which you do with Rowv and Colv = F parameters)

``` r
cols <-  colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pearsonCorrelations <- cor(normCounts,method = "pearson")
heatmap.2( pearsonCorrelations ,trace = "none",col = cols,Rowv = F, Colv = F,cexRow = .7, cexCol = .7)
```

![](2_BasicRNAseq_files/figure-markdown_github/Pearson-1.png)

### Plot a list of genes

Now we take an array of genes

``` r
  listOfGenes <- c("PAX6","ASCL1", "TUBB3", "DCX", "POU5F1","NANOG", "GAPDH")
  #Make plots 2 per row, 2 per column
  mypar(2,2)
  for(g in listOfGenes){
    plot(tps,normCounts[g,], main = g,xlab="Day",ylab = "Normalized counts" )
  }
```

![](2_BasicRNAseq_files/figure-markdown_github/ListOGenes-1.png)![](2_BasicRNAseq_files/figure-markdown_github/ListOGenes-2.png)

### Save Scatterplots to a PDF

We can also save this output as a multipage PDF by wrapping the code that creates a plot for each gene inside of the pdf(), which opens a pdf writer as the place your plots will be outputted to, and dev.off() which finalizes and saves the pdf.

``` r
  listOfGenes <- c("PAX6","ASCL1", "TUBB3", "DCX", "POU5F1","NANOG", "GAPDH")
  #Make plots 2 per row, 2 per column
  pdf(file = "~/Desktop/FavoriteGeneScatterplots.pdf", width = 10,height = 10)
  mypar(2,2)
  for(g in listOfGenes){
    plot(tps,normCounts[g,], main = g,xlab="Day",ylab = "Normalized counts" )
  }
  dev.off()
```

##### Quick basic computer tips and vocab review:

**directory** means folder in computer lingo

**command line**: the terminal app on mac (in utilities directory within Applications), a basic guide can be found here: <https://www.tjhsst.edu/~dhyatt/superap/unixcmd.html> . You don't need this unless you want to go deeper.

You can type "pwd ~" (pwd stands for print working directory) on the command line and it will tell you what your home directory is. Then you can write paths like ~/GhostOfMT/Data/FILENAME.txt if I have the GhostOfMT directory in my home folder (the one that contains ~/Downloads/ , ~/Documents/ and ~/Desktop/ etc)

**library**: An R package, a bunch of code and functions wrapped up which can be loaded together with library()

**vector** What you get when you string things together using c()

**list**: A data structure like a vector, get a single element using double brackets \[\[\]\] instead of single brackets \[\]

If the computer can't find a function, you probably haven't loaded the library you need! Paste error messages into google. You can solve 99.999% of your problems by copying and pasting to and from google or stackOverflow