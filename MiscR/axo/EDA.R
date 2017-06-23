library("getopt")
spec <- matrix( c("after",  'a', 1, "character",
                  "before", 'b', 1, "character",
                  "help",   'h', 0, "logical",
                  "out",    'o', 1, "character",
                  "res",    'r', 1, "integer",
                  "height", 't', 1, "integer",
                  "width",  'w', 1, "integer",
                  "ylim",   'y', 1, "double"),
                byrow=TRUE, ncol=4)
opt  <- getopt(spec)
if(is.null(opt$after))  { opt$after  <- "post-norm.png" }
if(is.null(opt$before)) { opt$before <- "pre-norm.png"  }
if(is.null(opt$height)) { opt$height <- 1200            }
if(is.null(opt$res))    { opt$res    <- 150             }
if(is.null(opt$width))  { opt$width  <- 1200            }
if(!is.null(opt$ylim))  { opt$ylim   <- c(0, opt$ylim)  }

# Load data, determine number of samples
data  <- read.table(file("stdin"), header=TRUE, sep="\t", quote="")
nsamp <- dim(data)[2] - 1
data  <- data[,1:nsamp+1]

# Plot distribution of expression values before normalization
png(opt$before, height=opt$height, width=opt$width, res=opt$res)
h <- hist(log(data[,1]), plot=FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="Log expression value", ylab="Proportion of molecules", ylim=opt$ylim)
for(i in 2:nsamp)
{
  h <- hist(log(data[,i]), plot=FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}
devnum <- dev.off()

# Normalize by median
size.factors <- MedianNorm(data.matrix(data))
data.norm <- t(apply(data, 1, function(x){ x / size.factors }))

# Plot distribution of normalized expression values
png(opt$after, height=opt$height, width=opt$width, res=opt$res)
h <- hist(log(data.norm[,1]), plot=FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="Log normalized expression value", ylab="Proportion of molecules",
     ylim=opt$ylim)
for(i in 2:nsamp)
{
  h <- hist(log(data.norm[,i]), plot=FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}



ECEmatch <- ECE[intersect(rownames(ECB),rownames(ECE)),]
ECBmatch <- ECB[intersect(rownames(ECB),rownames(ECE)),]
idx <- rowMeans(ECEmatch)>10 & rowMeans(ECBmatch)>10
ECEmatch <- ECEmatch[idx,]
ECBmatch <- ECBmatch[idx,]
unionnnorm <- cbind(ECEmatch,ECBmatch)
mnunion <-  MedianNorm(unionnnorm)
normmatch <- log(GetNormalizedMat(unionnnorm,mnunion),2)



allAmby002Data <- log(cbind(exprsGSE67118,exprsGSE35255,exprsGSE37198),2)
allAmby002Data <- normalizeBetweenArrays(exprsGSE67118)
nsamp <- dim(allAmby002Data)[2]
h <- hist(allAmby002Data[,1], plot=FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="Log normalized expression value", ylab="Proportion of molecules")
for(i in 2:nsamp){
  h <- hist(allAmby002Data[,i], plot=FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}
pData(GSE67118)

