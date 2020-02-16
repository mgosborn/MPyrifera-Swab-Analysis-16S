#This script utilizes the ShortRead library in R to read and write fastq files. 
#This is highly inefficient for demultiplexing large files and was rewritten in python.
#However, this is useful for demultiplexing smaller files or sampling larger files.
#This script is modified from https://github.com/levisimons/oysterGut/blob/master/demultiplex.py

#Note: "small_MO_1.fq.gz" is a smaller, identical version of MO_1.fq.gz. Same for rev ("_2") files

# Create a random sample (subset) fq file 
##require(ShortRead)
##setwd('/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads')
##wd <- c('/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads')
##path1 <- c('/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads/Sampler_1.fq.gz')
##path2 <- c('/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads/Sampler_2.fq.gz')
##readFastq(wd, "MO_")
##sample <- FastqSampler("MO_1.fq.gz", n = 1000)
##sFq <- yield(sample)
##writeFastq(object=sFq,file = "/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads/Sampler_1.fq.gz", mode = "w")
##sample <- FastqSampler("MO_2.fq.gz", n = 1000)
##sFq <- yield(sample)
##writeFastq(sample,path2, mode = "w")
##writeFastq(object=sFq,file = "/home/cmb-07/sn1/mgosborn/KelpMicrobiome/121319_SequenceReads/Sampler_2.fq.gz", mode = "w")

# Demultiplex 
print("Loading ShortRead library...")
library(ShortRead)
print("Library loaded.")

#Function to extract index from read identifier line
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x)-2)
}
print("Custom function written.")

barcode515F = c('TCAGC','GTATC','GCTAC','ACGCA','CGCGT','CTGGT','GCGTT','GGAAC','AAGCC','AGCTT','CCCTT','CGCAC','GGTGT','GGCAG','TGATA','TGTGC')
# barcode is found in line 2 of 'MO_1.fq.gz' (position 5-9)
index926R = c('ATCACG','CGATGT','TTAGGC','TGACCA','ACAGTG','GCCAAT','CAGATC','ACTTGA','ATGTCA','CCGTCC','GTCCGC','GTGAAA','CACCGG', 'CACGAT','CACTCA','CAGGCG')
# index is found in line 1 of 'MO_1.fq.gz' and 'MO_2.fq.gz' (position varies)
print("Barcode and index libraries uploaded")

leftFileSuffix = '_1.fastq'
rightFileSuffix = '_2.fastq'


setwd('/staging/sn1/mgosborn/fulltest')
print("Set working directory to: /staging/sn1/mgosborn/fulltest")

print("Attempting to read in fwd reads...")
fwd <- readFastq("small_MO_1.fq.gz") #Read in fwd reads. 
print("Fwd reads loaded.")
print("Attempting to read in rev reads...")
rev <- readFastq("small_MO_2.fq.gz") #Read in rev reads.
print("Rev reads loaded.")

print("Attempting to demultiplex...")
for(read in 1:length(rev)){
  
  fwdIdentifier = id(fwd[read])     #Read in the sequence identifiers.
  revIdentifier = id(rev[read])     
  
  
  index = substrRight(id(fwd[read]),8)                    #Extract index.
  barcode = substr(sread(fwd[read]), 5, 9)                #Extract barcode.
  
  if(identical(fwdIdentifier, revIdentifier) == TRUE){     #check if read idendifiers are the same
    if(index %in% index926R == TRUE){                      #check for index
      if(barcode %in% barcode515F == TRUE){                #check for barcode
        
        indexNum <- match(index, index926R)                #match to known index
        barcodeNum <- match(barcode,barcode515F)           #match to known barcode
        myFileName <- paste("F",barcodeNum,"R",indexNum,leftFileSuffix, sep="")       #file name specific to barcode/index pair
        myFilePath <- file.path("/staging","sn1","mgosborn","fulltest",myFileName)            #file path for file name

        if(file.exists(myFileName) == FALSE){
          writeFastq(fwd[read], myFilePath, mode = "w")          #write or append this read to appropriate file (i.e. demultiplex)
        }
        if(file.exists(myFileName) == TRUE){
          writeFastq(fwd[read], myFilePath, mode = "a")
        }
        
        
        myFileName <- paste("F",barcodeNum,"R",indexNum,rightFileSuffix, sep="")       #file name specific to barcode/index pair
        myFilePath <- file.path("/staging","sn1","mgosborn","fulltest",myFileName)            #file path for file name
        
        if(file.exists(myFileName) == FALSE){
          writeFastq(rev[read], myFilePath, mode = "w")           #write or append this read to appropriate file (i.e. demultiplex)
        }
        if(file.exists(myFileName) == TRUE){
          writeFastq(rev[read], myFilePath, mode = "a")
        }        
        
      }
      else{
        next 
        #If identifiers match but barcode not in barcode515F, skip to next iteration. 
      }
    }
    else{
      next
      #If identifiers match byt index not in index926R, skip to next iteration. 
    }
  }
  else{
    rev[read] <- rev[read+1]
    #If fwd and rev identifiers don't match, advance rev read by 1. 
  }
  if(read %% 100000 ==0) {
    print(paste("Still working... read", read, "of", length(rev),"."))
  }
}

print("Boom shakalaka!")

