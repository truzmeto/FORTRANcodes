#!/bin/Rscript

# read data
data <- read.table("tmp", sep="", header = FALSE)

#split data frame based on factor into many d.frames(V1 is n_esc)
dframes <- split(data, as.factor(data$V1))

# initilize empty list
my_list <- vector("list", 12)

# number of dframes in a list
N <- length(dframes)

uniq_frame <- 0
dq_thresh <- 3

#loop over list elements(data.frames)
for(i in 1:(N-1)){
    #drop first col from each since its not needed
    dframes[[i]] <- dframes[[i]][-1]
    
    #convert columns to rows by taking a transpose
    my_list[[i]]  <- t(dframes[[i]])
    
    #count total number of configs(frames) for each d.frame  
    colN <- ncol(my_list[[i]])
    
    for(j in 1:colN){    
        a <- colSums(abs(my_list[[i]] - my_list[[i]][,1])) > dq_thresh
        my_list[i] <- my_list[[i]][,a]
        uniq[i] <- uniq[i] + 1      
        
        a <- 0
        if(ncol(my_list[[i]]) <= 1) break
    }
}

#      uniq[j] <- sum(unique(colSums(abs(dframes[[i]] - dframes[[i]][,j]))) > 2) # too nested
#uniq_frame[i] <- sum(uniq)
