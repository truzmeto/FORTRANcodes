#!/bin/Rscript

#####################################################################################################
#---------------------------------- T. Ruzmetov, October 19, 2017 ----------------------------------#
#This code calculates number of unique configurations per binding attempt for bimolecular binding.
#It reads

#read data
input_file <- "./tmp"
data <- read.table(input_file, sep="", header = FALSE)

#split data frame based on factor into many d.frames(V1 is n_esc)
dframes <- split(data, as.factor(data$V1))

#number of dframes in a list
N <- length(unique(data$V1))

#initilize empty list
my_list <- vector("list", N)

dq_thresh <- 1 # threshold in q to distinguish configurations
uniq <- 0 

#loop over list elements(data.frames)
for(i in 1:N){
    #drop first col from each since its not needed
    dframes[[i]] <- dframes[[i]][-1]
    
    #convert columns to rows by taking a transpose
    my_list[[i]]  <- t(dframes[[i]])
    
    #count total number of configs(frames) for each d.frame  
    colN <- ncol(my_list[[i]])
    
    #loop over columns, compare each column to all other columns and subset
    #the matrix based on difference. This helps to count # of uniq configurations
    k <- 0  
    for(j in 1:colN){    
        #generating col.indicies for uniq columns
        sub_indx <- colSums(abs(my_list[[i]] - my_list[[i]][,1])) >= dq_thresh
  
        #subsetting matrix based on uniq index
        my_list[[i]] <- my_list[[i]][,sub_indx]
        
        #counting # of uniq cases
        k <- k + 1      
        
        #setting sub index to zero
        sub_indx <- 0
        
        #here I break out if matrix collapses to one column or less
        if(ncol(my_list[[i]]) <= 1) break
    }
    uniq[i] <- k 
}

#print out n_esc, sum of all uniq confs over all attampts, uniq confs for the last attempt

configs <- sum(uniq[-N])
last_att <- uniq[N]
output <- data.frame(N-1, configs, last_att)
write.table(output, file ="", row.names = FALSE, col.names = FALSE, sep = "  ")

