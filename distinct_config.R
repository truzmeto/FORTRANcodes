

# read data
data <- read.table("tmp", sep="", header = FALSE)

# split data frame based on factor into many d.frames(V1 is n_esc)
dframes <- split(data, as.factor(data$V1))

my_list <- vector("list", 12)

for(i in 1:length(dframes)){
    # drop first col from each
    dframes[[i]] <- dframes[[i]][-1]
    
    # convert columns to rows by taking a transpose
    dframes[[i]] <- t(dframes[[i]])
    
    colN <- ncol(dframes[[i]])
    
    for(j in 1:(colN-1)){
        for(k in j:colN){
            dframes[[i]][,j] - dframes[[i]][,k]
        }
    }
    
    # do crossproduct between all columns(configuration vectors) and get new matrix
    #my_list[[i]] <- crossprod(data.matrix(dframes[[i]]))
}
