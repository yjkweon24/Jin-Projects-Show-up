#######################################
# filter zeros and singletons this is not good enough ==> single occurance
# reads occur in only one sample should be treated as artifacts and removed


#I made some changes on this function below 
SanityCheckData <- function(datatype){
    
    if(datatype != "omics"){
        current.msg <<-"Not omics.";
        return(0);
    }
    current.msg <<- paste("No missing values for both omics datasets."); 
    dataSet$missingimp.msg <- current.msg;
    dataSet <<- dataSet;

    fircolnm <- colnames(dataSet$sample_data1)[1]; #metadata
    fircol <- dataSet$sample_data1[,1]; #metadata

    fircolnm2 <- colnames(dataSet$sample_data2)[1]; #metadata
    fircol2 <- dataSet$sample_data2[,1]; #metadata

    feat.sums1 <- apply(dataSet$sample_data1[,-1], 2, function(x){sum(x>0, na.rm=T)}); #first column is metatdata
    feat.sums2 <- apply(dataSet$sample_data2[,-1], 2, function(x){sum(x>0, na.rm=T)}); 
    
    gd.inx <- feat.sums1 > 1; # occur in at least 2 samples for the 1st omics
    gd.inx2  <- feat.sums2 > 1; # occur in at least 2 samples for the 2nd omics
    
    #gd.inxnew <- gd.inx & gd.inx2; #check at least one of them is False (omcis sample should always be the same)
    
    
    if(length(which(gd.inx=="TRUE"))==0 | length(which(gd.inx2=="TRUE"))==0){
        current.msg <<-"Reads occur in only one sample.  All these are considered as artifacts and have been removed from data. No data left after such processing.";
        return(0);
    }

    data.proc <- dataSet$sample_data1[,-1][, gd.inx];
    data.proc2 <- dataSet$sample_data2[,-1][, gd.inx2];

    # filtering the constant features here
    # check for columns with all constant (var=0)
    varCol <- apply(data.proc, 2, var, na.rm=T);
    constCol <- varCol == 0 | is.na(varCol);

    varCol2 <- apply(data.proc2, 2, var, na.rm=T);
    constCol2 <- varCol2 == 0 | is.na(varCol2);

    #constColnew <- constCol & constCol2; #check at least one of them is False (omcis sample should always be the same)

    data.proc <- data.proc[, !constCol];
    data.proc2 <- data.proc2[, !constCol2];
    if(length(data.proc)==0 | length(data.proc2)==0){
        current.msg <<-"All features are found to be constant and have been removed from data. No data left after such processing.";
        return(0);
    }

    data.proc <- cbind(factor(fircol), as.data.frame(data.proc));
    colnames(data.proc)[1] <- fircolnm;
    
    data.proc2 <- cbind(factor(fircol2), as.data.frame(data.proc2));
    colnames(data.proc2)[1] <- fircolnm2;
    

    saveRDS(data.proc, file="data.proc.orig"); # save an copy
    saveRDS(data.proc, file="data.prefilt"); # save an copy

    saveRDS(data.proc2, file="data.proc.orig2"); # save an copy for the 2nd omics
    saveRDS(data.proc2, file="data.prefilt2"); # save an copy for the 2nd omics

    saveRDS(dataSet$sample_data1, file = "data.sample_data1") 
    saveRDS(dataSet$sample_data2, file = "data.sample_data2") 

    
    
    dataSet$proc <- data.proc;
    dataSet$proc2 <- data.proc2;

    dataSet<<-dataSet; 

    sample_no <- nrow(dataSet$proc);
    sample_no2 <- nrow(dataSet$proc2);

    vari_no <- ncol(dataSet$proc)-1;
    vari_no2 <- ncol(dataSet$proc2)-1;

    smpl.sums <- apply(dataSet$proc[,-1], 1, sum);
    smpl.sums2 <- apply(dataSet$proc2[,-1], 1, sum);

    tot_size<-sum(smpl.sums);
    tot_size2<-sum(smpl.sums2);
    
    smean <- mean(smpl.sums);
    smean2 <- mean(smpl.sums2);

    smin <- min(smpl.sums);
    smin2 <- min(smpl.sums2);

    smax <- max(smpl.sums); 
    smax2 <- max(smpl.sums2); 

    if(identical(rownames(dataSet$proc), rownames(dataSet$proc2))){
        samanme_same <- 1
    }else{
        samanme_same <- 0
    }
    save.image("sanity.RData")
    return(c(1,sample_no, sample_no2, vari_no, vari_no2, tot_size, tot_size2, smean, smean2, smin, smin2, smax, smax2, samanme_same));
}


#This is when the data is already normalized and filtered well -> For example datasets - so no sanity check needed 
CheckData <- function(datatype){
    
    if(datatype != "omics"){
        current.msg <<-"Not omics.";
        return(0);
    }

    current.msg <<- paste("No missing values for both omics datasets."); 
    dataSet$missingimp.msg <- current.msg;
    dataSet <<- dataSet;

    data.proc <- dataSet$sample_data1;
    data.proc2 <- dataSet$sample_data2;

    saveRDS(data.proc, file="data.proc.orig"); # save an copy
    saveRDS(data.proc, file="data.prefilt"); # save an copy

    saveRDS(data.proc2, file="data.proc.orig2"); # save an copy for the 2nd omics
    saveRDS(data.proc2, file="data.prefilt2"); # save an copy for the 2nd omics

    saveRDS(dataSet$sample_data1, file = "data.sample_data1") 
    saveRDS(dataSet$sample_data2, file = "data.sample_data2") 

    

    dataSet$proc <- data.proc;
    dataSet$proc2 <- data.proc2;

    dataSet<<-dataSet; 

    sample_no <- nrow(dataSet$proc);
    sample_no2 <- nrow(dataSet$proc2);

    vari_no <- ncol(dataSet$proc)-1;
    vari_no2 <- ncol(dataSet$proc2)-1;

    smpl.sums <- apply(dataSet$proc[,-1], 1, sum);
    smpl.sums2 <- apply(dataSet$proc2[,-1], 1, sum);

    tot_size<-sum(smpl.sums);
    tot_size2<-sum(smpl.sums2);
    
    smean <- mean(smpl.sums);
    smean2 <- mean(smpl.sums2);

    smin <- min(smpl.sums);
    smin2 <- min(smpl.sums2);

    smax <- max(smpl.sums); 
    smax2 <- max(smpl.sums2); 

    if(identical(rownames(dataSet$proc), rownames(dataSet$proc2))){
        samanme_same <- 1
    }else{
        samanme_same <- 0
    }
    #save.image("sanity.RData")
    return(c(1,sample_no, sample_no2, vari_no, vari_no2, tot_size, tot_size2, smean, smean2, smin, smin2, smax, smax2, samanme_same));
}



############################



#Missing value check
MissingCheckData <- function(datatype){
    
    data1 <- dataSet$sample_data1; 
    data2 <- dataSet$sample_data2; 
    msg <- NULL;

    missingcheck <- c();
    miss_check <- function(data1, data2){
        if(sum(colSums(is.na(data1))) != 0 & sum(colSums(is.na(data2))) != 0){
          x <- 3 #both data have missing values
        }else{
          if(sum(colSums(is.na(data1))) == 0){
            if(sum(colSums(is.na(data2))) != 0){
              x <- 2 #only data2 has the missing values.
            }else{
              x <- 0 #no missing value...
            } 
          } else{
          x <- 1 #only data1 has the missing values
          }     
        }
        return(x)
    }

    ans <- miss_check(data1[,-1], data2[,-1]);

    saveRDS(data1, file="data.sample_data1"); # save an copy
    saveRDS(data2, file="data.sample_data2"); # save an copy for the 2nd omics

    #save.image("missingcheck.RData")
    return(ans);
}

#Missing value imputation omics 1  
ApplyMissing1 <- function(missImp){
    
    data <- readRDS("data.sample_data1"); 
    data2 <- readRDS("data.sample_data2"); 
    msg <- NULL;

    #this data is used for sample categorial comparision further
    rmn_feat <- ncol(data);
    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata

    #remove samples with 100% missing missing/missing rows (missing for every variable)
    NAs_in <- which(apply(data[,-1], 1, function(x) sum(is.na(x))/length(x))==1)
    if(length(NAs_in) > 0){
        data1_new <- data[, -1][-(NAs_in),];
        data2_new <- data2[-(NAs_in),];
    } else{
        data1_new <- data[,-1];
        data2_new <- data2;
    }


    if(missImp == "knn"){
        suppressMessages(require("DMwR"));
        data1_new <- knnImputation(data1_new);
        current.msg <<- paste("Omics1: Missing value imputed with KNN."); 
    }else if(missImp == "min"){
        suppressMessages(require("Hmisc"));
        
        data1_new <- apply(data1_new, 2, function(x) impute(x, fun = min))
        current.msg <<- paste("Omics1: Missing value imputed with the feature's minimum value."); 
    }else if(missImp == "med"){
        suppressMessages(require("Hmisc"));

        data1_new <- apply(data1_new, 2, impute)
        current.msg <<- paste("Omics1: Missing value imputed with the feature's median value."); 
    }else if(missImp == "mean"){
        suppressMessages(require("Hmisc"));

        data1_new <- apply(data1_new, 2, function(x) impute(x, fun = mean))
        current.msg <<- paste("Omics1: Missing value imputed with the feature's mean value."); 
    }else if(missImp == "min2"){
        suppressMessages(require("Hmisc"));

        min2func <- function(x){
          na_location <- unname(which(is.na(x)))
          x <- impute(x/2, fun = min)
          if (length(na_location) == 0){
            x <- 2*x
          }else{
            x[-na_location] <- 2*x[-na_location]
          }
          return(x)
        }

        data1_new <- apply(data1_new, 2, min2func)
        current.msg <<- paste("Omics1: Missing value imputed with the feature's minimum/2 value."); 
    }


    data <- cbind(factor(fircol), as.data.frame(data1_new));
    colnames(data)[1] <- fircolnm;
    
    data2 <- data2_new;
    
    data <- as.data.frame(data);
    data2 <- as.data.frame(data2);

    saveRDS(data, file="data.prefilt"); # save an copy
    saveRDS(data2, file="data.prefilt2"); # save an copy for the 2nd omics
    
    dataSet$sample_data1 <- data;
    dataSet <<- dataSet;


    dataSet$missingimp.msg <- current.msg;
    dataSet <<- dataSet;
    #save.image("missingcheck.RData")
    return(1);
}



#Missing value imputation omics 2  
ApplyMissing2 <- function(missImp){

    data <- readRDS("data.sample_data1"); 
    data2 <- readRDS("data.sample_data2"); 
    msg <- NULL;


    fircolnm <- colnames(data2)[1]; #metadata
    fircol <- data2[,1]; #metadata

    #remove samples with 100% missing missing/missing rows (missing for every variable)
    NAs_in <- which(apply(data2[,-1], 1, function(x) sum(is.na(x))/length(x))==1)
    if(length(NAs_in) > 0){
        data1_new <- data[-(NAs_in),];
        data2_new <- data2[, -1][-(NAs_in),];
    } else{
        data1_new <- data;
        data2_new <- data2[,-1];
    }


    if(missImp == "knn"){
        suppressMessages(require("DMwR"));
        data2_new <- knnImputation(data2_new);
        current.msg <<- paste("Omics2: Missing value imputed with KNN."); 
    }else if(missImp == "min"){
        suppressMessages(require("Hmisc"));
        
        data2_new <- apply(data2_new, 2, function(x) impute(x, fun = min))
        current.msg <<- paste("Omics2: Missing value imputed with the feature's minimum value."); 
    }else if(missImp == "med"){
        suppressMessages(require("Hmisc"));

        data2_new <- apply(data2_new, 2, impute)
        current.msg <<- paste("Omics2: Missing value imputed with the feature's median value."); 
    }else if(missImp == "mean"){
        suppressMessages(require("Hmisc"));

        data2_new <- apply(data2_new, 2, function(x) impute(x, fun = mean))
        current.msg <<- paste("Omics2: Missing value imputed with the feature's mean value."); 
    }else if(missImp == "min2"){
        suppressMessages(require("Hmisc"));
        
        min2func <- function(x){
          na_location <- unname(which(is.na(x)))
          x <- impute(x/2, fun = min)
          if (length(na_location) == 0){
            x <- 2*x
          }else{
            x[-na_location] <- 2*x[-na_location]
          }
          return(x)
        }

        data2_new <- apply(data2_new, 2, min2func)
        current.msg <<- paste("Omics2: Missing value imputed with the feature's minimum/2 value."); 
    }


    data2 <- cbind(factor(fircol), as.data.frame(data2_new));
    colnames(data2)[1] <- fircolnm;
    
    data <- data1_new;

    data <- as.data.frame(data);
    data2 <- as.data.frame(data2);

    saveRDS(data, file="data.prefilt"); # save an copy
    saveRDS(data2, file="data.prefilt2"); # save an copy for the 2nd omics
    
    dataSet$sample_data2 <- data2;
    dataSet <<- dataSet;
    
    dataSet$missingimp.msg <- current.msg;
    dataSet <<- dataSet;
    #save.image("missingcheck.RData")
    return(1);
}


#Missing value imputation both omics 1 and 2  
ApplyMissing3 <- function(missImp, missImp2){

    data <- readRDS("data.sample_data1"); 
    data2 <- readRDS("data.sample_data2"); 
    msg <- NULL;


    fircolnm <- colnames(data2)[1]; #metadata
    fircol <- data2[,1]; #metadata

    #remove samples with 100% missing missing/missing rows (missing for every variable)
    NAs_in1 <- which(apply(data[,-1], 1, function(x) sum(is.na(x))/length(x))==1)
    NAs_in2 <- which(apply(data2[,-1], 1, function(x) sum(is.na(x))/length(x))==1)

    unions <- union(NAs_in1, NAs_in2)
    if(length(unions) > 0){
      data1_new <- data[,-1][-(unions),]
      data2_new <- data2[,-1][-(unions),]
    } else{
        data1_new <- data[,-1];
        data2_new <- data2[,-1];
    }


    #For omics 1
    if(missImp == "knn"){
        suppressMessages(require("DMwR"));
        data1_new <- knnImputation(data1_new);
        current.msg <<- paste("Omics1: Missing value imputed with KNN."); 
    }else if(missImp == "min"){
        suppressMessages(require("Hmisc"));
        
        data1_new <- apply(data1_new, 2, function(x) impute(x, fun = min))
        current.msg <<- paste("Omics1: Missing value imputed with the feature's minimum value."); 
    }else if(missImp == "med"){
        suppressMessages(require("Hmisc"));

        data1_new <- apply(data1_new, 2, impute)
        current.msg <<- paste("Omics1: Missing value imputed with the feature's median value."); 
    }else if(missImp == "mean"){
        suppressMessages(require("Hmisc"));

        data1_new <- apply(data1_new, 2, function(x) impute(x, fun = mean))
        current.msg <<- paste("Omics1: Missing value imputed with the feature's mean value."); 
    }else if(missImp == "min2"){
        suppressMessages(require("Hmisc"));

        min2func <- function(x){
          na_location <- unname(which(is.na(x)))
          x <- impute(x/2, fun = min)
          if (length(na_location) == 0){
            x <- 2*x
          }else{
            x[-na_location] <- 2*x[-na_location]
          }
          return(x)
        }    

        data1_new <- apply(data1_new, 2, min2func)
        current.msg <<- paste("Omics1: Missing value imputed with the feature's minimum/2 value."); 
    }


    #For omics 2
    if(missImp2 == "knn"){
        suppressMessages(require("DMwR"));
        data2_new <- knnImputation(data2_new);
        current.msg <<- c(current.msg, paste(" Omics2: Missing value imputed with KNN.")); 
    }else if(missImp2 == "min"){
        suppressMessages(require("Hmisc"));
        
        data2_new <- apply(data2_new, 2, function(x) impute(x, fun = min))
        current.msg <<- c(current.msg, paste(" Omics2: Missing value imputed with the feature's minimum value."));
    }else if(missImp2 == "med"){
        suppressMessages(require("Hmisc"));

        data2_new <- apply(data2_new, 2, impute)
        current.msg <<- c(current.msg, paste(" Omics2: Missing value imputed with the feature's median value."));
    }else if(missImp2 == "mean"){
        suppressMessages(require("Hmisc"));

        data2_new <- apply(data2_new, 2, function(x) impute(x, fun = mean))
        current.msg <<- c(current.msg, paste(" Omics2: Missing value imputed with the feature's mean value."));
    }else if(missImp2 == "min"){
        suppressMessages(require("Hmisc"));
        
        min2func <- function(x){
          na_location <- unname(which(is.na(x)))
          x <- impute(x/2, fun = min)
          if (length(na_location) == 0){
            x <- 2*x
          }else{
            x[-na_location] <- 2*x[-na_location]
          }
          return(x)
        }

        data2_new <- apply(data2_new, 2, min2func)
        current.msg <<- c(current.msg, paste(" Omics2: Missing value imputed with the feature's minimum/2 value."));
    }


    data <- cbind(factor(fircol), as.data.frame(data1_new));
    colnames(data)[1] <- fircolnm;

    data2 <- cbind(factor(fircol), as.data.frame(data2_new));
    colnames(data2)[1] <- fircolnm;
    
    data <- as.data.frame(data);
    data2 <- as.data.frame(data2);

    saveRDS(data, file="data.prefilt"); # save an copy
    saveRDS(data2, file="data.prefilt2"); # save an copy for the 2nd omics
    
    dataSet$sample_data1 <- data;
    dataSet$sample_data2 <- data2;
    dataSet <<- dataSet;
     
    dataSet$missingimp.msg <- current.msg;
    dataSet <<- dataSet;
    save.image("missingcheck.RData")
    return(1);
}



SanityCheckData3 <- function(datatype){
    
    data <- readRDS("data.prefilt"); 
    data2 <- readRDS("data.prefilt2"); 

    dataSet$sample_data1 <- data;
    dataSet$sample_data2 <- data2;

    if(datatype != "omics"){
        current.msg <<-"Not omics.";
        return(0);
    }

    fircolnm <- colnames(dataSet$sample_data1)[1]; #metadata
    fircol <- dataSet$sample_data1[,1]; #metadata

    fircolnm2 <- colnames(dataSet$sample_data2)[1]; #metadata
    fircol2 <- dataSet$sample_data2[,1]; #metadata

    feat.sums1 <- apply(dataSet$sample_data1[,-1], 2, function(x){sum(x>0, na.rm=T)}); #first column is metatdata
    feat.sums2 <- apply(dataSet$sample_data2[,-1], 2, function(x){sum(x>0, na.rm=T)}); 
    
    gd.inx <- feat.sums1 > 1; # occur in at least 2 samples for the 1st omics
    gd.inx2  <- feat.sums2 > 1; # occur in at least 2 samples for the 2nd omics
    
    #gd.inxnew <- gd.inx & gd.inx2; #check at least one of them is False (omcis sample should always be the same)
    
    
    if(length(which(gd.inx=="TRUE"))==0 | length(which(gd.inx2=="TRUE"))==0){
        current.msg <<-"Reads occur in only one sample.  All these are considered as artifacts and have been removed from data. No data left after such processing.";
        return(0);
    }

    data.proc <- dataSet$sample_data1[,-1][, gd.inx];
    data.proc2 <- dataSet$sample_data2[,-1][, gd.inx2];

    # filtering the constant features here
    # check for columns with all constant (var=0)
    varCol <- apply(data.proc, 2, var, na.rm=T);
    constCol <- varCol == 0 | is.na(varCol);

    varCol2 <- apply(data.proc2, 2, var, na.rm=T);
    constCol2 <- varCol2 == 0 | is.na(varCol2);

    #constColnew <- constCol & constCol2; #check at least one of them is False (omcis sample should always be the same)

    # making copy of data.proc and proc.phyobj(phyloseq)
    data.proc <- data.proc[, !constCol];
    data.proc2 <- data.proc2[, !constCol2];
    if(length(data.proc)==0 | length(data.proc2)==0){
        current.msg <<-"All features are found to be constant and have been removed from data. No data left after such processing.";
        return(0);
    }

    data.proc <- cbind(factor(fircol), as.data.frame(data.proc));
    colnames(data.proc)[1] <- fircolnm;
    
    data.proc2 <- cbind(factor(fircol2), as.data.frame(data.proc2));
    colnames(data.proc2)[1] <- fircolnm2;
    

    saveRDS(data.proc, file="data.proc.orig"); # save an copy
    saveRDS(data.proc, file="data.prefilt"); # save an copy

    saveRDS(data.proc2, file="data.proc.orig2"); # save an copy for the 2nd omics
    saveRDS(data.proc2, file="data.prefilt2"); # save an copy for the 2nd omics

    saveRDS(dataSet$sample_data1, file = "data.sample_data1") 
    saveRDS(dataSet$sample_data2, file = "data.sample_data2") 

    

    dataSet$proc <- data.proc;
    dataSet$proc2 <- data.proc2;

    dataSet<<-dataSet; 

    sample_no <- nrow(dataSet$proc);
    sample_no2 <- nrow(dataSet$proc2);

    vari_no <- ncol(dataSet$proc)-1;
    vari_no2 <- ncol(dataSet$proc2)-1;

    smpl.sums <- apply(dataSet$proc[,-1], 1, sum);
    smpl.sums2 <- apply(dataSet$proc2[,-1], 1, sum);

    tot_size<-sum(smpl.sums);
    tot_size2<-sum(smpl.sums2);
    
    smean <- mean(smpl.sums);
    smean2 <- mean(smpl.sums2);

    smin <- min(smpl.sums);
    smin2 <- min(smpl.sums2);

    smax <- max(smpl.sums); 
    smax2 <- max(smpl.sums2); 

    if(identical(rownames(dataSet$proc), rownames(dataSet$proc2))){
        samanme_same <- 1
    }else{
        samanme_same <- 0
    }
    #save.image("sanity.RData")
    return(c(1,sample_no, sample_no2, vari_no, vari_no2, tot_size, tot_size2, smean, smean2, smin, smin2, smax, smax2, samanme_same));
}







############################



# filter data based on low counts in high percentage samples
# note, first is abundance, followed by variance
# This is the for the first omics

ApplyAbundanceFilter <- function(filt.opt, count, smpl.perc){

    data <- readRDS("data.prefilt"); 
    
    msg <- NULL;
    #this data is used for sample categorial comparision further
    rmn_feat <- ncol(data);
    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata

    dataSet$metadata <- factor(fircol);
    dataSet$firstcolname <- fircolnm;
    
    if(count==0){# no low-count filtering
        rmn_feat <- ncol(data)-1;
        kept.inx <- rep(TRUE, rmn_feat);
    }else{
        if(filt.opt == "prevalence"){
            rmn_feat <- ncol(data)-1;
            minLen <- smpl.perc*nrow(data);
            kept.inx <- apply(data[,-1], MARGIN = 2,function(x) {sum(x >= count) >= minLen});  
        }else if (filt.opt == "mean"){
            filter.val <- apply(data[,-1], 2, mean, na.rm=T);
            kept.inx <- filter.val >= count;
        }else if (filt.opt == "median"){
            filter.val <- apply(data[,-1], 2, median, na.rm=T);
            kept.inx <- filter.val >= count;
        }else if (filt.opt == "sum"){
            filter.val <- apply(data[,-1], 2, sum, na.rm=T);
            kept.inx <- filter.val >= count;
        }
    }
    
    data <- cbind(factor(fircol), as.data.frame(data[,-1][, kept.inx]));
    colnames(data)[1] <- fircolnm;

    
    dataSet$filt.data <- data;
    saveRDS(dataSet$filt.data, file="filt.data.orig"); # save an copy
    dataSet <<- dataSet;
    current.msg <<- paste("Omics1: A total of ", sum(!kept.inx), " low abundance features were removed based on ", filt.opt, ".", sep=""); 
    dataSet <<- dataSet;
    return(1);
}

# filter data based on low abundace or variance
# note this is applied after abundance filter
# This is the for the first omics

ApplyVarianceFilter <- function(filtopt, filtPerct){

    data <- dataSet$filt.data;
    msg <- NULL;
    
    rmn_feat <- ncol(data);
    filter.val <- nm <- NULL;
    
    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata

    if(filtPerct==0){# no low-count filtering
        rmn_feat <- ncol(data)-1;
        remain <- rep(TRUE, rmn_feat);
    }else{
        if (filtopt == "iqr"){
            filter.val <- apply(data[,-1], 2, IQR, na.rm=T);
            nm <- "IQR";
        }else if (filtopt == "sd"){
            filter.val <- apply(data[,-1], 2, sd, na.rm=T);
            nm <- "standard deviation";
        }else if (filtopt == "cov"){
            sds <- apply(data[,-1], 2, sd, na.rm=T);
            mns <- apply(data[,-1], 2, mean, na.rm=T);
            filter.val <- abs(sds/mns);
            nm <- "Coeffecient of variation";
        }else if (filtopt == "abs"){
          filter.val <- apply(data[,-1], 2, function(x) {max(abs(x))});
          nm <- "Maximum of absolute"
        }
        # get the rank
        rk <- rank(-filter.val, ties.method='random');
        var.num <- ncol(data)-1;
        remain <- rk < var.num*(1-filtPerct);
    }
    
    data <- cbind(factor(fircol), as.data.frame(data[,-1][, remain]));
    colnames(data)[1] <- fircolnm;
        
    dataSet$filt.data <- data;
    saveRDS(dataSet$filt.data, file="filt.data.orig"); # save an copy
    rm.msg1 <- paste("Omics1: A total of ", sum(!remain), " low variance features were removed based on ", filtopt, ".", sep="");
    rm.msg2 <- paste("The number of features remains after the data filtering step:", ncol(data));   
    current.msg <<- paste(c(current.msg,rm.msg1,rm.msg2), collapse=" ");
    dataSet$filt.msg <- current.msg;
    dataSet <<- dataSet;
    #save.image("filter1.RData")
    return(1);
}



# filter data based on low counts in high percentage samples
# note, first is abundance, followed by variance
# This is the for the second omics

ApplyAbundanceFilter2 <- function(filt.opt, count, smpl.perc){

    data <- readRDS("data.prefilt2"); 
    msg <- NULL;
    #this data is used for sample categorial comparision further
    rmn_feat <- ncol(data);
    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata
    
    if(count==0){# no low-count filtering
        rmn_feat <- ncol(data)-1;
        kept.inx <- rep(TRUE, rmn_feat);
    }else{
        if(filt.opt == "prevalence"){
            rmn_feat <- ncol(data)-1;
            minLen <- smpl.perc*nrow(data);
            kept.inx <- apply(data[,-1], MARGIN = 2,function(x) {sum(x >= count) >= minLen});  
        }else if (filt.opt == "mean"){
            filter.val <- apply(data[,-1], 2, mean, na.rm=T);
            kept.inx <- filter.val >= count;
        }else if (filt.opt == "median"){
            filter.val <- apply(data[,-1], 2, median, na.rm=T);
            kept.inx <- filter.val >= count;
        }else if (filt.opt == "sum"){
            filter.val <- apply(data[,-1], 2, sum, na.rm=T);
            kept.inx <- filter.val >= count;
        }
    }
    data <- cbind(factor(fircol), as.data.frame(data[,-1][, kept.inx]));
    colnames(data)[1] <- fircolnm;
    
    dataSet$filt.data2 <- data;
    saveRDS(dataSet$filt.data2, file="filt.data.orig2"); # save an copy
    dataSet <<- dataSet;
    current.msg <<- paste("Omics2: A total of ", sum(!kept.inx), " low abundance features were removed based on ", filt.opt, ".", sep=""); 
    dataSet <<- dataSet;
    return(1);
}


# filter data based on low abundace or variance
# note this is applied after abundance filter
# This is the for the second omics

ApplyVarianceFilter2 <- function(filtopt, filtPerct){

    data <- dataSet$filt.data2;
    msg <- NULL;
    
    rmn_feat <- ncol(data);
    filter.val <- nm <- NULL;
    
    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata

    if(filtPerct==0){# no low-count filtering
        rmn_feat <- ncol(data)-1;
        remain <- rep(TRUE, rmn_feat);
    }else{
        if (filtopt == "iqr"){
            filter.val <- apply(data[,-1], 2, IQR, na.rm=T);
            nm <- "IQR";
        }else if (filtopt == "sd"){
            filter.val <- apply(data[,-1], 2, sd, na.rm=T);
            nm <- "standard deviation";
        }else if (filtopt == "cov"){
            sds <- apply(data[,-1], 2, sd, na.rm=T);
            mns <- apply(data[,-1], 2, mean, na.rm=T);
            filter.val <- abs(sds/mns);
            nm <- "Coeffecient of variation";
        }else if (filtopt == "abs"){
          filter.val <- apply(data[,-1], 2, function(x) {max(abs(x))});
          nm <- "Maximum of absolute"
        }
        # get the rank
        rk <- rank(-filter.val, ties.method='random');
        var.num <- ncol(data)-1;
        remain <- rk < var.num*(1-filtPerct);
    }
        
    data <- cbind(factor(fircol), as.data.frame(data[,-1][, remain]));
    colnames(data)[1] <- fircolnm;
    
    dataSet$filt.data2 <- data;
    saveRDS(dataSet$filt.data2, file="filt.data.orig2"); # save an copy
    rm.msg1 <- paste("Omics2: A total of ", sum(!remain), " low variance features were removed based on ", filtopt, ".", sep="");
    rm.msg2 <- paste("The number of features remains after the data filtering step:", ncol(data));   
    current.msg <<- paste(c(current.msg,rm.msg1,rm.msg2), collapse=" ");
    dataSet$filt.msg2 <- current.msg;
    dataSet <<- dataSet;
    #save.image("filter2.RData")
    return(1);
}











###################


# note, here also update data type array/count
# This is for the first omics.
PerformNormalization <- function(scale.opt,transform.opt, lambda1, lambda1b, lambda2b){
    
    data <- readRDS("filt.data.orig"); 
    msg <- NULL;

    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata
    lambda1 <- as.double(lambda1);
    lambda1b <- as.double(lambda1b);
    lambda2b <- as.double(lambda2b);  
  
    if(scale.opt != "none"){
        if(scale.opt=="rowsum"){    
            data <- sweep(data[,-1], 1, rowSums(data[,-1]), FUN="/")
            data <- data*10000000;  
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics1: Performed total sum normalization."));
        }else if(scale.opt=="upperquartile"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuUQ <- edgeRnorm(t(data[,-1]),method="upperquartile");
            data <- cbind(fircol, as.data.frame(t(otuUQ$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics1: Performed upper quartile normalization."));
        }else if(scale.opt=="CSS"){
            suppressMessages(require("metagenomeSeq"));
            #has to be in class(matrix only not in phyloseq:otu_table)
            data1 <- as(data[,-1],"matrix");
            dataMR <- newMRexperiment(t(data1));
            data <- cumNorm(dataMR,p=cumNormStat(dataMR));
            data <- MRcounts(data,norm = T);
            data <- cbind(fircol, as.data.frame(t(data)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics1: Performed cumulative sum scaling normalization."));
        }else if(scale.opt=="as"){
            data <- scale(data[,-1])
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics1: Performed auto scaling (centered and standardized)."));
        }else if(scale.opt=="mms"){
            min_max <- function(x){
              return ((x - min(x)) / (max(x) - min(x)));
            }
            data <- apply(data[,-1], 2, min_max);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics1: Performed min-max scaling."));
        }else if(scale.opt=="qn"){
            suppressMessages(require("preprocessCore"));
            data <- normalize.quantiles(t(data[,-1]), copy=TRUE);
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics1: Performed quantile normalization."));
        }else{
            print(paste("Unknown scaling parameter:", scale.opt));
        }
    }else{
        msg <- c(msg, paste("Omics1: No data scaling was performed."));
    }

  

   if(transform.opt != "none"){
        if(transform.opt=="rle"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuRLE <- edgeRnorm(t(data[,-1]),method="RLE");
            data <- cbind(fircol, as.data.frame(t(otuRLE$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed RLE Normalization."));
        }else if(transform.opt=="TMM"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuTMM <- edgeRnorm(t(data[,-1]),method="TMM");
            data <- cbind(fircol, as.data.frame(t(otuTMM$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed TMM Normalization."));
        }else if(transform.opt=="clr"){ 
            data <- apply(t(data[,-1]), 2, clr_transform);
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed centered-log-ratio normalization.")); 
        }else if(transform.opt=="glt"){
            suppressMessages(require("LMGene"));
            data <- glog(data[,-1], lambda1);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed generalized log transformation, with ", lambda1, ".", sep="")); 
        }else if(transform.opt=="tpbct"){
            if(all(data[,-1] > -lambda2b)){
              if(lambda1b == 0){
                data <- log(data[,-1]+lambda2b)
              }else{
                data <- ((data[,-1]+lambda2b)^(lambda1b) - 1)/lambda1b
              }
            }else{
              stop(paste0("some of the data <= -", lambda2))
            }
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed two-parameter box-cox transformation, with ", lambda1b, " and ", lambda2b, ".", sep="")); 
        }else if(transform.opt=="vst"){
            suppressMessages(require("limma"));
            data <- normalizeVSN(t(data[,-1]));
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed variance stabilizing transformation.")); 
        }else if(transform.opt=="l2t"){
            min.val <- min(data[,-1][data[,-1]>0], na.rm=T)/10;
            data[,-1][data[,-1]<=0] <- min.val;
            data <- log2(data[,-1]);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed log2 transformation."));
        }else if(transform.opt=="l2cpm"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            nf <- edgeR::calcNormFactors(as(t(data[,-1]), "matrix"));
            y <- voom(t(data[,-1]),plot=F,lib.size=colSums(t(data))*nf);
            data <- y$E;
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed log2-counts per million."));
        }else{
            print(paste("Unknown scaling parameter:", transform.opt));
        }
    }else{
        msg <- c(msg, paste("No data transformation was performed."));  
    }

    metadata <- dataSet$metadata;
    metadataname <- dataSet$firstcolname;;

    data <- cbind(metadata, data[,-1]);
    colnames(data)[1] <- metadataname;

    #using this object for plotting
    dataSet$norm1 <- data; 
    current.msg <<- paste(msg, collapse=" ");
    dataSet$norm.msg <- current.msg;
    dataSet <<- dataSet;
    #save.image("norm1.RData")
    saveRDS(dataSet$norm1, file="norm.data1"); # save an copy
    return(1);
}





# note, here also update data type array/count
# This is for the second omics.
PerformNormalization2 <- function(scale.opt,transform.opt, lambda1, lambda1b, lambda2b){
    
    data <- readRDS("filt.data.orig2"); 
    msg <- NULL;

    fircolnm <- colnames(data)[1]; #metadata
    fircol <- data[,1]; #metadata
    lambda1 <- as.double(lambda1);
    lambda1b <- as.double(lambda1b);
    lambda2b <- as.double(lambda2b);  
  
    if(scale.opt != "none"){
        if(scale.opt=="rowsum"){    
            data <- sweep(data[,-1], 1, rowSums(data[,-1]), FUN="/")
            data <- data*10000000;  
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics2: Performed total sum normalization."));
        }else if(scale.opt=="upperquartile"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuUQ <- edgeRnorm(t(data[,-1]),method="upperquartile");
            data <- cbind(fircol, as.data.frame(t(otuUQ$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics2: Performed upper quartile normalization."));
        }else if(scale.opt=="CSS"){
            suppressMessages(require("metagenomeSeq"));
            #has to be in class(matrix only not in phyloseq:otu_table)
            data1 <- as(data[,-1],"matrix");
            dataMR <- newMRexperiment(t(data1));
            data <- cumNorm(dataMR,p=cumNormStat(dataMR));
            data <- MRcounts(data,norm = T);
            data <- cbind(fircol, as.data.frame(t(data)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Omics2: Performed cumulative sum scaling normalization."));
        }else if(scale.opt=="as"){
            data <- scale(data[,-1])
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics2: Performed auto scaling (centered and standardized)."));
        }else if(scale.opt=="mms"){
            min_max <- function(x){
              return ((x - min(x)) / (max(x) - min(x)));
            }
            data <- apply(data[,-1], 2, min_max);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics2: Performed min-max scaling."));
        }else if(scale.opt=="qn"){
            suppressMessages(require("preprocessCore"));
            data <- normalize.quantiles(t(data[,-1]), copy=TRUE);
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Omics2: Performed quantile normalization."));
        }else{
            print(paste("Unknown scaling parameter:", scale.opt));
        }
    }else{
        msg <- c(msg, paste("Omics2: No data scaling was performed."));
    }

  

   if(transform.opt != "none"){
        if(transform.opt=="rle"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuRLE <- edgeRnorm(t(data[,-1]),method="RLE");
            data <- cbind(fircol, as.data.frame(t(otuRLE$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed RLE Normalization."));
        }else if(transform.opt=="TMM"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            otuTMM <- edgeRnorm(t(data[,-1]),method="TMM");
            data <- cbind(fircol, as.data.frame(t(otuTMM$counts)));
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed TMM Normalization."));
        }else if(transform.opt=="clr"){ 
            data <- apply(t(data[,-1]), 2, clr_transform);
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed centered-log-ratio normalization."));
        }else if(transform.opt=="glt"){
            suppressMessages(require("LMGene"));
            data <- glog(data[,-1], lambda1);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed generalized log transformation, with ", lambda1, ".", sep=""));
        }else if(transform.opt=="tpbct"){
            if(all(data[,-1] > -lambda2b)){
              if(lambda1b == 0){
                data <- log(data[,-1]+lambda2b)
              }else{
                data <- ((data[,-1]+lambda2b)^(lambda1b) - 1)/lambda1b
              }
            }else{
              stop(paste0("some of the data <= -", lambda2))
            }
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed two-parameter box-cox transformation, with ", lambda1b, " and ", lambda2b, ".", sep=""));
        }else if(transform.opt=="vst"){
            suppressMessages(require("limma"));
            data <- normalizeVSN(t(data[,-1]));
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed variance stabilizing transformation."));
        }else if(transform.opt=="l2t"){
            min.val <- min(data[,-1][data[,-1]>0], na.rm=T)/10;
            data[,-1][data[,-1]<=0] <- min.val;
            data <- log2(data[,-1]);
            data <- cbind(fircol, data);
            colnames(data)[1] <- fircolnm;
            msg <- c(msg, paste("Performed log2 transformation."));
        }else if(transform.opt=="l2cpm"){
            suppressMessages(require("edgeR"));
            data[,-1][data[,-1]<0] <- 0
            nf <- edgeR::calcNormFactors(as(t(data[,-1]), "matrix"));
            y <- voom(t(data[,-1]),plot=F,lib.size=colSums(t(data))*nf);
            data <- y$E;
            data <- cbind(fircol, t(data));
            colnames(data)[1] <- fircolnm;
            data <- as.data.frame(data);
            msg <- c(msg, paste("Performed log2-counts per million."));
        }else{
            print(paste("Unknown scaling parameter:", transform.opt));
        }
    }else{
        msg <- c(msg, paste("No data transformation was performed."));  
    }

    metadata <- dataSet$metadata;
    metadataname <- dataSet$firstcolname;;

    data <- cbind(metadata, data[,-1]);
    colnames(data)[1] <- metadataname;

    #using this object for plotting
    dataSet$norm2 <- data; #instead of dataSet$norm.phyobj <- phy.obj;
    current.msg <<- paste(msg, collapse=" ");
    dataSet$norm.msg2 <- current.msg;
    dataSet <<- dataSet;
    #save.image("norm2.RData")
    saveRDS(dataSet$norm2, file="norm.data2"); # save an copy
    return(1);
}




###################


SanityCheckData2 <- function(){
    
    sample_no <- nrow(dataSet$proc);
    sample_no2 <- nrow(dataSet$proc2);

    vari_no <- ncol(dataSet$proc)-1;
    vari_no2 <- ncol(dataSet$proc2)-1;

    

    sample_no_fil <- nrow(dataSet$filt.data);
    sample_no2_fil <- nrow(dataSet$filt.data2);

    vari_no_fil <- ncol(dataSet$filt.data)-1;
    vari_no2_fil <- ncol(dataSet$filt.data2)-1;


    sample_no_norm <- nrow(dataSet$norm1);
    sample_no2_norm <- nrow(dataSet$norm2);

    vari_no_norm <- ncol(dataSet$norm1)-1;
    vari_no2_norm <- ncol(dataSet$norm2)-1;


    return(c(0, sample_no, sample_no2, vari_no, vari_no2, 
sample_no_fil, sample_no2_fil, vari_no_fil, vari_no2_fil, 
sample_no_norm, sample_no2_norm, vari_no_norm, vari_no2_norm));
}








###################
#For simple plots - summaryview.xhtml


qc.boxplot <- function(dat, imgNm, dpi=72, format){
   
  dpi = as.numeric(dpi);
  library('ggplot2');
  library('lattice');

  imgNm = paste(imgNm, ".", format, sep="");
  
  #Sampling for the features -- too many features not shown on the plot
  subgene <- 10000;
  if (ncol(dat)>subgene) {
    set.seed(28051968);
    sg  <- sample(ncol(dat), subgene);
    Mss <- dat[,sg,drop=FALSE];
  } else {
    Mss <- dat;
  }
  
  #Sampling for the samples -- too many samples not shown on the plot
  subsmpl <- 100;
  if (nrow(Mss)>subsmpl) {
    set.seed(28051968);
    ss  <- sample(nrow(Mss), subsmpl)
    Mss <- Mss[ss,,drop=FALSE]
  } else {
    Mss <- Mss
  }
  
  sample_id <- rep(seq_len(nrow(Mss)), each = ncol(Mss)); 
  values  <- as.numeric(as.matrix(t(Mss)));
  
  df = cbind(values, sample_id)
  
  df = data.frame(df)
  df$sample_id = factor(df$sample_id)
  xlower = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[1])
  xupper = unname(quantile(df$values, probs = c(0.01, 0.99), na.rm=TRUE)[2])
  bp = ggplot(df, aes(sample_id, values)) +
    ylab("Counts") + xlab("Samples") + scale_x_discrete(labels=colnames(dataSet$data.norm)) + ylim(xlower, xupper) + stat_boxplot(geom = "errorbar", color="black")+ geom_boxplot(outlier.size=0.5, outlier.alpha=0.4)
  bp = bp + coord_flip();
  
  Cairo (file=imgNm, width=8, height=6, unit="in",dpi=dpi, type=format, bg="white");
  print(bp);
  dev.off();
}

#first omics
PlotDataBox <- function(boxplotName, dpi, format){
  qc.boxplot(dataSet$norm1[,-1], boxplotName, dpi, format);
}

#second omics 
PlotDataBox2 <- function(boxplotName, dpi, format){
  qc.boxplot(dataSet$norm2[,-1], boxplotName, dpi, format);
}








qc.pcaplot <- function(x, imgNm, dpi=72, format, Factor){
  dpi = as.numeric(dpi);
  imgNm = paste(imgNm, ".", format, sep="");
  library('lattice');
  library('ggplot2');
  library('ggrepel');
  
  pca <- prcomp(na.omit(x));
  names <- rownames(x);
  pca.res <- as.data.frame(pca$x);
  
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  
  if(all(Factor == 0)){
      pcafig <- ggplot(pca.res, aes(x=PC1, y=PC2, label=rownames(pca.res))) +
        geom_point(size=4) + xlim(xlim)+ ylim(ylim) + geom_text_repel(force=1.5)
  }else{
      pca.res <- cbind(rep(0, nrow(pca.res)), pca.res)
      colnames(pca.res)[1] <- "Factor"
      pca.res$Factor <- Factor
      pcafig <- ggplot(pca.res, aes(x=PC1, y=PC2,  color=factor(Factor), label=rownames(pca.res))) +
    geom_point(size=4) + xlim(xlim)+ ylim(ylim) + geom_text_repel(force=1.5)
  }
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
}


PlotDataPCA <- function(pcaName, dpi, format){
  qc.pcaplot(dataSet$norm1[,-1], pcaName, dpi, format, dataSet$norm1[,1]);
}

PlotDataPCA2 <- function(pcaName, dpi, format){
  qc.pcaplot(dataSet$norm2[,-1], pcaName, dpi, format, dataSet$norm2[,1]);
}







qc.density<- function(dat, imgNm, dpi=72, format, Factor){
  library("ggplot2")
  imgNm = paste(imgNm, ".", format, sep="");
  dpi = as.numeric(dpi);
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");  
  df = data.frame(t(dat), stringsAsFactors = FALSE);
  colnames(df) <- rownames(dat);  
  df2 = utils::stack(df); #Really useful function to melt...
  
  if(all(factor == 0)){
    g <- ggplot(df2, aes(x=values)) + geom_density()
  }else{
    conv = data.frame(ind=rownames(dat), class=factor(Factor))
    #conv$ind=gsub("[^[:alnum:]]", ".", conv$ind) #Remove everything but numbers and letters
    df3 = merge(df2, conv, by="ind")
    g <- ggplot(df3, aes(x=values, color=class)) + geom_density()
  }
  
  print(g);
  dev.off();
}


PlotDataDensity <- function(densityName, dpi, format){
  qc.density(dataSet$norm1[,-1], densityName, dpi, format, dataSet$norm1[,1]);
}

PlotDataDensity2 <- function(densityName, dpi, format){
  qc.density(dataSet$norm2[,-1], densityName, dpi, format, dataSet$norm2[,1]);
}







qc.meanstd <- function(dat, imgNm,dpi=72, format="png"){
  suppressMessages(require("vsn"));
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, ".", format, sep="");
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  plot <- suppressMessages(meanSdPlot(as.matrix(t(dat)), ranks=FALSE)); 
  print(plot);
  dev.off();
}



PlotDataMeanStd <- function(meanStdName, dpi,format){
  qc.meanstd(dataSet$norm1[,-1], meanStdName, dpi, format);
}

PlotDataMeanStd2 <- function(meanStdName, dpi,format){
  qc.meanstd(dataSet$norm2[,-1], meanStdName, dpi, format);
}






qc.norm <- function(dat, imgNm,dpi=72, format="png"){
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgNm, ".", format, sep="");
  
  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
  
  plot <- qqnorm(c(as.matrix(t(dat))), main = "QQnorm"); qqline(c(as.matrix(t(dat)))); 
  print(plot);
  dev.off();
}


PlotDataNorm <- function(normName, dpi,format){
  qc.norm(dataSet$norm1[,-1], normName, dpi, format);
}

PlotDataNorm2 <- function(normName, dpi,format){
  qc.norm(dataSet$norm2[,-1], normName, dpi, format);
}



####################

#O2PLS variable pick

#adj R^2

crossval_o2m_adjR2_fixed <- function(X, Y, a, ax, ay, nr_folds, nr_cores = 1,
                               stripped = TRUE, p_thresh = 3000, 
                               q_thresh = p_thresh, tol = 1e-10, max_iterations = 100)
{
  tic = proc.time()
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...")}
  kcv = nr_folds
  stopifnot(ncol(X) > max(a)+max(ax) , ncol(Y) > max(a)+max(ay) , nrow(X) >= kcv)
  stopifnot(nr_cores == abs(round(nr_cores)))
  if(nr_folds==1){stop("Cross-validation with 1 fold does not make sense, use 2 folds or more")}
  cl_crossval_o2m <- NULL
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  parms = data.frame(a = a)
  parms = apply(parms,1,as.list)
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(O2PLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      parms = data.frame(nx = ax)
      parms = merge(parms,data.frame(ny = ay))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, ax, ay,
                                                      stripped = stripped, p_thresh = p_thresh, 
                                                      q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))), 
                      nrow = length(ay), byrow=TRUE)
      nxny = which(R2grid == max(R2grid), arr.ind = TRUE)[1,]
      a_mse = suppressMessages(loocv_combi(X,Y,e$a,ax[nxny[2]],ay[nxny[1]],app_err=F,func=o2m,kcv=kcv,
                          stripped = stripped, p_thresh = p_thresh, 
                          q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
      c(a_mse, e$a, ax[nxny[2]],ay[nxny[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      parms = data.frame(nx = ax)
      parms = merge(parms,data.frame(ny = ay))
      parms = apply(parms,1,as.list)
      R2grid = matrix(colMeans(suppressMessages(adjR2(Y, X, e$a, ax, ay,
                                     stripped = stripped, p_thresh = p_thresh, 
                                     q_thresh = q_thresh, tol = tol, max_iterations = max_iterations))), 
                      nrow = length(ay), byrow=TRUE)
      nxny = which(R2grid == max(R2grid), arr.ind = TRUE)[1,]
      a_mse = suppressMessages(loocv_combi(X,Y,e$a,ax[nxny[2]],ay[nxny[1]],app_err=F,func=o2m,kcv=kcv,
                          stripped = stripped, p_thresh = p_thresh, 
                          q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
      c(a_mse, e$a, ax[nxny[2]],ay[nxny[1]])
    })
  }
  outp2 = matrix(unlist(outp), nrow = length(a), byrow = T)
  outp2 <- as.data.frame(outp2)
  outp2 <- cbind(outp2, rep(round((proc.time() - tic)[3],2), nrow(outp2)))
  names(outp2) <- c("MSE", "n", "nx", "ny", "time")
  message("minimum is at n = ", outp2[,2][which.min(outp2[,1])], sep = ' ')
  message("Elapsed time: ", round((proc.time() - tic)[3],2), " sec")
  return(outp2)
}


PlotDatatable <- function(tableName, dpi, format, jointc, orth1c, orth2c){
    library(OmicsPLS);
    library(gridExtra);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1")
    data2 <- readRDS("norm.data2")

    jointt <- as.integer(jointc);
    orth1t <- seq.int(0, as.integer(orth1c), 2);
    orth2t <- seq.int(0, as.integer(orth2c), 2);

    table_cv <- crossval_o2m_adjR2_fixed(data1[,-1], data2[,-1], 1:jointt, orth1t, orth2t, nr_folds = 4, nr_cores = 4)
    table_cv <- table_cv[order(table_cv$MSE), ]

    if(dim(table_cv)[1] >= 5){
      for(i in 1:5){ #output the best optimal 5 models based on MSE
          name <- paste0("fit", i)
          assign(name, o2m(data1[,-1], data2[,-1], table_cv$n[i], table_cv$nx[i], table_cv$ny[i]))
        }
    }else{
      for(i in 1:dim(table_cv)[1]){ #output the best optimal 5 models based on MSE
          name <- paste0("fit", i)
          assign(name, o2m(data1[,-1], data2[,-1], table_cv$n[i], table_cv$nx[i], table_cv$ny[i]))
        }
    }        


    if(dim(table_cv)[1] >= 5){
      summary_variance <- data.frame(matrix(ncol = 6, nrow = 5))
        x <- c("Joint#", "Orthogonal_omics1#", "Orthogonal_omics2", "omics1_joint_explained_by_omics2_joint", "omics2_joint_explained_by_omics1_joint", "MSE")
        colnames(summary_variance) <- x
    }else{
      summary_variance <- data.frame(matrix(ncol = 6, nrow = dim(table_cv)[1]))
        x <- c("Joint#", "Orthogonal_omics1#", "Orthogonal_omics2", "omics1_joint_explained_by_omics2_joint", "omics2_joint_explained_by_omics1_joint", "MSE")
        colnames(summary_variance) <- x
    }   

    if(dim(table_cv)[1] >= 5){
          for(i in 1:5){ #output the best optimal 5 models based on MSE
          name <- paste0("fit", i)
          assign(name, o2m(data1[,-1], data2[,-1], table_cv$n[i], table_cv$nx[i], table_cv$ny[i]))
          getting <- get(name)
          vec <- c(table_cv$n[i],
                   table_cv$nx[i],
                   table_cv$ny[i],
                   round(getting$R2Xhat/getting$R2Xcorr * 100, 3),
                   round(getting$R2Yhat/getting$R2Ycorr * 100, 3),
                   table_cv$MSE[i]
                   )
          summary_variance[i,] <- vec
        }     
    }else{
        for(i in 1:dim(table_cv)[1]){ #output the best optimal 5 models based on MSE
          name <- paste0("fit", i)
          assign(name, o2m(data1[,-1], data2[,-1], table_cv$n[i], table_cv$nx[i], table_cv$ny[i]))
          getting <- get(name)
          vec <- c(table_cv$n[i],
                   table_cv$nx[i],
                   table_cv$ny[i],
                   round(getting$R2Xhat/getting$R2Xcorr * 100, 3),
                   round(getting$R2Yhat/getting$R2Ycorr * 100, 3),
                   table_cv$MSE[i]
                   )
          summary_variance[i,] <- vec
        }       
    }    
    
    Cairo(file=imgNm, width=13, height=2, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(round(summary_variance, 5));
    dev.off();
    
    current.msg<<-"The number of variables are chosen with Cross validating adjusted R square.";
    dataSet$o2pls.varchoicemsg <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.adjR2 <- TRUE;
    analSet$varcriteria <- "adjR2";

    imgSet$o2pls.adjR2table<- summary_variance;
    imgSet<<-imgSet;
    analSet<<-analSet;
    
    return(1);
}



PlotDatatable2 <- function(tableName, dpi, format){
  library(OmicsPLS);
  dpi <- as.numeric(dpi)
  imgNm <- paste(tableName, ".", format, sep="");
  
  summary_variance <- imgSet$o2pls.adjR2table
  
  Cairo(file=imgNm, width=6, height=3, type=format, bg="white", dpi=dpi, unit="in");
  plot(summary_variance$MSE, main = "MSE", ylab = "MSE");
  dev.off();
}

PlotDatatable3 <- function(tableName, dpi, format){
  library(OmicsPLS);
  dpi <- as.numeric(dpi)
  imgNm <- paste(tableName, ".", format, sep="");
  
  summary_variance <- imgSet$o2pls.adjR2table
  
  Cairo(file=imgNm, width=6, height=3, type=format, bg="white", dpi=dpi, unit="in");
  plot(summary_variance$omics1_joint_explained_by_omics2_joint,
       main = "omics1_joint_explained_by_omics2_joint", ylab = "%");
  dev.off();
}
     
PlotDatatable4 <- function(tableName, dpi, format){
  library(OmicsPLS);
  dpi <- as.numeric(dpi)
  imgNm <- paste(tableName, ".", format, sep="");
  
  summary_variance <- imgSet$o2pls.adjR2table
  
  Cairo(file=imgNm, width=6, height=3, type=format, bg="white", dpi=dpi, unit="in");
  plot(summary_variance$omics2_joint_explained_by_omics1_joint,
       main = "omics2_joint_explained_by_omics1_joint", ylab = "%");
  dev.off();
}   






#prediction MSE

crossval_o2m_fixed <- function(X, Y, a, ax, ay, nr_folds, nr_cores = 1, 
                         stripped = TRUE, p_thresh = 3000, 
                         q_thresh = p_thresh, tol = 1e-10, max_iterations = 100) {
  tic = proc.time()
  if(any(abs(colMeans(X)) > 1e-5)){message("Data is not centered, proceeding...")}
  kcv = nr_folds
  stopifnot(ncol(X) > max(a)+max(ax) , ncol(Y) > max(a)+max(ay) , nrow(X) >= kcv)
  stopifnot(nr_cores == abs(round(nr_cores)))
  if(nr_folds==1){stop("Cross-validation with 1 fold does not make sense, use 2 folds or more")}
  
  parms = data.frame(nx = ax)
  parms = merge(parms,data.frame(ny = ay))
  parms = merge(parms,data.frame(a = a))
  parms = apply(parms,1,as.list)
  cl_crossval_o2m <- NULL
  
  on.exit({if(!is.null(cl_crossval_o2m)) stopCluster(cl_crossval_o2m)})
  
  if(Sys.info()[["sysname"]] == "Windows" && nr_cores > 1){
    cl_crossval_o2m <- makePSOCKcluster(nr_cores)
    clusterEvalQ(cl_crossval_o2m, library(O2PLS))
    clusterExport(cl_crossval_o2m, varlist = ls(), envir = environment())
    outp=parLapply(cl_crossval_o2m,parms,function(e){
      suppressMessages(loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                  stripped = stripped, p_thresh = p_thresh, 
                  q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
    })
  } else {
    outp=mclapply(mc.cores=nr_cores,parms,function(e){
      suppressMessages(loocv_combi(X,Y,e$a,e$nx,e$ny,app_err=F,func=o2m,kcv=kcv,
                  stripped = stripped, p_thresh = p_thresh, 
                  q_thresh = q_thresh, tol = tol, max_iterations = max_iterations)[[1]])
    })
  }
  dnams = list(paste("ax=",ax,sep=""),paste("ay=",ay,sep=""),paste("a=",a,sep=""))
  outp1 = array(unlist(outp),dim=c(length(ax),length(ay),length(a)),
                dimnames=dnams)
  outp = aperm(outp1,order(dim(outp1),decreasing=TRUE))
  dnams = dimnames(outp)
  
  if(dim(outp)[3]==1){
    dim(outp) = dim(outp)[-3]
    dimnames(outp) = dnams[1:2]
  }
  
  outpt = list(Original=outp1,Sorted=outp,kcv=kcv)
  #class(outpt) <- "cvo2m"
  toc = proc.time() - tic
  outpt$time = round(toc[3],2)
  x <- outpt
  wmCV = which(min(x$Or)==x$Or,TRUE,FALSE)
  dnams = dimnames(x$Or)
  dnams1 = dnams[[1]][wmCV[1]]
  dnams2 = dnams[[2]][wmCV[2]]
  dnams3 = dnams[[3]][wmCV[3]]
  
  
  summary_variance <- data.frame(matrix(ncol = 6, nrow = 1))
  zz <- c("Joint#", "Orthogonal_omics1#", "Orthogonal_omics2#", "omics1_joint_explained_by_omics2_joint", "omics2_joint_explained_by_omics1_joint", "MSE")
  colnames(summary_variance) <- zz
  
  getting <- o2m(X, Y, as.numeric(strsplit(dnams3, "=")[[1]][2]),
                   as.numeric(strsplit(dnams1, "=")[[1]][2]),
                   as.numeric(strsplit(dnams2, "=")[[1]][2]))

  vec <- c(strsplit(dnams3, "=")[[1]][2],
           strsplit(dnams1, "=")[[1]][2],
           strsplit(dnams2, "=")[[1]][2],
           round(getting$R2Xhat/getting$R2Xcorr * 100, 3),
           round(getting$R2Yhat/getting$R2Ycorr * 100, 3),
           round(min(x$Sor), 5)
           )
  summary_variance[1,] <- vec
  
  return(summary_variance)
}


PlotDatatableprediction <- function(tableName, dpi, format, jointc, orth1c, orth2c){
    library(OmicsPLS);
    library(gridExtra);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");

    jointt <- as.integer(jointc);
    orth1t <- seq.int(0, as.integer(orth1c), 2);
    orth2t <- seq.int(0, as.integer(orth2c), 2);

    
    table_cv <- crossval_o2m_fixed(data1[,-1], data2[,-1], 1:jointt, orth1t, orth2t, nr_folds = 4, nr_cores = 4);
    
    Cairo(file=imgNm, width=13, height=1, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(table_cv);
    dev.off();
    
    current.msg<<-"The number of variables are chosen with Cross validating prediction MSE.";
    dataSet$o2pls.varchoicemsg <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.predmse <- TRUE;

    analSet$varcriteria <- "predmse";

    imgSet$o2pls.predmsetable<- table_cv;
    imgSet<<-imgSet;
    analSet<<-analSet;
    
    return(1);
}



#pick one of joint and orth by users (option 4)

PlotDatatableuserpick <- function(tableName, dpi, format, joint, orth1, orth2){
    library(OmicsPLS);
    library(gridExtra);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    

    getting <- o2m(data1[,-1], data2[,-1], joint, orth1, orth2);
    summary_variance <- data.frame(matrix(ncol = 5, nrow = 1))
    x <- c("Joint#", "Orthogonal_omics1#", "Orthogonal_omics2", "omics1_joint_explained_by_omics2_joint", "omics2_joint_explained_by_omics1_joint")
    colnames(summary_variance) <- x
    
    vec <- c(joint,
               orth1,
               orth2,
               round(getting$R2Xhat/getting$R2Xcorr * 100, 3),
               round(getting$R2Yhat/getting$R2Ycorr * 100, 3)
               )
      summary_variance[1,] <- vec

    Cairo(file=imgNm, width=12, height=1, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(summary_variance);
    dev.off();
    
    current.msg<<-"The number of variables are picked by the user.";
    dataSet$o2pls.varchoicemsg <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.user <- TRUE;

    analSet$varcriteria <- "user";

    imgSet$o2pls.usertable<- summary_variance;
    imgSet<<-imgSet;
    analSet<<-analSet;
    
    return(1);
}





#(option3) fixed orth + joint

PlotDatatablefix1 <- function(tableName, dpi, format, joint, orth1, orth2, loops){
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    #joint <- as.integer(joint);
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    loops <- as.integer(loops);
    
    vec <- c()
    for (i in 1:loops){
      name <- "fit"
      assign(name, o2m(data1[,-1], data2[,-1], i, orth1, orth2))
      getting <- get(name)
      vec[i] <- ssq(getting$Tt %*% t(getting$C.))
    }

    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    pics <- plot(vec, main = "Scree plot: SSQ of prediction of Omics2 based on joint Omics1 scores", 
                 xlab = "# of joint component", ylab = "SSQ of Omics2 hat", col = "red", pch = "o", lty = 1);
    print(pics);
    dev.off();
    
    bases <- paste("The number of variables are picked by fixing joint and orthogonals based on joint: ",
          joint, "orth1: ", orth1, "orth2: ", orth2, ".", sep="");
    
    fixplot1 <- paste("Plot SSQ of prediction of Omics2 based on joint Omics1 scores - joint: ",
          loops, "orth1: ", orth1, "orth2: ", orth2, ".", sep="");
    
    current.msg<<-paste(c(bases,fixplot1), collapse=" ");
    dataSet$o2pls.varchoicemsg1 <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.fix <- TRUE;

    analSet$varcriteria <- "fix";
    
    analSet<<-analSet;
    
    return(1);
}


PlotDatatablefix2 <- function(tableName, dpi, format, joint, orth1, orth2, loops){
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    loops <- as.integer(loops);
    
    vec <- c()
    for (i in 1:loops){
      name <- "fit"
      assign(name, o2m(data1[,-1], data2[,-1], i, orth1, orth2))
      getting <- get(name)
      vec[i] <- ssq(getting$U %*% t(getting$W.))
    }

    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    pics <- plot(vec, main = "Scree plot: SSQ of prediction of Omics1 based on joint Omics2 scores", 
                 xlab = "# of joint component", ylab = "SSQ of Omics1 hat", col = "red", pch = "o", lty = 1);
    print(pics);
    dev.off();
    
    fixplot2 <- paste("Plot SSQ of prediction of Omics1 based on joint Omics2 scores - joint: ",
          loops, "orth1: ", orth1, "orth2: ", orth2, ".", sep="");
    
    current.msg<<-fixplot2;
    dataSet$o2pls.varchoicemsg2 <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.fix <- TRUE;
    analSet$varcriteria <- "fix";
    analSet<<-analSet;
    
    return(1);
}


PlotDatatablefix3 <- function(tableName, dpi, format, joint, orth1, orth2, loops){
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    orth1 <- as.integer(orth1);
    #orth2 <- as.integer(orth2);
    loops <- as.integer(loops);
    
    vec <- c()
    for (i in 1:loops){
      name <- "fit"
      assign(name, o2m(data1[,-1], data2[,-1], joint, orth1, i))
      getting <- get(name)
      vec[i] <- ssq(getting$U_Xosc %*% t(getting$P_Xosc.))
    }

    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    pics <- plot(vec, main = "Scree plot: SSQ for Orthogonal part in Omics2 per component", 
                 xlab = "# of orth2 component", ylab = "SSQ", col = "red", pch = "o", lty = 1);
    print(pics);
    dev.off();
    
    fixplot3 <- paste("Plot SSQ of orthogonal part in omics2 per ortho2 component changes - joint: ",
          joint, "orth1: ", orth1, "orth2: ", loops, ".", sep="");
    
    current.msg<<-fixplot3;
    dataSet$o2pls.varchoicemsg3 <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.fix <- TRUE;
    
    analSet$varcriteria <- "fix";

    analSet<<-analSet;
    
    return(1);
}


PlotDatatablefix4 <- function(tableName, dpi, format, joint, orth1, orth2, loops){
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    #orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    loops <- as.integer(loops);
    
    vec <- c()
    for (i in 1:loops){
      name <- "fit"
      assign(name, o2m(data1[,-1], data2[,-1], joint, i, orth2))
      getting <- get(name)
      vec[i] <- ssq(getting$T_Yosc %*% t(getting$P_Yosc.))
    }

    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    pics <- plot(vec, main = "Scree plot: SSQ for Orthogonal part in Omics1 per component", 
                 xlab = "# of orth1 component", ylab = "SSQ", col = "red", pch = "o", lty = 1);
    print(pics);
    dev.off();
    
    fixplot4 <- paste("Plot SSQ of orthogonal part in omics1 per ortho1 component changes - joint: ",
          joint, "orth1: ", loops, "orth2: ", orth2, ".", sep="");
    
    current.msg<<-fixplot4;
    dataSet$o2pls.varchoicemsg4 <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.fix <- TRUE;
    analSet$varcriteria <- "fix";
    analSet<<-analSet;
    
    return(1);
}


PlotDatatablefix4 <- function(tableName, dpi, format, joint, orth1, orth2, loops){
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    #orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    loops <- as.integer(loops);
    
    vec <- c()
    for (i in 1:loops){
      name <- "fit"
      assign(name, o2m(data1[,-1], data2[,-1], joint, i, orth2))
      getting <- get(name)
      vec[i] <- ssq(getting$T_Yosc %*% t(getting$P_Yosc.))
    }

    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    pics <- plot(vec, main = "Scree plot: SSQ for Orthogonal part in Omics1 per component", 
                 xlab = "# of orth1 component", ylab = "SSQ", col = "red", pch = "o", lty = 1);
    print(pics);
    dev.off();
    
    fixplot4 <- paste("Plot SSQ of orthogonal part in omics1 per ortho1 component changes - joint: ",
          joint, "orth1: ", loops, "orth2: ", orth2, ".", sep="");
    
    current.msg<<-fixplot4;
    dataSet$o2pls.varchoicemsg4 <- current.msg;
    dataSet <<- dataSet;
    analSet$o2pls <- TRUE;
    analSet$o2pls.fix <- TRUE;
    analSet$varcriteria <- "fix";
    analSet<<-analSet;
    
    return(1);
}



####################

#O2PLS main

#table 1
PlotMaintable1 <- function(tableName, dpi, format, joint, orth1, orth2){
    library(gridExtra);
    library(OmicsPLS);

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    joint <- as.integer(joint);
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    fircolnm <- colnames(data1)[1];

    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    
 

    fit <- o2m(data1[,-1], data2[,-1], joint, orth1, orth2);

    summary_score_noise <- data.frame(matrix(ncol = 6, nrow = 2))
    x <- c("omics1_joint_score", "omics1_ortho_score", 
           "omics2_joint_score", "omics2_ortho_score", 
           "omics1_noise", "omics2_noise"
           )
    colnames(summary_score_noise) <- x
    rownames(summary_score_noise) <- c("Absolute", "Relative")
    
    summary_score_noise[1,] <- c(sum(summary(fit)[12]$flags$varXjoint),
                                 sum(summary(fit)[12]$flags$varXorth),
                                 sum(summary(fit)[12]$flags$varYjoint),
                                 sum(summary(fit)[12]$flags$varYorth),
                                 (summary(fit)[12]$flags$ssqX - sum(summary(fit)[12]$flags$varXjoint) - sum(summary(fit)[12]$flags$varXorth)),
                                 (summary(fit)[12]$flags$ssqY - sum(summary(fit)[12]$flags$varYjoint) - sum(summary(fit)[12]$flags$varYorth))
                                 )
    
    summary_score_noise[2,] <- c(round(fit$R2Xcorr * 100, 3), 
                                 (round(fit$R2X * 100, 3) - round(fit$R2Xcorr * 100, 3)),
                                  round(fit$R2Ycorr * 100, 3),
                                  (round(fit$R2Y * 100, 3) - round(fit$R2Ycorr * 100, 3)),
                                 (100 - round(fit$R2X * 100, 3)),
                                 (100 - round(fit$R2Y * 100, 3))
                                 )

    
    Cairo(file=imgNm, width=13, height=2, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(round(summary_score_noise, 5));
    dev.off();
    
    current.msg<<-"Absolute and relative variation summary.";
    dataSet$o2pls.main1msg <- current.msg;
    dataSet$o2pls.fitmsg <- paste("Final Fit for O2PLS - joint: ",
          joint, "orth1: ", orth1, "orth2: ", orth2, ".", sep="");

    dataSet$fit <- fit;
    dataSet$fit.joint <- joint;
    dataSet$fit.orth1 <- orth1;
    dataSet$fit.orth2 <- orth2;
    dataSet <<- dataSet;
    analSet$o2plsrun <- TRUE;

    imgSet$o2pls.main1<- summary_score_noise;
    imgSet<<-imgSet;
    
    return(1);
}


#Main - table2 
PlotMaintable2 <- function(tableName, dpi, format, joint, orth1, orth2){

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    #data1 <- readRDS("norm.data1");
    #data2 <- readRDS("norm.data2");
    #joint <- as.integer(joint);
    #orth1 <- as.integer(orth1);
    #orth2 <- as.integer(orth2);
    
    fit <- dataSet$fit;

    model_summary <- data.frame(matrix(ncol = 3, nrow = summary(fit)[12]$flags$n + 1))
    x <- c("Component", "Omics1", "Omics2")
    colnames(model_summary) <- x
    y <- c(1:summary(fit)[12]$flags$n, "sum")
    model_summary[,1] <- y
    
    first <- c()
    second <- c()
    for (i in 1:summary(fit)[12]$flags$n){
      first[i] <- summary(fit)[12]$flags$varXjoint[i]/summary(fit)[12]$flags$ssqX
      second[i] <- summary(fit)[12]$flags$varYjoint[i]/summary(fit)[12]$flags$ssqY
    }
    
    first <- c(first, sum(first))
    second <- c(second, sum(second))
    
    model_summary[,2] <- round(first, 4) * 100
    model_summary[,3] <- round(second, 4) * 100
    
    Cairo(file=imgNm, width=5, height=3, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(model_summary);
    dev.off();
    
    current.msg<<-"Variation for each joint component for each omics.";
    dataSet$o2pls.main2msg <- current.msg;
    dataSet <<- dataSet;

    imgSet$o2pls.main2<- model_summary;
    imgSet<<-imgSet;
    
    return(1);
}



#Main - plot1
PlotMainplot1 <- function(tableName, dpi, format, orth1){
    library(reshape2);
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    #data1 <- readRDS("norm.data1");
    #data2 <- readRDS("norm.data2");
    #fircolnm <- colnames(data1)[1];
    orth1 <- as.integer(orth1);
    
    fit <- dataSet$fit;
    if(orth1 == 0){
      xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                 apply((fit$E)^2, 2, sum));
      xsq <- t(xsq);
      xsq_new <- xsq;

      colnames(xsq) <- c("joint", "noise");
      sum_xsq <- apply(xsq_new, 1, sum);
      xsq <- xsq[order(sum_xsq, decreasing = T)[1:20],]  #default plot before update (pick only top 20 entire SSQ)

      xsq <- melt(xsq);
      xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "noise"));      
      
      
    }else{
      xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                 apply((fit$T_Yosc %*% t(fit$P_Yosc.))^2, 2, sum), 
                 apply((fit$E)^2, 2, sum));
      xsq <- t(xsq);
      xsq_new <- xsq;
      
      colnames(xsq) <- c("joint", "orth", "noise");
      sum_xsq <- apply(xsq_new, 1, sum);
      xsq <- xsq[order(sum_xsq, decreasing = T)[1:20],] #default plot before update (pick only top 20 entire SSQ)

      xsq <- melt(xsq);
      xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "orth", "noise"));
      
    }

    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    g <- ggplot(xsq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
      labs(title = "Omics1: SSQ per variable", y = "SSQ", x = "Variable") +
       theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      guides(fill=guide_legend(title="Types"))
        
    print(g);
    dev.off();
    
    current.msg<<-paste("Here are the top 10 SSQ genes:", 
                        rownames(melt(sort(sum_xsq, decreasing = T)[1:10])),
                        sep="");
    dataSet$o2pls.main3msgall <- current.msg;
    
    dataSet$o2pls.main3msgjoint <- paste("Here are the top 10 joint SSQ genes:", 
                        rownames(melt(sort(xsq_new[,1], decreasing = T)[1:10])),
                        sep="");
    
    dataSet <<- dataSet;

    imgSet$o2pls.main3<- g;
    imgSet<<-imgSet;
    
    return(1);
}


#Main - plot2
PlotMainplot2 <- function(tableName, dpi, format, orth1){

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    #data1 <- readRDS("norm.data1");
    #data2 <- readRDS("norm.data2");
    #fircolnm <- colnames(data1)[1];
    
    orth1 <- as.integer(orth1);
    
    fit <- dataSet$fit;
    if(orth1 == 0){
      xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                 apply((fit$E)^2, 2, sum));
      xsq <- t(xsq);
      xsq_new <- xsq;

      colnames(xsq) <- c("joint", "noise");
      sum_xsq <- apply(xsq_new, 1, sum);
      xsq <- xsq[order(sum_xsq, decreasing = T)[1:20],]  #default plot before update (pick only top 20 entire SSQ)

      xsq <- melt(xsq);
      xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "noise"));      
      
      
    }else{
      xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                 apply((fit$T_Yosc %*% t(fit$P_Yosc.))^2, 2, sum), 
                 apply((fit$E)^2, 2, sum));
      
      xsq <- t(xsq);
      xsq_new <- xsq;

      colnames(xsq) <- c("joint", "orth", "noise");
      sum_xsq <- apply(xsq_new, 1, sum);
      xsq <- xsq[order(sum_xsq, decreasing = T)[1:20],]  #default plot before update (pick only top 20 entire SSQ)

      xsq <- melt(xsq);
      xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "orth", "noise"));
      
    }

    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    g <- ggplot(xsq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = "fill") +
  labs(title = "Omics1: SSQ per variable (normalized)", y = "SSQ", x = "Variable") +
   theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(fill=guide_legend(title="Types"))
        
    print(g);
    dev.off();
    
    imgSet$o2pls.main4<- g;
    imgSet<<-imgSet;
    
    return(1);
}




#Main - plot1 table update
PlotSSQTable <- function(imgName, dpi, format, omicsChoice, tableN, tableChoice, orth1, orth2){
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
  
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables
    tablechoice <- as.integer(tableChoice); #1: joint  2: orthogonal 3: entire SSQ
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    
    fit <- dataSet$fit;
    
    if(omicschoice == 1){

        if(tablen > dim(fit$W.)[1]){
            return(1);
        }

        if(orth1 == 0){
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          xsq <- t(xsq);
          xsq_new <- xsq;

          colnames(xsq) <- c("joint", "noise");
          sum_xsq <- apply(xsq_new, 1, sum);
          
        }else{
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$T_Yosc %*% t(fit$P_Yosc.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          
          xsq <- t(xsq);
          xsq_new <- xsq;

          colnames(xsq) <- c("joint", "orth", "noise");
          sum_xsq <- apply(xsq, 1, sum);
          
        }

        if(tablechoice == 1){
          tableans <- melt(sort(xsq_new[,1], decreasing = T)[1:tablen])
          colnames(tableans) <- "Joint"
        }else if(tablechoice == 2){
          tableans <- melt(sort(xsq_new[,2], decreasing = T)[1:tablen])
          colnames(tableans) <- "Orthogonal"
        }else{
          tableans <- melt(sort(sum_xsq, decreasing = T)[1:tablen])
          colnames(tableans) <- "Entire"
        }
      
    }else{

        if(tablen > dim(fit$C.)[1]){
            return(1);
        }
            
        if(orth2 == 0){
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum),
               apply((fit$Ff)^2, 2, sum))
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "noise");
        sum_ysq <- apply(ysq_new, 1, sum)
        
      }else{
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum), 
               apply((fit$U_Xosc %*% t(fit$P_Xosc.))^2, 2, sum), 
               apply((fit$Ff)^2, 2, sum))
        
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "orth", "noise")
        sum_ysq <- apply(ysq_new, 1, sum)

      }
        if(tablechoice == 1){
          tableans <- melt(sort(ysq_new[,1], decreasing = T)[1:tablen])
          colnames(tableans) <- "Joint"
        }else if(tablechoice == 2){
          tableans <- melt(sort(ysq_new[,2], decreasing = T)[1:tablen])
          colnames(tableans) <- "Orthogonal"
        }else{
          tableans <- melt(sort(sum_ysq, decreasing = T)[1:tablen])
          colnames(tableans) <- "Entire"
        }
      
    }
    
    tableans <- round(tableans, 5);
    
    Cairo(file=imgNm, width=8, height=nrow(tableans)/2, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();

    return(1);
}



#Main - plot1 update
PlotMainplotupdate1 <- function(tableName, dpi, format, orth1, orth2, omicsChoice, tablen, tablechoice){

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    randoms <- sample.int(28051968, 1);
    set.seed(randoms); 
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    tablen <- as.integer(tablen); #how many top n to choose
    tablechoice <- as.integer(tablechoice); #1. joint, 2. orth 3. entire

    dataSet$random.num <- randoms;
    dataSet <<- dataSet;
    fit <- dataSet$fit;
    
    if(omicschoice == 1){
        if(orth1 == 0){
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          xsq <- t(xsq);
          xsq_new <- xsq;

          colnames(xsq) <- c("joint", "noise");
          sum_xsq <- apply(xsq_new, 1, sum);
          
          if(tablechoice == 1){
            xsq <- xsq[order(xsq[,1], decreasing = T)[1:tablen],]
            
          }else if(tablechoice == 2){
            xsq <- xsq[order(xsq[,2], decreasing = T)[1:tablen],]

          }else{
            xsq <- xsq[order(sum_xsq, decreasing = T)[1:tablen],]
           
          }
          xsq <- melt(xsq);
          xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "noise"));  
          
        }else{
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$T_Yosc %*% t(fit$P_Yosc.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          
          xsq <- t(xsq);
          xsq_new <- xsq;

          colnames(xsq) <- c("joint", "orth", "noise");
          sum_xsq <- apply(xsq_new, 1, sum);
          
          if(tablechoice == 1){
            xsq <- xsq[order(xsq[,1], decreasing = T)[1:tablen],]
            
          }else if(tablechoice == 2){
            xsq <- xsq[order(xsq[,2], decreasing = T)[1:tablen],]

          }else{
            xsq <- xsq[order(sum_xsq, decreasing = T)[1:tablen],]
           
          }
          xsq <- melt(xsq);
          xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "orth", "noise"));

        }
        g <- ggplot(xsq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
        labs(title = "Omics1: SSQ per variable", y = "SSQ", x = "Variable") +
         theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        guides(fill=guide_legend(title="Types"))
    }else{
        if(orth2 == 0){
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum),
               apply((fit$Ff)^2, 2, sum))
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "noise");
        sum_ysq <- apply(ysq_new, 1, sum)

        if(tablechoice == 1){
           ysq <- ysq[order(ysq[,1], decreasing = T)[1:tablen],]

        }else if(tablechoice == 2){
          ysq <- ysq[order(ysq[,2], decreasing = T)[1:tablen],]

        }else{
          ysq <- ysq[order(sum_ysq, decreasing = T)[1:tablen],]

        } 
        ysq <- melt(ysq)
        ysq$Var2 <- factor(ysq$Var2, labels = c("joint", "noise"))
        
      }else{
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum), 
               apply((fit$U_Xosc %*% t(fit$P_Xosc.))^2, 2, sum), 
               apply((fit$Ff)^2, 2, sum))
        
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "orth", "noise")
        sum_ysq <- apply(ysq_new, 1, sum)

        if(tablechoice == 1){
            ysq <- ysq[order(ysq[,1], decreasing = T)[1:tablen],]
            
        }else if(tablechoice == 2){
          ysq <- ysq[order(ysq[,2], decreasing = T)[1:tablen],]

        }else{
          ysq <- ysq[order(sum_ysq, decreasing = T)[1:tablen],]

        }

        ysq <- melt(ysq)
        ysq$Var2 <- factor(ysq$Var2, labels = c("joint", "orth", "noise"))
        
      }
      g <- ggplot(ysq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
      labs(title = "Omics2: SSQ per variable", y = "SSQ", x = "Variable") +
       theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      guides(fill=guide_legend(title="Types"))
    }
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    print(g);
    dev.off();

    return(1);
}



#Main - plot2 update
PlotMainplotupdate2 <- function(tableName, dpi, format, orth1, orth2, omicsChoice, tablen, tablechoice){

    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  

    set.seed(dataSet$random.num);
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    orth1 <- as.integer(orth1);
    orth2 <- as.integer(orth2);
    tablen <- as.integer(tablen); #how many top n to choose
    tablechoice <- as.integer(tablechoice); #1. joint, 2. orth 3. entire
    
    fit <- dataSet$fit;
    
    if(omicschoice == 1){
        if(orth1 == 0){
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          xsq <- t(xsq);
          xsq_new <- xsq;
          
          colnames(xsq) <- c("joint", "noise");
          sum_xsq <- apply(xsq_new, 1, sum);
          
          if(tablechoice == 1){
            xsq <- xsq[order(xsq[,1], decreasing = T)[1:tablen],]
            
          }else if(tablechoice == 2){
            xsq <- xsq[order(xsq[,2], decreasing = T)[1:tablen],]

          }else{
            xsq <- xsq[order(sum_xsq, decreasing = T)[1:tablen],]
           
          }
          xsq <- melt(xsq);
          xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "noise")); 
          
          
        }else{
          xsq <- rbind(apply((fit$Tt %*% t(fit$W.))^2, 2, sum), 
                     apply((fit$T_Yosc %*% t(fit$P_Yosc.))^2, 2, sum), 
                     apply((fit$E)^2, 2, sum));
          
          xsq <- t(xsq);
          xsq_new <- xsq;
          
          colnames(xsq) <- c("joint", "orth", "noise");
          sum_xsq <- apply(xsq_new, 1, sum);
          
          if(tablechoice == 1){
            xsq <- xsq[order(xsq[,1], decreasing = T)[1:tablen],]
            
          }else if(tablechoice == 2){
            xsq <- xsq[order(xsq[,2], decreasing = T)[1:tablen],]

          }else{
            xsq <- xsq[order(sum_xsq, decreasing = T)[1:tablen],]
           
          }
          xsq <- melt(xsq);
          xsq$Var2 <- factor(xsq$Var2, labels = c("joint", "orth", "noise"));
        }
      
        g <- ggplot(xsq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = "fill") +
                labs(title = "Omics1: SSQ per variable (normalized)", y = "SSQ", x = "Variable") +
                 theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()) +
                guides(fill=guide_legend(title="Types"))
    }else{
        if(orth2 == 0){
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum),
               apply((fit$Ff)^2, 2, sum))
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "noise");
        sum_ysq <- apply(ysq_new, 1, sum)

        if(tablechoice == 1){
           ysq <- ysq[order(ysq[,1], decreasing = T)[1:tablen],]

        }else if(tablechoice == 2){
          ysq <- ysq[order(ysq[,2], decreasing = T)[1:tablen],]

        }else{
          ysq <- ysq[order(sum_ysq, decreasing = T)[1:tablen],]

        } 
        ysq <- melt(ysq)
        ysq$Var2 <- factor(ysq$Var2, labels = c("joint", "noise"))
        
        
      }else{
        ysq <- rbind(apply((fit$U %*% t(fit$C.))^2, 2, sum), 
               apply((fit$U_Xosc %*% t(fit$P_Xosc.))^2, 2, sum), 
               apply((fit$Ff)^2, 2, sum))
        
        ysq <- t(ysq)
        ysq_new <- ysq;

        colnames(ysq) <- c("joint", "orth", "noise")
        sum_ysq <- apply(ysq_new, 1, sum)

        if(tablechoice == 1){
            ysq <- ysq[order(ysq[,1], decreasing = T)[1:tablen],]
            
        }else if(tablechoice == 2){
          ysq <- ysq[order(ysq[,2], decreasing = T)[1:tablen],]

        }else{
          ysq <- ysq[order(sum_ysq, decreasing = T)[1:tablen],]

        }

        ysq <- melt(ysq)
        ysq$Var2 <- factor(ysq$Var2, labels = c("joint", "orth", "noise"))
        
      }
      g <- ggplot(ysq, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity", position = "fill") +
                labs(title = "Omics2: SSQ per variable (normalized)", y = "SSQ", x = "Variable") +
                 theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()) +
                guides(fill=guide_legend(title="Types"))
    }
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    print(g);
    dev.off();

    return(1);
}





#Main - plot3
PlotMainplot3 <- function(tableName, dpi, format){
    library(ggplot2);
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    options(scipen = 999);

    
    fit <- dataSet$fit;

    #works well even though there is no orthogonal -> removing factor orthogonal can give harder time for the users, so i will just leave it like this!
    load1x <- rbind(fit$W.[,1], fit$P_Yosc.[,1]);
    rownames(load1x) <- c("joint", "orth");
    load1x <- t(load1x);
    
    abssum_load1x <- apply(abs(load1x), 1, sum)
    load1x <- melt(load1x[order(abssum_load1x, decreasing = T)[1:20],]); #based on top 20 sum of absolute joint and orth

    load1x$Var2 <- factor(load1x$Var2, labels = c("joint", "orth"));
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    g <- ggplot(load1x, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
      labs(title = "Omics1: Loadings", y = "Loading value", x = "Variable") +
       theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      guides(fill=guide_legend(title="Types"))  
    print(g);
    dev.off();

    return(1);
}


#Main - plot3 table
PlotLoadingTable <- function(imgName, dpi, format, omicsChoice, tableN, tableChoice, orth1, orth2, joint, jointlimit, orth1limit, orth2limit){
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
  
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables
    tablechoice <- as.integer(tableChoice); #1. joint 2. orth 3. joint+orth
    orth1 <- as.integer(orth1); #which orth1 loading to show -> not working when omicschoice == 2
    orth2 <- as.integer(orth2);  #which orth2 loading to show -> not working when omicschoice == 1
    joint <- as.integer(joint); #which joint loading to show
    
    fit <- dataSet$fit;

    jointlimit <- as.integer(jointlimit);
    orth1limit <- as.integer(orth1limit);
    orth2limit <- as.integer(orth2limit);

    if( (jointlimit < joint) | (orth1limit < orth1) | (orth2limit < orth2)){
        return(1);
    }    
    
    if(omicschoice == 1){

        if(tablen > dim(fit$W.)[1]){
            return(1);
        }

        if(orth1 == 0){
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1+1]); #prevent any dimension error by adding 1, but its basically adding zero vector
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);
            
            abssum_load1x <- apply(abs(load1x), 1, sum)
            load1xall <- melt(sort(abssum_load1x, decreasing = T)[1:tablen]) #sum of absolute joint and orth 
            colnames(load1xall) <- "Absolute sum"
            load1xall <- cbind(load1xall, load1x[,1][order(abssum_load1x, decreasing = T)[1:tablen]])
            load1xall <- cbind(load1xall, load1x[,2][order(abssum_load1x, decreasing = T)[1:tablen]])
            colnames(load1xall)[c(2,3)] <- c("Original joint", "Original orth")
            
            load1xorth<- melt(sort(abs(load1x[,2]), decreasing = T)[1:tablen]) #sum of absolute orth
            colnames(load1xorth) <- "Absolute value" 
            load1xorth <- cbind(load1xorth, load1x[,2][order(abs(load1x[,2]), decreasing = T)[1:tablen]])
            colnames(load1xorth)[2] <- "Original orth"
            
            load1xjoint <- melt(sort(abs(load1x[,1]), decreasing = T)[1:tablen]) #sum of joint
            colnames(load1xjoint) <- "Absolute value" 
            load1xjoint <- cbind(load1xjoint, load1x[,1][order(abs(load1x[,1]), decreasing = T)[1:tablen]])
            colnames(load1xjoint)[2] <- "Original joint"
            
            if(tablechoice == 1){
              tableans <- load1xjoint;
            }else if(tablechoice == 2){
              tableans <- load1xorth;
            }else{
              tableans <- load1xall;
            }          
        }else{
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1]);
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);
            
            abssum_load1x <- apply(abs(load1x), 1, sum)
            load1xall <- melt(sort(abssum_load1x, decreasing = T)[1:tablen]) #sum of absolute joint and orth 
            colnames(load1xall) <- "Absolute sum"
            load1xall <- cbind(load1xall, load1x[,1][order(abssum_load1x, decreasing = T)[1:tablen]])
            load1xall <- cbind(load1xall, load1x[,2][order(abssum_load1x, decreasing = T)[1:tablen]])
            colnames(load1xall)[c(2,3)] <- c("Original joint", "Original orth")
            
            load1xorth<- melt(sort(abs(load1x[,2]), decreasing = T)[1:tablen]) #sum of absolute orth
            colnames(load1xorth) <- "Absolute value"
            load1xorth <- cbind(load1xorth, load1x[,2][order(abs(load1x[,2]), decreasing = T)[1:tablen]])
            colnames(load1xorth)[2] <- "Original orth"
            
            load1xjoint <- melt(sort(abs(load1x[,1]), decreasing = T)[1:tablen]) #sum of joint
            colnames(load1xjoint) <- "Absolute value" 
            load1xjoint <- cbind(load1xjoint, load1x[,1][order(abs(load1x[,1]), decreasing = T)[1:tablen]])
            colnames(load1xjoint)[2] <- "Original joint"
            
            if(tablechoice == 1){
              tableans <- load1xjoint;
            }else if(tablechoice == 2){
              tableans <- load1xorth;
            }else{
              tableans <- load1xall;
            }
        }

    }else{

        if(tablen > dim(fit$C.)[1]){
            return(1);
        }

        if(orth2 == 0){
            load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2+1])
            rownames(load1y) <- c("joint", "orth")
            load1y <- t(load1y)

            abssum_load1y <- apply(abs(load1y), 1, sum)
            load1yall <-melt(sort(abssum_load1y, decreasing = T)[1:tablen])
            colnames(load1yall) <- "Absolute sum"
            load1yall <- cbind(load1yall, load1y[,1][order(abssum_load1y, decreasing = T)[1:tablen]])
            load1yall <- cbind(load1yall, load1y[,2][order(abssum_load1y, decreasing = T)[1:tablen]])
            colnames(load1yall)[c(2,3)] <- c("Original joint", "Original orth")
            
            load1yorth<- melt(sort(abs(load1y[,2]), decreasing = T)[1:tablen]) 
            colnames(load1yorth) <- "Absolute value" 
            load1yorth <- cbind(load1yorth, load1y[,2][order(abs(load1y[,2]), decreasing = T)[1:tablen]])
            colnames(load1yorth)[2] <- "Original orth"
            
            load1yjoint <-melt(sort(abs(load1y[,1]), decreasing = T)[1:tablen])
            colnames(load1yjoint) <- "Absolute value"
            load1yjoint <- cbind(load1yjoint, load1y[,1][order(abs(load1y[,1]), decreasing = T)[1:tablen]])
            colnames(load1yjoint)[2] <- "Original joint"

            if(tablechoice == 1){
              tableans <- load1yjoint;
            }else if(tablechoice == 2){
              tableans <- load1yorth;
            }else{
              tableans <- load1yall;
            }        
        }else{
            load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2])
            rownames(load1y) <- c("joint", "orth")
            load1y <- t(load1y)

            abssum_load1y <- apply(abs(load1y), 1, sum)
            load1yall <-melt(sort(abssum_load1y, decreasing = T)[1:tablen])
            colnames(load1yall) <- "Absolute sum"
            load1yall <- cbind(load1yall, load1y[,1][order(abssum_load1y, decreasing = T)[1:tablen]])
            load1yall <- cbind(load1yall, load1y[,2][order(abssum_load1y, decreasing = T)[1:tablen]])
            colnames(load1yall)[c(2,3)] <- c("Original joint", "Original orth")
            
            load1yorth<- melt(sort(abs(load1y[,2]), decreasing = T)[1:tablen]) 
            colnames(load1yorth) <- "Absolute value" 
            load1yorth <- cbind(load1yorth, load1y[,2][order(abs(load1y[,2]), decreasing = T)[1:tablen]])
            colnames(load1yorth)[2] <- "Original orth"
            
            load1yjoint <-melt(sort(abs(load1y[,1]), decreasing = T)[1:tablen])
            colnames(load1yjoint) <- "Absolute value"
            load1yjoint <- cbind(load1yjoint, load1y[,1][order(abs(load1y[,1]), decreasing = T)[1:tablen]])
            colnames(load1yjoint)[2] <- "Original joint"

            if(tablechoice == 1){
              tableans <- load1yjoint;
            }else if(tablechoice == 2){
              tableans <- load1yorth;
            }else{
              tableans <- load1yall;
            }
      }
    }
    
    tableans <- round(tableans, 5);

    Cairo(file=imgNm, width=9, height=nrow(tableans)/2.7, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();
    return(1);
}
 




#Main - plot3 update
PlotMainplotupdate3 <- function(tableName, dpi, format, joint, orth1, orth2, omicsChoice, jointlimit, orth1limit, orth2limit, tablen, tablechoice){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    options(scipen = 999);
    
    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    orth1 <- as.integer(orth1); #which orth1 loading to show -> not working when omicschoice == 2
    orth2 <- as.integer(orth2);  #which orth2 loading to show -> not working when omicschoice == 1
    joint <- as.integer(joint); #which joint loading to show

    jointlimit <- as.integer(jointlimit);
    orth1limit <- as.integer(orth1limit);
    orth2limit <- as.integer(orth2limit);

    tablen <- as.integer(tablen); #top n variables
    tablechoice <- as.integer(tablechoice); #1. joint 2. orth 3. joint+orth


    if( (jointlimit < joint) | (orth1limit < orth1) | (orth2limit < orth2)){
        return(1);
    }

    fit <- dataSet$fit;
    
    
    if(omicschoice == 1){
        if(orth1 == 0){
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1+1]);
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);

            if(tablechoice == 1){
              load1x <- melt(load1x[order(abs(load1x[,1]), decreasing = T)[1:tablen],]);
            }else if(tablechoice == 2){
              load1x <- melt(load1x[order(abs(load1x[,2]), decreasing = T)[1:tablen],]);
            }else{
              abssum_load1x <- apply(abs(load1x), 1, sum)
              load1x <- melt(load1x[order(abssum_load1x, decreasing = T)[1:tablen],]);
            }  
            
            load1x$Var2 <- factor(load1x$Var2, labels = c("joint", "orth"));
        }else{
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1]);
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);

            if(tablechoice == 1){
              load1x <- melt(load1x[order(abs(load1x[,1]), decreasing = T)[1:tablen],]);
            }else if(tablechoice == 2){
              load1x <- melt(load1x[order(abs(load1x[,2]), decreasing = T)[1:tablen],]);
            }else{
              abssum_load1x <- apply(abs(load1x), 1, sum)
              load1x <- melt(load1x[order(abssum_load1x, decreasing = T)[1:tablen],]);
            }  
            
            load1x$Var2 <- factor(load1x$Var2, labels = c("joint", "orth"));
        }
        g <- ggplot(load1x, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
            labs(title = "Omics1: Loadings", y = "Loading value", x = "Variable") +
             theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            guides(fill=guide_legend(title="Types"))  
    }else{
      if(orth2 == 0){
          load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2+1])
          rownames(load1y) <- c("joint", "orth")
          load1y <- t(load1y)
 
          if(tablechoice == 1){
              load1y <- melt(load1y[order(abs(load1y[,1]), decreasing = T)[1:tablen],]);
            }else if(tablechoice == 2){
              load1y <- melt(load1y[order(abs(load1y[,2]), decreasing = T)[1:tablen],]);
            }else{
              abssum_load1y <- apply(abs(load1y), 1, sum)
              load1y <- melt(load1y[order(abssum_load1y, decreasing = T)[1:tablen],]);
            }  
            
          load1y$Var2 <- factor(load1y$Var2, labels = c("joint", "orth"));
          
      }else{
          load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2])
          rownames(load1y) <- c("joint", "orth")
          load1y <- t(load1y)

          if(tablechoice == 1){
              load1y <- melt(load1y[order(abs(load1y[,1]), decreasing = T)[1:tablen],]);
            }else if(tablechoice == 2){
              load1y <- melt(load1y[order(abs(load1y[,2]), decreasing = T)[1:tablen],]);
            }else{
              abssum_load1y <- apply(abs(load1y), 1, sum)
              load1y <- melt(load1y[order(abssum_load1y, decreasing = T)[1:tablen],]);
            }  
            
          load1y$Var2 <- factor(load1y$Var2, labels = c("joint", "orth"));
      }
      g <- ggplot(load1y, aes(x = Var1, y = value, fill = Var2)) + geom_bar(stat = "identity") +
            labs(title = "Omics2: Loadings", y = "Loading value", x = "Variable") +
             theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            guides(fill=guide_legend(title="Types"))  
    }
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    print(g);
    dev.off();

    return(1);
}




#Loading sig features correlations
PlotUnivCor <- function(imgName, dpi, format, omicsChoice, tableN, tableChoice, orth1, orth2, joint, jointlimit, orth1limit, orth2limit, correlationoption){
    library(pheatmap);
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
  
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables
    tablechoice <- as.integer(tableChoice); #1. joint 2. orth 3. joint+orth
    orth1 <- as.integer(orth1); #which orth1 loading to show -> not working when omicschoice == 2
    orth2 <- as.integer(orth2);  #which orth2 loading to show -> not working when omicschoice == 1
    joint <- as.integer(joint); #which joint loading to show
    
    fit <- dataSet$fit;

    jointlimit <- as.integer(jointlimit);
    orth1limit <- as.integer(orth1limit);
    orth2limit <- as.integer(orth2limit);
    
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    correlationoption <- as.integer(correlationoption); #1. pearson 2. spearman
    
    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
      if (length(delete.columns) > 0) {
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;

    if( (jointlimit < joint) | (orth1limit < orth1) | (orth2limit < orth2)){
        return(1);
    }    
    
    if(omicschoice == 1){

        if(tablen > dim(fit$W.)[1]){
            return(1);
        }

        if(orth1 == 0){
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1+1]); #prevent any dimension error by adding 1, but its basically adding zero vector
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);
            
            abssum_load1x <- apply(abs(load1x), 1, sum)
            
            if(tablechoice == 1){
              #add 1 since the first column is sample meta-data
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,1]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,1]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }else if(tablechoice == 2){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,2]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,2]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }  
            }else{
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abssum_load1x, decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abssum_load1x, decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }          
        }else{
            load1x <- rbind(fit$W.[,joint], fit$P_Yosc.[,orth1]);
            rownames(load1x) <- c("joint", "orth");
            load1x <- t(load1x);
            
            abssum_load1x <- apply(abs(load1x), 1, sum)
            
            if(tablechoice == 1){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,1]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,1]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }else if(tablechoice == 2){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,2]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abs(load1x[,2]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }  
            }else{
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abssum_load1x, decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data1[,order(abssum_load1x, decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }
        }
    }else{
        if(tablen > dim(fit$C.)[1]){
            return(1);
        }
        if(orth2 == 0){
            load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2+1])
            rownames(load1y) <- c("joint", "orth")
            load1y <- t(load1y)

            abssum_load1y <- apply(abs(load1y), 1, sum)

            if(tablechoice == 1){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,1]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,1]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }else if(tablechoice == 2){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,2]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,2]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }  
            }else{
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abssum_load1y, decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abssum_load1y, decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }        
        }else{
            load1y <- rbind(fit$C.[,joint], fit$P_Xosc.[,orth2])
            rownames(load1y) <- c("joint", "orth")
            load1y <- t(load1y)

            abssum_load1y <- apply(abs(load1y), 1, sum)
            
            if(tablechoice == 1){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,1]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,1]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }else if(tablechoice == 2){
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,2]), decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abs(load1y[,2]), decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }  
            }else{
              if(correlationoption == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abssum_load1y, decreasing = T)[1:tablen]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                pheatmap(cor(data2[,order(abssum_load1y, decreasing = T)[1:tablen]+1], method = "spearman"), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
                dev.off();
              }
            }
      }
    }
    return(1);
}



#Main - plot4 
PlotMainplot4 <- function(tableName, dpi, format, joint){
    library(dplyr);
    library(magrittr);
    library(stringr);
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    joint <- as.integer(joint);
    
    fit <- dataSet$fit;
    if(joint == 1){
        return(1);
    } 

    xjoints <- sqrt(loadings(fit, "Xjoint", 1:2) %>% raise_to_power(2) %>% rowSums) 
    xjoints[-(order(xjoints,decreasing=T)[1:10])] = 0
    xjoints <- sign(xjoints)
    

    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    g <- plot(fit, loading_name="Xj", i=1, j=2, label = "c", use_ggplot2 = TRUE, 
         alpha = xjoints,
         aes(label = c(stringr::str_sub(names(xjoints), start = 0))),
         hjust = rep(c(0, 1), length.out = length(xjoints)), size = 3) + 
      theme_bw() +
      geom_point(alpha = 0.2+0.7*xjoints, col = 'grey') + 
      labs(title = "Omics1: joint loadings",
           x = "Joint Loading 1", y = "Joint Loading 2") + 
      theme(plot.title = element_text(face='bold')) 
    print(g);
    dev.off();
    
    return(1);
}


#Main - plot4 table
PlotCompLoadingTable <- function(imgName, dpi, format, omicsChoice, tableN, loadingx, loadingy, joint){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables
    loadingx <- as.integer(loadingx); #choose variable for the x-axis loading 
    loadingy <- as.integer(loadingy); #choose variable for the y-axis loading 
    
    fit <- dataSet$fit;

    joint <- as.integer(joint);
    if( (joint == 1) | (joint < loadingx) | (joint < loadingy)){
        return(1);
    } 

    
    if(omicschoice == 1){

        if(tablen > dim(fit$W.)[1]){
            return(1);
        }

        xjoints <- sqrt(loadings(fit, "Xjoint", c(loadingx, loadingy)) %>% raise_to_power(2) %>% rowSums) 

        tablex <- melt(sort(xjoints,decreasing=T)[1:tablen])
        colnames(tablex) <- "Euclidean_dist"
        tablex <- cbind(tablex, loadings(fit, "Xjoint", c(loadingx, loadingy))[order(xjoints,decreasing=T)[1:tablen],])
        colnames(tablex)[c(2,3)] <- c("Loading_x-axis", "Loading_y-axis") 
        
        tableans <- tablex;
    }else{

        if(tablen > dim(fit$C.)[1]){
            return(1);
        }

        yjoints <- sqrt(loadings(fit, "Yjoint", c(loadingx, loadingy)) %>% raise_to_power(2) %>% rowSums)

        tabley <- melt(sort(yjoints,decreasing=T)[1:tablen])
        colnames(tabley) <- "Euclidean_dist"
        tabley <- cbind(tabley, loadings(fit, "Yjoint", c(loadingx, loadingy))[order(yjoints,decreasing=T)[1:tablen],])
        colnames(tabley)[c(2,3)] <- c("Loading_x-axis", "Loading_y-axis") 
        
        tableans <- tabley;
    }

    tableans <- round(tableans, 5);
    

    Cairo(file=imgNm, width=9.5, height=nrow(tableans)/2.8,
          type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();

    return(1);
}



#Main - plot4 update
PlotMainplotupdate4 <- function(tableName, dpi, format, loadingx, loadingy, omicsChoice, joint, tableN){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables
    loadingx <- as.integer(loadingx); #choose variable for the x-axis loading 
    loadingy <- as.integer(loadingy); #choose variable for the y-axis loading 
    
    joint <- as.integer(joint);
    if( (joint == 1) | (joint < loadingx) | (joint < loadingy)){
        return(1);
    }  
    
    fit <- dataSet$fit;

    
    if(omicschoice == 1){

        if(tablen > dim(fit$W.)[1]){
            return(1);
        }

        xjoints <- sqrt(loadings(fit, "Xjoint", c(loadingx, loadingy)) %>% raise_to_power(2) %>% rowSums) 
        
        xjoints[-(order(xjoints,decreasing=T)[1:tablen])] = 0
        xjoints <- sign(xjoints)
        
        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        g <- plot(fit, loading_name="Xj", i=loadingx, j=loadingy, label = "c", use_ggplot2 = TRUE, 
               alpha = xjoints,
               aes(label = c(stringr::str_sub(names(xjoints), start = 0))),
               hjust = rep(c(0, 1), length.out = length(xjoints)), size = 3) + 
            theme_bw() +
            geom_point(alpha = 0.2+0.7*xjoints, col = 'grey') + 
            labs(title = "Omics1: joint loadings",
                 x = paste0("Joint Loading ", loadingx), y = paste0("Joint Loading", loadingy)) + 
            theme(plot.title = element_text(face='bold')) 
        print(g);
        dev.off();
        
    }else{      
        if(tablen > dim(fit$C.)[1]){
            return(1);
        }

        yjoints <- sqrt(loadings(fit, "Yjoint", c(loadingx, loadingy)) %>% raise_to_power(2) %>% rowSums) 

        yjoints[-(order(yjoints,decreasing=T)[1:tablen])] = 0
        yjoints <- sign(yjoints)

        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        g <- plot(fit, loading_name="Yj", i=loadingx, j=loadingy, label = "c", use_ggplot2 = TRUE, 
               alpha = yjoints,
               aes(label = c(stringr::str_sub(names(yjoints), start = 0))),
               hjust = rep(c(0, 1), length.out = length(yjoints)), size = 3) + 
            theme_bw() +
            geom_point(alpha = 0.2+0.7*yjoints, col = 'grey') + 
            labs(title = "Omics2: joint loadings",
                 x = paste0("Joint Loading ", loadingx), y = paste0("Joint Loading", loadingy)) + 
            theme(plot.title = element_text(face='bold')) 
        print(g);
        dev.off();
    }
    
    return(1);
}






#Main - plot5 
PlotMainplot5 <- function(tableName, dpi, format){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    #options(scipen = 999);
    
    fit <- dataSet$fit;
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    g <- plot(fit$Tt[,1], fit$U[,1], main = "Least square: Joint1",
              xlab = "Joint Omics1",
              ylab = "Joint Omics2");abline(lm(fit$U[,1] ~ fit$Tt[,1]), col = "red")
    print(g);
    dev.off();
    
    return(1);
}




#Main - plot6 (table) 
PlotMainplot6 <- function(tableName, dpi, format){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    
    fit <- dataSet$fit;
    
    getting <- lm(fit$U[,1] ~ fit$Tt[,1])
    tableans <- summary(getting)[4]$coefficients
    rownames(tableans)[2] <- "Omics1 Jointscore1"

    tableans <- round(tableans, 5);
    
    Cairo(file=imgNm, width=8, height=1,
          type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();
    
    return(1);
}


#Main - plot6 (statistics string)
PerformQualityStat<-function(){
    suppressMessages(require(OmicsPLS));

    fit <- dataSet$fit;
    getting <- lm(fit$U[,1] ~ fit$Tt[,1]);
    stat.info <- paste("Adjusted R^2: ", signif(summary(getting)[9]$adj.r.squared, 5), " on sample: ", dim(fit$Tt)[1], sep=""); 

    analSet$quality.stat.info <- stat.info;
    analSet<<-analSet;

    return(stat.info);
}




#Main - plot 5 and 6 table
PlotPredictionTable <- function(imgName, dpi, format, omicsChoice, tableN, jointPred, jointb){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    
    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    jointpred <- as.integer(jointPred); #what joint component to compare
    joint <- as.integer(jointb); #joint component limit
    
    if(joint < jointpred){ #if user wants to compare joint outside of the joint dimensions
        return(1);
    } 
    
    fit <- dataSet$fit;

    if(tablen > nrow(fit$U)){ #if user wants to bring out more samples than we have
      return(1);
    }
    
    if(omicschoice == 1){
        getting <- lm(fit$U[,jointpred] ~ fit$Tt[,jointpred])
        
        tableans <- melt(sort(cooks.distance(getting), decreasing = F)[1:tablen])
        colnames(tableans) <- "Cook's distance"
        
        tableans <- cbind(tableans, fit$Tt[,jointpred][order(cooks.distance(getting), decreasing=F)[1:tablen]])
        colnames(tableans)[2] <- "X-axis(Omics1) score"
        
        tableans <- cbind(tableans, fit$U[,jointpred][order(cooks.distance(getting), decreasing=F)[1:tablen]])
        colnames(tableans)[3] <- "Y-axis(Omics2) score"
        
    }else{
        getting <- lm(fit$Tt[,jointpred] ~ fit$U[,jointpred])
        
        tableans <- melt(sort(cooks.distance(getting), decreasing = F)[1:tablen])
        colnames(tableans) <- "Cook's distance"
        
        tableans <- cbind(tableans, fit$U[,jointpred][order(cooks.distance(getting), decreasing=F)[1:tablen]])
        colnames(tableans)[2] <- "X-axis(Omics2) score"  
        
        tableans <- cbind(tableans, fit$Tt[,jointpred][order(cooks.distance(getting), decreasing=F)[1:tablen]])
        colnames(tableans)[3] <- "Y-axis(Omics1) score"
    
    }

    tableans <- round(tableans, 5);
    
    Cairo(file=imgNm, width=9.5, height=nrow(tableans)/2.2,
          type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();
    
    return(1);
}


#Main - plot 5 update
PlotMainplotupdate5 <- function(imgName, dpi, format, omicsChoice, tableN, jointPred, jointb, samplemeta){

    data1 <- readRDS("norm.data1");
    
    metadata <- data1[,1];


    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    data2 <- readRDS("norm.data2");
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    

    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    
    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    jointpred <- as.integer(jointPred); #what joint component to compare
    joint <- as.integer(jointb); #joint component limit
    metaoption <- as.integer(samplemeta); #1: none 2: apply 
    
    if(joint < jointpred){ #if user wants to compare joint outside of the joint dimensions
        return(1);
    } 
    
    fit <- dataSet$fit;
    if(tablen > nrow(fit$U)){ #if user wants to bring out more samples than we have
      return(1);
    }
    
    
    if(omicschoice == 1){
        getting <- lm(fit$U[,jointpred] ~ fit$Tt[,jointpred])
        
        datas <- data.frame(Omics1 = fit$Tt[,jointpred], Omics2 = fit$U[,jointpred])
        datas_copy <- datas
        datas_copy[-(order(cooks.distance(getting), decreasing=F)[1:tablen]), ] = 0
        datas_copy <- sign(datas_copy)
        datas_copy <- abs(datas_copy)
        
        if(metaoption == 1){
            g <- ggplot(datas, aes(x = Omics1, y = Omics2)) + 
              geom_point(alpha = 0.2+0.7*datas_copy[,1], col = 'grey')
            g <- g + geom_text(aes(label = c(stringr::str_sub(rownames(datas_copy), start = 0))),
                               hjust = rep(c(0, 1), length.out = nrow(datas_copy)), size = 3,
                               alpha = datas_copy[,1])
            g <- g + theme_bw() + labs(title = paste0("Least square: Joint", jointpred),
                     x = "Joint Omics1", y = "Joint Omics2") +
              theme(plot.title = element_text(face='bold')) 
            g <- g + stat_smooth(method="lm", se=FALSE, col = "red", size = 0.3)
        }else{
            datas <- cbind(metadata, datas)
            colnames(datas)[1] <- "sample_meta"

            g <- ggplot(datas, aes(x = Omics1, y = Omics2, col = sample_meta)) +
              geom_point(alpha = 0.2+0.7*datas_copy[,1])
            g <- g + geom_text(aes(label = c(stringr::str_sub(rownames(datas_copy), start = 0))),
                               hjust = rep(c(0, 1), length.out = nrow(datas_copy)), size = 3,
                               alpha = datas_copy[,1])
            g <- g + theme_bw() + labs(title = paste0("Least square: Joint", jointpred),
                     x = "Joint Omics1", y = "Joint Omics2", fill = "sample-meta") +
              theme(plot.title = element_text(face='bold')) 
            g <- g + stat_smooth(method="lm", se=FALSE, col = "red", size = 0.3)
        }
    }else{
        getting <- lm(fit$Tt[,jointpred] ~ fit$U[,jointpred])
        
        datas <- data.frame(Omics2 = fit$U[,jointpred], Omics1 = fit$Tt[,jointpred])
        datas_copy <- datas
        datas_copy[-(order(cooks.distance(getting), decreasing=F)[1:tablen]), ] = 0
        datas_copy <- sign(datas_copy)
        datas_copy <- abs(datas_copy)
        
        if(metaoption == 1){
            g <- ggplot(datas, aes(x = Omics2, y = Omics1)) +
              geom_point(alpha = 0.2+0.7*datas_copy[,1], col = 'grey')
            g <- g + geom_text(aes(label = c(stringr::str_sub(rownames(datas_copy), start = 0))),
                               hjust = rep(c(0, 1), length.out = nrow(datas_copy)), size = 3,
                               alpha = datas_copy[,1])
            g <- g + theme_bw() + labs(title = paste0("Least square: Joint", jointpred),
                     x = "Joint Omics2", y = "Joint Omics1") +
              theme(plot.title = element_text(face='bold')) 
            g <- g + stat_smooth(method="lm", se=FALSE, col = "red", size = 0.3)    
        }else{
            datas <- cbind(metadata, datas)
            colnames(datas)[1] <- "sample_meta"
            
            g <- ggplot(datas, aes(x = Omics2, y = Omics1, col = sample_meta)) +
              geom_point(alpha = 0.2+0.7*datas_copy[,1])
            g <- g + geom_text(aes(label = c(stringr::str_sub(rownames(datas_copy), start = 0))),
                               hjust = rep(c(0, 1), length.out = nrow(datas_copy)), size = 3,
                               alpha = datas_copy[,1])
            g <- g + theme_bw() + labs(title = paste0("Least square: Joint", jointpred),
                     x = "Joint Omics2", y = "Joint Omics1") +
              theme(plot.title = element_text(face='bold')) 
            g <- g + stat_smooth(method="lm", se=FALSE, col = "red", size = 0.3)   
        }  

    }
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    print(g);
    dev.off();
    
    return(1);
}


#Main - plot 6 (table) update
PlotMainplotupdate6 <- function(imgName, dpi, format, omicsChoice, tableN, jointPred, jointb){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    
    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    jointpred <- as.integer(jointPred); #what joint component to compare
    joint <- as.integer(jointb); #joint component limit
    
    if(joint < jointpred){ #if user wants to compare joint outside of the joint dimensions
        return(1);
    } 
    
    fit <- dataSet$fit;
    if(tablen > nrow(fit$U)){ #if user wants to bring out more samples than we have
      return(1);
    }
  
    
    if(omicschoice == 1){
        getting <- lm(fit$U[,jointpred] ~ fit$Tt[,jointpred])
        getting2 <- lm(fit$U[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])] ~
                        fit$Tt[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])])
        tableans <- summary(getting2)[4]$coefficients;
        rownames(tableans)[2] <- paste0("Omics1 Jointscore", jointpred,"-top:", tablen);

    }else{
        getting <- lm(fit$Tt[,jointpred] ~ fit$U[,jointpred])
        getting2 <- lm(fit$Tt[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])] ~
                        fit$U[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])])
        tableans <- summary(getting2)[4]$coefficients;
        rownames(tableans)[2] <- paste0("Omics2 Jointscore", jointpred, "-top:", tablen);     
    }

    tableans <- round(tableans, 5);

    Cairo(file=imgNm, width=8, height=1, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();
    
    return(1);
}


#Main - plot6 (statistics string) update
PerformQualityStatUpdate<-function(omicsChoice, tableN, jointPred, jointb){

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    jointpred <- as.integer(jointPred); #what joint component to compare
    joint <- as.integer(jointb); #joint component limit
    
    if(joint < jointpred){ #if user wants to compare joint outside of the joint dimensions
        return(1);
    } 
    
    fit <- dataSet$fit;
    if(tablen > nrow(fit$U)){ #if user wants to bring out more samples than we have
      return(1);
    }
    
    if(omicschoice == 1){
        getting <- lm(fit$U[,jointpred] ~ fit$Tt[,jointpred])
        getting2 <- lm(fit$U[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])] ~
                        fit$Tt[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])])
        
        stat.info <- paste("UPDATE!!! Adjusted R^2: ", 
                           signif(summary(getting2)[9]$adj.r.squared, 5),
                           " on newly selected samples: ", tablen, sep=""); 

    }else{
        getting <- lm(fit$Tt[,jointpred] ~ fit$U[,jointpred])
        getting2 <- lm(fit$Tt[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])] ~
                        fit$U[,jointpred][(order(cooks.distance(getting), decreasing=F)[1:tablen])])
        
        stat.info <- paste("UPDATE!!! Adjusted p-value: ", 
                           signif(summary(getting2)[9]$adj.r.squared, 5),
                           " on newly selected samples: ", tablen, sep="");  
    }    
    analSet$quality.stat.info.update <- stat.info;
    analSet<<-analSet;

    return(stat.info);
}



#Main - plot7
PlotMainplot7 <- function(imgName, dpi, format, jointb){
    library(plotrix);
    
    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    
    joint <- as.integer(jointb); #joint component limit
    
    if(joint == 1){ #if there is only one joint component, it does not work!!!
        return(1);
    } 
    
    fit <- dataSet$fit;
    #The max of the length can be sqrt(2), so i normalized by dividing it. 
    cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,1], 
                                pc2 = t(cor(fit$Tt, data1[,-1]))[,2], 
                    length = sqrt((t(cor(fit$Tt, data1[,-1]))[,1])^2 + (t(cor(fit$Tt, data1[,-1]))[,2])^2)/sqrt(2)
                      )
                    )

    if(dim(cors)[1] > 200){
        tt <- sample(nrow(cors), 200);
        cors <- cors[tt,,drop=FALSE]
    }
    
    
    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
      plot(cors$pc1, cors$pc2, 
       main = "Correlations - variables & joint scores (Omics1)",
       xlab = "Joint score1", ylab = "Joint score2", 
       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
      abline(h = 0, v = 0, lty = 1)
      arrows(0, 0, cors$pc1, cors$pc2,
             length = 0.07, angle = 30, code = 3)
      draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
    dev.off();
    
    return(1);
}


#Main - plot7 table update
PlotCorVarTable <- function(imgName, dpi, format, omicsChoice, tableN, xaxisjoint, yaxisjoint, option, jointb){
    
    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");


    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    
    
    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    xaxisjoint <- as.integer(xaxisjoint); #what joint component on xaxis
    yaxisjoint <- as.integer(yaxisjoint); #what joint component on yaxis
    
    option <- as.integer(option); 
    #options: 1. n variables of Euclidean dist  2. n absolute variables of both X,Y-axis comp   3. n absolute variables of X-axis comp   4. n absolute variables of Y-axis comp  
    joint <- as.integer(jointb); #joint component limit
    
    if( (joint == 1) | (xaxisjoint > joint) | (yaxisjoint > joint)){ #if there is only one joint component, it does not work!!!
        return(1);
    } 
    
    fit <- dataSet$fit;
      
    if(omicschoice == 1){
        if(tablen > dim(data1[,-1])[2]){ #if user wants to bring out more samples than we have
          return(1);
        }
        
        if(option == 1){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint], 
                      length = sqrt((t(cor(fit$Tt, data1[,-1]))[,xaxisjoint])^2 +
                                      (t(cor(fit$Tt, data1[,-1]))[,yaxisjoint])^2)/sqrt(2)
                        )
                      )
          
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Euclidean")
          
          cors <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- round(cors, 5)
        }else if(option == 2){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.x = abs(cors$pc1), abs.y = abs(cors$pc2)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors3 <- cors[order(cors[,4], decreasing = T)[1:tablen],]
          cors <- as.data.frame(rbind(cors2, cors3))
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:X", "Abs value:Y")
          cors <- round(cors, 5)
        }else if(option == 3){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.x = abs(cors$pc1)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- cors2
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:X")
          
          cors <- round(cors, 5)
        }else{
         cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.y = abs(cors$pc2)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- cors2
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:Y")
          
          cors <- round(cors, 5)
        }
      tableans <- cors;
    }else{
       if(tablen > dim(data2[,-1])[2]){ #if user wants to bring out more samples than we have
          return(1);
        }
        
        if(option == 1){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint], 
                      length = sqrt((t(cor(fit$U, data2[,-1]))[,xaxisjoint])^2 +
                                      (t(cor(fit$U, data2[,-1]))[,yaxisjoint])^2)/sqrt(2)
                        )
                      )
          
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Euclidean")
          
          cors <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- round(cors, 5)
        }else if(option == 2){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.x = abs(cors$pc1), abs.y = abs(cors$pc2)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors3 <- cors[order(cors[,4], decreasing = T)[1:tablen],]
          cors <- as.data.frame(rbind(cors2, cors3))
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:X", "Abs value:Y")
          cors <- round(cors, 5)
        }else if(option == 3){
          cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.x = abs(cors$pc1)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- cors2
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:X")
          
          cors <- round(cors, 5)
        }else{
         cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                  pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                        )
                      )
          
          cors <- as.data.frame(cbind(cors, abs.y = abs(cors$pc2)
                    )
                  )
          
          cors2 <- cors[order(cors[,3], decreasing = T)[1:tablen],]
          cors <- cors2
          colnames(cors) <- c(paste0("Score", xaxisjoint),
                              paste0("Score", yaxisjoint),
                              "Abs value:Y")
          
          cors <- round(cors, 5)
        }
      tableans <- cors;
    }
    
    Cairo(file=imgNm, width=10, height=nrow(tableans)/3,
          type=format, bg="white", dpi=dpi, unit="in");
    grid.table(tableans);
    dev.off();  
    
    return(1);
}






#Main - plot7 plot update
PlotMainplotupdate7 <- function(imgName, dpi, format, omicsChoice, tableN, xaxisjoint, yaxisjoint, option, visualChoice, jointb){
    
    dpi <- as.numeric(dpi);
    imgNm <- paste(imgName, ".", format, sep="");
  
    #set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");


    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    

    omicschoice <- as.integer(omicsChoice); #1: omics1, 2: omics2
    tablen <- as.integer(tableN); #how many top variables    
    xaxisjoint <- as.integer(xaxisjoint); #what joint component on xaxis
    yaxisjoint <- as.integer(yaxisjoint); #what joint component on yaxis
    visual <- as.integer(visualChoice); #1. random 200 including top n 2. only top n
    option <- as.integer(option); 
    #options: 1. n variables of Euclidean dist  2. n absolute variables of both X,Y-axis comp   3. n absolute variables of X-axis comp   4. n absolute variables of Y-axis comp  
    joint <- as.integer(jointb); #joint component limit
    
    
    if( (joint == 1) | (xaxisjoint > joint) | (yaxisjoint > joint)){ #if there is only one joint component, it does not work!!!
        return(1);
    } 
    
    fit <- dataSet$fit;
    
    set.seed(sample.int(28051968, 1)) #for random samplings!!!
    if(omicschoice == 1){
        if(tablen > dim(data1[,-1])[2]){ #if user wants to bring out more samples than we have
          return(1);
        }
        if(visual == 1){
            if(option == 1){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint], 
                          length = sqrt((t(cor(fit$Tt, data1[,-1]))[,xaxisjoint])^2 +
                                          (t(cor(fit$Tt, data1[,-1]))[,yaxisjoint])^2)/sqrt(2)
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                    lens <- (1:nrow(cors))[-(order(cors[,3], decreasing = T)[1:tablen])]
                    tt <- sample(lens, 190);
                    cors <- cors[c(order(cors[,3], decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");    
                    plot(cors$pc1, cors$pc2, 
                    main = "Correlations - variables & joint scores (Omics1)",
                    xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                    xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(cors$length, decreasing = T)[1:tablen]],
                         cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(cors$length, decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();     
              }

          }else if(option == 2){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )
              xcomplist <- order(abs(cors$pc1), decreasing = T)[1:tablen]
              ycomplist <- order(abs(cors$pc2), decreasing = T)[1:tablen]
              intersec <- intersect(xcomplist, ycomplist)
              if(length(intersec) == 0){
                xcomplist <- xcomplist 
                ycomplist <- ycomplist
              }else{
                xcomplist <- xcomplist[-(match(intersec, xcomplist))]
                ycomplist <- ycomplist[-(match(intersec, ycomplist))]
              }
              
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-c(intersec, xcomplist, ycomplist)]
                  tt <- sample(lens, (200 - length(c(intersec, xcomplist, ycomplist))));
                  cors <- cors[c(c(intersec, xcomplist, ycomplist) ,tt),,drop=FALSE]
              
              
                  if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         cors$pc2[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         labels = rownames(cors)[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:length(intersec)], #intersected
                         cors$pc2[1:length(intersec)],
                         labels = rownames(cors)[1:length(intersec)],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))], #x axis score
                         cors$pc2[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                         labels = rownames(cors)[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                         cors$pc2[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist)) ],
                         labels = rownames(cors)[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }else{
                  if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[xcomplist],
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[intersec], #intersected
                         cors$pc2[intersec],
                         labels = rownames(cors)[intersec],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[xcomplist], #x axis score
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }
          }else if(option == 3){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )

              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc1), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc1), decreasing = T)[1:tablen] ,tt),,drop=FALSE]
                  
                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1, cors$pc2,
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();
              }
              
              
          }else{
             cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc2), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc2), decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1, cors$pc2,
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();    
              }
              
              
          }
      }else{ #only top n
          if(option == 1){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint], 
                          length = sqrt((t(cor(fit$Tt, data1[,-1]))[,xaxisjoint])^2 +
                                          (t(cor(fit$Tt, data1[,-1]))[,yaxisjoint])^2)/sqrt(2)
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(cors[,3], decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(cors[,3], decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[order(cors$length, decreasing = T)[1:tablen]], cors$pc2[order(cors$length, decreasing = T)[1:tablen]], 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[order(cors$length, decreasing = T)[1:tablen]], cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(cors$length, decreasing = T)[1:tablen]],
                         cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(cors$length, decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }
          }else if(option == 2){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )
              xcomplist <- order(abs(cors$pc1), decreasing = T)[1:tablen]
              ycomplist <- order(abs(cors$pc2), decreasing = T)[1:tablen]
              intersec <- intersect(xcomplist, ycomplist)
              if(length(intersec) == 0){
                xcomplist <- xcomplist 
                ycomplist <- ycomplist
              }else{
                xcomplist <- xcomplist[-(match(intersec, xcomplist))]
                ycomplist <- ycomplist[-(match(intersec, ycomplist))]
              }
              
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-c(intersec, xcomplist, ycomplist)]
                  tt <- sample(lens, (200 - length(c(intersec, xcomplist, ycomplist))));
                  cors <- cors[c(c(intersec, xcomplist, ycomplist) ,tt),,drop=FALSE]
                  
                  if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:as.numeric(tablen+tablen)], cors$pc2[1:as.numeric(tablen+tablen)], 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:as.numeric(tablen+tablen)], cors$pc2[1:as.numeric(tablen+tablen)],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         cors$pc2[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         labels = rownames(cors)[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))], 
                         cors$pc2[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))], 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           cors$pc2[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:length(intersec)], #intersected
                         cors$pc2[1:length(intersec)],
                         labels = rownames(cors)[1:length(intersec)],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))], #x axis score
                         cors$pc2[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                         labels = rownames(cors)[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                         cors$pc2[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist)) ],
                         labels = rownames(cors)[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }else{
                 if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[c(xcomplist, ycomplist)], cors$pc2[c(xcomplist, ycomplist)], 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[c(xcomplist, ycomplist)], cors$pc2[c(xcomplist, ycomplist)],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[xcomplist],
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[c(intersec, xcomplist, ycomplist)], 
                         cors$pc2[c(intersec, xcomplist, ycomplist)], 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[c(intersec, xcomplist, ycomplist)],
                           cors$pc2[c(intersec, xcomplist, ycomplist)],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[intersec], #intersected
                         cors$pc2[intersec],
                         labels = rownames(cors)[intersec],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[xcomplist], #x axis score
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }
          }else if(option == 3){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )

              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc1), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc1), decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                 Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]], 
                       main = "Correlations - variables & joint scores (Omics1)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();   

              }
          }else{
             cors <- as.data.frame(cbind(pc1 = t(cor(fit$Tt, data1[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$Tt, data1[,-1]))[,yaxisjoint]
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc2), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc2), decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]], 
                         main = "Correlations - variables & joint scores (Omics1)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off(); 
                }    
          }
      }
    }else{
       if(tablen > dim(data2[,-1])[2]){ #if user wants to bring out more samples than we have
          return(1);
        }
        if(visual == 1){
            if(option == 1){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint], 
                          length = sqrt((t(cor(fit$U, data2[,-1]))[,xaxisjoint])^2 +
                                          (t(cor(fit$U, data2[,-1]))[,yaxisjoint])^2)/sqrt(2)
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(cors[,3], decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(cors[,3], decreasing = T)[1:tablen] ,tt),,drop=FALSE]
              
              
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
               }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");    
                    plot(cors$pc1, cors$pc2, 
                    main = "Correlations - variables & joint scores (Omics2)",
                    xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                    xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(cors$length, decreasing = T)[1:tablen]],
                         cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(cors$length, decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();     
              }     

          }else if(option == 2){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )
              xcomplist <- order(abs(cors$pc1), decreasing = T)[1:tablen]
              ycomplist <- order(abs(cors$pc2), decreasing = T)[1:tablen]
              intersec <- intersect(xcomplist, ycomplist)
              if(length(intersec) == 0){
                xcomplist <- xcomplist 
                ycomplist <- ycomplist
              }else{
                xcomplist <- xcomplist[-(match(intersec, xcomplist))]
                ycomplist <- ycomplist[-(match(intersec, ycomplist))]
              }
              
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-c(intersec, xcomplist, ycomplist)]
                  tt <- sample(lens, (200 - length(c(intersec, xcomplist, ycomplist))));
                  cors <- cors[c(c(intersec, xcomplist, ycomplist) ,tt),,drop=FALSE]

                    if(length(intersec) == 0){
                      Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                      plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                      abline(h = 0, v = 0, lty = 1)
                      arrows(0, 0, cors$pc1, cors$pc2,
                             length = 0.07, angle = 30, code = 3)
                      draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                      text(cors$pc1[1:tablen],
                           cors$pc2[1:tablen],
                           labels = rownames(cors)[1:tablen],
                           pos = 2,
                           col = "blue", cex = 0.5)
                      text(cors$pc1[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           cors$pc2[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           labels = rownames(cors)[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           pos = 2,
                           col = "red", cex = 0.5)
                      dev.off();
                    }else{
                      Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                      plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                      abline(h = 0, v = 0, lty = 1)
                      arrows(0, 0, cors$pc1, cors$pc2,
                             length = 0.07, angle = 30, code = 3)
                      draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                      text(cors$pc1[1:length(intersec)], #intersected
                           cors$pc2[1:length(intersec)],
                           labels = rownames(cors)[1:length(intersec)],
                           pos = 2,
                           col = "green", cex = 0.5)
                      text(cors$pc1[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))], #x axis score
                           cors$pc2[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                           labels = rownames(cors)[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                           pos = 2,
                           col = "blue", cex = 0.5)
                      text(cors$pc1[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           cors$pc2[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist)) ],
                           labels = rownames(cors)[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           pos = 2,
                           col = "red", cex = 0.5)
                      dev.off();
                    }
                }else{
                  if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[xcomplist],
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[intersec], #intersected
                         cors$pc2[intersec],
                         labels = rownames(cors)[intersec],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[xcomplist], #x axis score
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }
          }else if(option == 3){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )

              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc1), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc1), decreasing = T)[1:tablen] ,tt),,drop=FALSE]
              
              
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
               }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1, cors$pc2,
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();
              }

          }else{
             cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc2), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc2), decreasing = T)[1:tablen] ,tt),,drop=FALSE]
              
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1, cors$pc2, 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1, cors$pc2,
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1, cors$pc2, 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1, cors$pc2,
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();    
              }
          }
      }else{ #only top n
          if(option == 1){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint], 
                          length = sqrt((t(cor(fit$U, data2[,-1]))[,xaxisjoint])^2 +
                                          (t(cor(fit$U, data2[,-1]))[,yaxisjoint])^2)/sqrt(2)
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(cors[,3], decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(cors[,3], decreasing = T)[1:tablen] ,tt),,drop=FALSE]

                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[order(cors$length, decreasing = T)[1:tablen]], cors$pc2[order(cors$length, decreasing = T)[1:tablen]], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[order(cors$length, decreasing = T)[1:tablen]], cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(cors$length, decreasing = T)[1:tablen]],
                         cors$pc2[order(cors$length, decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(cors$length, decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
              }      
              
          }else if(option == 2){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )
              xcomplist <- order(abs(cors$pc1), decreasing = T)[1:tablen]
              ycomplist <- order(abs(cors$pc2), decreasing = T)[1:tablen]
              intersec <- intersect(xcomplist, ycomplist)
              if(length(intersec) == 0){
                xcomplist <- xcomplist 
                ycomplist <- ycomplist
              }else{
                xcomplist <- xcomplist[-(match(intersec, xcomplist))]
                ycomplist <- ycomplist[-(match(intersec, ycomplist))]
              }
              
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-c(intersec, xcomplist, ycomplist)]
                  tt <- sample(lens, (200 - length(c(intersec, xcomplist, ycomplist))));
                  cors <- cors[c(c(intersec, xcomplist, ycomplist) ,tt),,drop=FALSE]
              
                    if(length(intersec) == 0){
                      Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                      plot(cors$pc1[1:as.numeric(tablen+tablen)], cors$pc2[1:as.numeric(tablen+tablen)], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                      abline(h = 0, v = 0, lty = 1)
                      arrows(0, 0, cors$pc1[1:as.numeric(tablen+tablen)], cors$pc2[1:as.numeric(tablen+tablen)],
                             length = 0.07, angle = 30, code = 3)
                      draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                      text(cors$pc1[1:tablen],
                           cors$pc2[1:tablen],
                           labels = rownames(cors)[1:tablen],
                           pos = 2,
                           col = "blue", cex = 0.5)
                      text(cors$pc1[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           cors$pc2[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           labels = rownames(cors)[as.numeric(tablen+1):as.numeric(tablen+tablen)],
                           pos = 2,
                           col = "red", cex = 0.5)
                      dev.off();
                    }else{
                      Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                      plot(cors$pc1[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))], 
                           cors$pc2[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                      abline(h = 0, v = 0, lty = 1)
                      arrows(0, 0, cors$pc1[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                             cors$pc2[1:as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                             length = 0.07, angle = 30, code = 3)
                      draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                      text(cors$pc1[1:length(intersec)], #intersected
                           cors$pc2[1:length(intersec)],
                           labels = rownames(cors)[1:length(intersec)],
                           pos = 2,
                           col = "green", cex = 0.5)
                      text(cors$pc1[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))], #x axis score
                           cors$pc2[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                           labels = rownames(cors)[as.numeric(length(intersec)+1): as.numeric(length(intersec)+length(xcomplist))],
                           pos = 2,
                           col = "blue", cex = 0.5)
                      text(cors$pc1[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           cors$pc2[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist)) ],
                           labels = rownames(cors)[as.numeric(length(intersec)+length(xcomplist)+1):as.numeric(length(intersec)+length(xcomplist)+length(ycomplist))],
                           pos = 2,
                           col = "red", cex = 0.5)
                      dev.off();
                    }
                }else{
                 if(length(intersec) == 0){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[c(xcomplist, ycomplist)], cors$pc2[c(xcomplist, ycomplist)], 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[c(xcomplist, ycomplist)], cors$pc2[c(xcomplist, ycomplist)],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[xcomplist],
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[c(intersec, xcomplist, ycomplist)], 
                         cors$pc2[c(intersec, xcomplist, ycomplist)], 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[c(intersec, xcomplist, ycomplist)],
                           cors$pc2[c(intersec, xcomplist, ycomplist)],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[intersec], #intersected
                         cors$pc2[intersec],
                         labels = rownames(cors)[intersec],
                         pos = 2,
                         col = "green", cex = 0.5)
                    text(cors$pc1[xcomplist], #x axis score
                         cors$pc2[xcomplist],
                         labels = rownames(cors)[xcomplist],
                         pos = 2,
                         col = "blue", cex = 0.5)
                    text(cors$pc1[ycomplist],
                         cors$pc2[ycomplist],
                         labels = rownames(cors)[ycomplist],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
                  }
              }    
          }else if(option == 3){
              cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )

              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc1), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc1), decreasing = T)[1:tablen] ,tt),,drop=FALSE]
              
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
               }else{
                 Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]], 
                       main = "Correlations - variables & joint scores (Omics2)",
                       xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                       xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                  abline(h = 0, v = 0, lty = 1)
                  arrows(0, 0, cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                         length = 0.07, angle = 30, code = 3)
                  draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                  text(cors$pc1[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       cors$pc2[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       labels = rownames(cors)[order(abs(cors$pc1), decreasing = T)[1:tablen]],
                       pos = 2,
                       col = "red", cex = 0.5)
                  dev.off();   

              }
     
          }else{
             cors <- as.data.frame(cbind(pc1 = t(cor(fit$U, data2[,-1]))[,xaxisjoint], 
                                      pc2 = t(cor(fit$U, data2[,-1]))[,yaxisjoint]
                            )
                          )
              
              if(dim(cors)[1] > 200){ #for visual purpose
                  lens <- (1:nrow(cors))[-(order(abs(cors$pc2), decreasing = T)[1:tablen])]
                  tt <- sample(lens, 190);
                  cors <- cors[c(order(abs(cors$pc2), decreasing = T)[1:tablen] ,tt),,drop=FALSE]
              
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[1:tablen], cors$pc2[1:tablen], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[1:tablen], cors$pc2[1:tablen],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[1:tablen],
                         cors$pc2[1:tablen],
                         labels = rownames(cors)[1:tablen],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off();
               }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]], 
                         main = "Correlations - variables & joint scores (Omics2)",
                         xlab = paste0("Joint score", xaxisjoint), ylab = paste0("Joint score", yaxisjoint), 
                         xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0), asp = 1)
                    abline(h = 0, v = 0, lty = 1)
                    arrows(0, 0, cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]], cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                           length = 0.07, angle = 30, code = 3)
                    draw.circle(0, 0, 1, border = "gray", lty = 1, lwd = 1)
                    text(cors$pc1[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         cors$pc2[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         labels = rownames(cors)[order(abs(cors$pc2), decreasing = T)[1:tablen]],
                         pos = 2,
                         col = "red", cex = 0.5)
                    dev.off(); 
                }       
          }
      }
    }
    return(1);
}    



#Main - plot8
PlotMainplot8 <- function(tableName, dpi, format){
    library(corrplot);
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- dataSet$norm1;
    #save.image("correlation.RData")
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;
    

    if(is.character(levels(data1[,1]))){ #if the factor has levels of characters
        xx <- factor(as.numeric(data1[,1]), levels = as.numeric((1:length(levels(data1[,1])))))
        xx <- as.numeric(as.character(xx))
        xx <- as.data.frame(xx)
        colnames(xx) <- colnames(data1)[1]
        
        M1 <- cor(fit$Tt, xx)
        x <- c()
        for(i in 1:nrow(M1)){
          x[i] <- paste0("Score", i, "-omics1(", model_summary$Omics1[i], "%)")
        }
        rownames(M1) <- x
        
        M2 <- cor(fit$U, xx)
        y <- c()
        for(i in 1:nrow(M2)){
          y[i] <- paste0("Score", i, "-omics2(", model_summary$Omics2[i], "%)")
        }
        rownames(M2) <- y
    }else{
        xx <- as.numeric(as.character(data1[,1]))
        xx <- as.data.frame(xx)
        colnames(xx) <- colnames(data1)[1]
        
        M1 <- cor(fit$Tt, xx)
        x <- c()
        for(i in 1:nrow(M1)){
          x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
        }
        rownames(M1) <- x
        
        M2 <- cor(fit$U, xx)
        y <- c()
        for(i in 1:nrow(M2)){
          y[i] <- paste0("C", i, "-protein(", model_summary$Omics2[i], "%)")
        }
        rownames(M2) <- y
    }
    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    
    Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
    corrplot(rbind(M1, M2), method = "number", col = col1(100),
             tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
    dev.off();
    #save.image("correlation.RData")
    return(1);
}



#Main - plot9
PlotMainplot9 <- function(tableName, dpi, format){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;

    if(ncol(fit$Tt) == 1){
        return(1);
    }
    
    jointx <- fit$Tt
    x <- c()
    for(i in 1:ncol(jointx)){
      x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
    }
    colnames(jointx) <- x

    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    
    
    
        
        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number"
                 )
        dev.off();
   
    
    return(1);
}

#Main - plot10
PlotMainplot10 <- function(tableName, dpi, format){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;
    
    jointy <- fit$U
    y <- c()
    for(i in 1:ncol(jointy)){
      y[i] <- paste0("C", i, "-omics2(", model_summary$Omics2[i], "%)")
    }
    colnames(jointy) <- y

    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    

    if(ncol(fit$U) == 1){
        return(1);
    }                       

    
        
        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number"
                 )
        dev.off();
    
    return(1);
}



#Main - plot11
PlotMainplot11 <- function(tableName, dpi, format){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;
    
    jointx <- fit$Tt
    jointy <- fit$U
    x <- c()
    for(i in 1:ncol(jointx)){
      x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
    }
    colnames(jointx) <- x
    
    y <- c()
    for(i in 1:ncol(jointy)){
      y[i] <- paste0("C", i, "-omics2(", model_summary$Omics2[i], "%)")
    }
    colnames(jointy) <- y

    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    
    
    if(ncol(fit$Tt)==1){
        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx, jointy), 
                 insig = "pch", addrect = 2, col = col1(100), method = "number"
                 )
        dev.off();
    }else{
        Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx, jointy), 
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number"
                 )
        dev.off();
    }

    return(1);
}






#Main - plot8 update
PlotMainplotupdate8 <- function(tableName, dpi, format, color, shape){
    library(OmicsPLS);
    library(corrplot);
    library(RColorBrewer);
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    color <- as.integer(color);
    shape <- as.integer(shape);
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;

    if(is.character(levels(data1[,1]))){ #if the factor has levels of characters
        xx <- factor(as.numeric(data1[,1]), levels = as.numeric((1:length(levels(data1[,1])))))
        xx <- as.numeric(as.character(xx))
        xx <- as.data.frame(xx)
        colnames(xx) <- colnames(data1)[1]
        
        M1 <- cor(fit$Tt, xx)
        x <- c()
        for(i in 1:nrow(M1)){
          x[i] <- paste0("Score", i, "-omics1(", model_summary$Omics1[i], "%)")
        }
        rownames(M1) <- x
        
        M2 <- cor(fit$U, xx)
        y <- c()
        for(i in 1:nrow(M2)){
          y[i] <- paste0("Score", i, "-omics2(", model_summary$Omics2[i], "%)")
        }
        rownames(M2) <- y
    }else{
        xx <- as.numeric(as.character(data1[,1]))
        xx <- as.data.frame(xx)
        colnames(xx) <- colnames(data1)[1]
        
        M1 <- cor(fit$Tt, as.numeric(as.character(xx)))
        x <- c()
        for(i in 1:nrow(M1)){
          x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
        }
        rownames(M1) <- x
        
        M2 <- cor(fit$U, as.numeric(as.character(xx)))
        y <- c()
        for(i in 1:nrow(M2)){
          y[i] <- paste0("C", i, "-protein(", model_summary$Omics2[i], "%)")
        }
        rownames(M2) <- y
    }
    
    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
    col3 <- colorRampPalette(c("red", "white", "blue")) 
    col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                               "cyan", "#007FFF", "blue", "#00007F"))
    whiteblack <- c("white", "black")

    if(shape == 1){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "number",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else if(shape == 2){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "circle",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else if(shape == 3){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "square",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else if(shape == 4){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "color",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else if(shape == 5){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "ellipse",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else if(shape == 6){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "shade",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }else{
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie", col = col1(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie", col = col2(50),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = col3(20),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = col4(10),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = whiteblack, bg = "gold2",
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = heat.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = terrain.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = cm.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = gray.colors(100),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = brewer.pal(n = 8, name = "RdBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = brewer.pal(n = 8, name = "RdYlBu"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(rbind(M1, M2), method = "pie",  col = brewer.pal(n = 8, name = "PuOr"),
                 tl.col = "blue", cl.pos = "r", cl.ratio = 0.5)
        dev.off();
      }
    }
    
    return(1);
}






#Main - plot9 update
PlotMainplotupdate9 <- function(tableName, dpi, format, color, shape){
    library(corrplot);
    library(RColorBrewer);
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    color <- as.integer(color);
    shape <- as.integer(shape);
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;

    if(ncol(fit$Tt) == 1){
        return(1);
    }

    
    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
    col3 <- colorRampPalette(c("red", "white", "blue")) 
    col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                               "cyan", "#007FFF", "blue", "#00007F"))
    whiteblack <- c("white", "black")

    jointx <- fit$Tt
    x <- c()
    for(i in 1:ncol(jointx)){
      x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
    }
    colnames(jointx) <- x


    
    if(shape == 1){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "number")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "number")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "number")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "number")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "number")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "number")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "number")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "number")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "number")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "number")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "number")
        dev.off();
      }
    }else if(shape == 2){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "circle")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "circle")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "circle")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "circle")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "circle")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "circle")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "circle")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "circle")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "circle")
        dev.off();
      }
    }else if(shape == 3){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "square")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "square")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "square")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "square")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "square")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "square")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "square")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "square")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "square")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "square")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "square")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "square")
        dev.off();
      }
    }else if(shape == 4){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "color")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "color")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "color")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "color")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "color")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "color")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "color")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "color")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "color")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "color")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "color")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "color")
        dev.off();
      }
    }else if(shape == 5){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "ellipse")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "ellipse")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "ellipse")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "ellipse")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "ellipse")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "ellipse")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "ellipse")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "ellipse")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "ellipse")
        dev.off();
      }
    }else if(shape == 6){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "shade")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "shade")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "shade")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "shade")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "shade")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "shade")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),  
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), p.mat = res1$p, 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "shade")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "shade")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx), 
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "shade")
        dev.off();
      }
    }else{
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "pie")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "pie")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "pie")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "pie")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "pie")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "pie")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "pie")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "pie")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointx),   
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "pie")
        dev.off();
      }
    }
    
    return(1);
}






#Main - plot10 update
PlotMainplotupdate10 <- function(tableName, dpi, format, color, shape){
    library(corrplot);
    library(RColorBrewer);
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    color <- as.integer(color);
    shape <- as.integer(shape);
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;
    
    if(ncol(fit$U) == 1){
        return(1);
    }
    
    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
    col3 <- colorRampPalette(c("red", "white", "blue")) 
    col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                               "cyan", "#007FFF", "blue", "#00007F"))
    whiteblack <- c("white", "black")

    jointy <- fit$U
    y <- c()
    for(i in 1:ncol(jointy)){
      y[i] <- paste0("C", i, "-omics2(", model_summary$Omics2[i], "%)")
    }
    colnames(jointy) <- y

    
    
    if(shape == 1){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "number")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "number")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "number")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "number")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "number")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "number")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "number")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "number")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "number")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "number")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),    
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "number")
        dev.off();
      }
    }else if(shape == 2){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "circle")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "circle")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "circle")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "circle")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "circle")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "circle")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "circle")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "circle")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "circle")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "circle")
        dev.off();
      }
    }else if(shape == 3){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "square")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "square")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "square")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "square")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "square")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "square")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "square")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "square")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "square")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "square")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "square")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "square")
        dev.off();
      }
    }else if(shape == 4){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "color")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "color")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "color")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "color")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "color")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "color")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "color")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "color")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "color")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "color")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "color")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "color")
        dev.off();
      }
    }else if(shape == 5){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "ellipse")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "ellipse")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "ellipse")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "ellipse")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "ellipse")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "ellipse")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "ellipse")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "ellipse")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "ellipse")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "ellipse")
        dev.off();
      }
    }else if(shape == 6){
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "shade")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "shade")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "shade")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "shade")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "shade")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "shade")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "shade")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "shade")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "shade")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "shade")
        dev.off();
      }
    }else{
      if(color == 1){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "pie")
        dev.off();
      }else if(color == 2){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "pie")
        dev.off();
      }else if(color == 3){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "pie")
        dev.off();
      }else if(color == 4){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "pie")
        dev.off();
      }else if(color == 5){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "pie")
        dev.off();
      }else if(color == 6){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
             method = "pie")
        dev.off();
      }else if(color == 7){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 8){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 9){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
             method = "pie")
        dev.off();
      }else if(color == 10){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
             method = "pie")
        dev.off();
      }else if(color==11){
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
             method = "pie")
        dev.off();
      }else{
        Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
        corrplot(cor(jointy),     
             order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
             method = "pie")
        dev.off();
      }
    }
    
    return(1);
}


#Main - plot11 update
PlotMainplotupdate11 <- function(tableName, dpi, format, color, shape){
  
    dpi <- as.numeric(dpi);
    imgNm <- paste(tableName, ".", format, sep="");
  
    set.seed(28051968);
    options(scipen = 999);
    data1 <- readRDS("norm.data1");
    color <- as.integer(color);
    shape <- as.integer(shape);
    
    
    model_summary <- imgSet$o2pls.main2 #model_summary
    fit <- dataSet$fit;

    
    col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
    col3 <- colorRampPalette(c("red", "white", "blue")) 
    col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                               "cyan", "#007FFF", "blue", "#00007F"))
    whiteblack <- c("white", "black")

    jointx <- fit$Tt
    jointy <- fit$U
    x <- c()
    for(i in 1:ncol(jointx)){
      x[i] <- paste0("C", i, "-omics1(", model_summary$Omics1[i], "%)")
    }
    colnames(jointx) <- x
    
    y <- c()
    for(i in 1:ncol(jointy)){
      y[i] <- paste0("C", i, "-omics2(", model_summary$Omics2[i], "%)")
    }
    colnames(jointy) <- y
    
    
        if(ncol(fit$Tt) == 1){
           if(shape == 1){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "number")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "number")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "number")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "number")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "number")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), method = "number")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "number")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "number")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "number")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "number")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "number")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "number")
                dev.off();
              }
            }else if(shape == 2){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "circle")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "circle")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "circle")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "circle")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "circle")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), method = "circle")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "circle")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "circle")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "circle")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "circle")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "circle")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "circle")
                dev.off();
              }
            }else if(shape == 3){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "square")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "square")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "square")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "square")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "square")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), method = "square")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "square")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "square")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "square")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "square")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "square")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "square")
                dev.off();
              }
            }else if(shape == 4){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "color")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "color")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "color")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "color")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "color")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), method = "color")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "color")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "color")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "color")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "color")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "color")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "color")
                dev.off();
              }
            }else if(shape == 5){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "ellipse")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "ellipse")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "ellipse")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "ellipse")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "ellipse")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), 
                     method = "ellipse")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "ellipse")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "ellipse")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "ellipse")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "ellipse")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "ellipse")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "ellipse")
                dev.off();
              }
            }else if(shape == 6){
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col1(100), method = "shade")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col2(50), method = "shade")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col3(20), method = "shade")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = col4(10), method = "shade")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = whiteblack, method = "shade")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = heat.colors(100), 
                     method = "shade")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "shade")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                       insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "shade")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "shade")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "shade")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "shade")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "shade")
                dev.off();
              }
            }else{
              if(color == 1){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = col1(100), method = "pie")
                dev.off();
              }else if(color == 2){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = col2(50), method = "pie")
                dev.off();
              }else if(color == 3){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = col3(20), method = "pie")
                dev.off();
              }else if(color == 4){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = col4(10), method = "pie")
                dev.off();
              }else if(color == 5){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = whiteblack, method = "pie")
                dev.off();
              }else if(color == 6){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = heat.colors(100), 
                     method = "pie")
                dev.off();
              }else if(color == 7){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = terrain.colors(100),
                     method = "pie")
                dev.off();
              }else if(color == 8){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = cm.colors(100),
                     method = "pie")
                dev.off();
              }else if(color == 9){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = gray.colors(100),
                     method = "pie")
                dev.off();
              }else if(color == 10){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                     method = "pie")
                dev.off();
              }else if(color==11){
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                     method = "pie")
                dev.off();
              }else{
                Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
                corrplot(cor(jointx, jointy),     
                      insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                     method = "pie")
                dev.off();
              }
            }
            return(1);
        }    
        if(shape == 1){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "number")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "number")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "number")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "number")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "number")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "number")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "number")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "number")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "number")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "number")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "number")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "number")
            dev.off();
          }
        }else if(shape == 2){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "circle")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "circle")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "circle")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "circle")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "circle")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "circle")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "circle")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "circle")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "circle")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "circle")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "circle")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "circle")
            dev.off();
          }
        }else if(shape == 3){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "square")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "square")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "square")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "square")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "square")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "square")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "square")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "square")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "square")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "square")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "square")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "square")
            dev.off();
          }
        }else if(shape == 4){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "color")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "color")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "color")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "color")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "color")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), method = "color")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "color")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "color")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "color")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "color")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "color")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "color")
            dev.off();
          }
        }else if(shape == 5){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "ellipse")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "ellipse")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "ellipse")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "ellipse")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "ellipse")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
                 method = "ellipse")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "ellipse")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "ellipse")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "ellipse")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "ellipse")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "ellipse")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "ellipse")
            dev.off();
          }
        }else if(shape == 6){
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "shade")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "shade")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "shade")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "shade")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "shade")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
                 method = "shade")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "shade")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "shade")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "shade")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "shade")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "shade")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "shade")
            dev.off();
          }
        }else{
          if(color == 1){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col1(100), method = "pie")
            dev.off();
          }else if(color == 2){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col2(50), method = "pie")
            dev.off();
          }else if(color == 3){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col3(20), method = "pie")
            dev.off();
          }else if(color == 4){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = col4(10), method = "pie")
            dev.off();
          }else if(color == 5){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = whiteblack, method = "pie")
            dev.off();
          }else if(color == 6){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = heat.colors(100), 
                 method = "pie")
            dev.off();
          }else if(color == 7){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = terrain.colors(100),
                 method = "pie")
            dev.off();
          }else if(color == 8){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = cm.colors(100),
                 method = "pie")
            dev.off();
          }else if(color == 9){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = gray.colors(100),
                 method = "pie")
            dev.off();
          }else if(color == 10){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdBu"),
                 method = "pie")
            dev.off();
          }else if(color==11){
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "RdYlBu"),
                 method = "pie")
            dev.off();
          }else{
            Cairo(file=imgNm, width=8, height=7, type=format, bg="white", dpi=dpi, unit="in");
            corrplot(cor(jointx, jointy),     
                 order = "hclust", insig = "pch", addrect = 2, col = brewer.pal(n = 8, name = "PuOr"),
                 method = "pie")
            dev.off();
          }
        }
        
    return(1);
}




PlotMainplotsummary <- function(tableName, dpi, format, choice, samples, jointcomp, jointb){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
    options(scipen = 999);
    fit <- dataSet$fit;
    model_summary <- imgSet$o2pls.main2 #model_summary
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");

    joint <- as.integer(jointb); #how many joint components users choose to have
    choice <- as.integer(choice); #1: on all samples 2. on selected samples
    samples <- as.integer(samples);#disabled if choice == 1
    jointcomp <- as.integer(jointcomp); #what joint score to show as the summary table


    fircolnm <- colnames(data1)[1];

    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;


    if(length(samples)==1){
      if(samples > dim(fit$Tt)[1]){
        return(1);
      }
    }

    if(jointcomp > joint){
      return(1);
    }

    samplenm <- c(); #sample names
    metanm <- c(); #sample metadata names
    if(choice == 1){ #all samples
          samplenm <- rownames(fit$U)
          metanm <- data1[,1]
    }else{ #select samples
          getting <- lm(fit$U[,jointcomp] ~ fit$Tt[,jointcomp])

          samplenm <- rownames(melt(sort(cooks.distance(getting), decreasing = F)[1:samples]))
          metanm <- data1[order(cooks.distance(getting), decreasing = F)[1:samples],1]
    }
    
    databrick <- data.frame.na(metadata_names = as.character(metanm), 
                               sample_name = samplenm
                               )

    Cairo(file=imgNm, width=7, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
    grid.table(databrick);
    dev.off();

    return(1);
}



PlotMainplotsummary2 <- function(tableName, dpi, format, tableNomics1, tableNomics2, jointcomp, jointb){
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
    options(scipen = 999);
    fit <- dataSet$fit;
    model_summary <- imgSet$o2pls.main2 #model_summary
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    tablenomics1 <- as.integer(tableNomics1); #how many top variables -> "Loading"
    tablenomics2 <- as.integer(tableNomics2); #how many top variables -> "Loading"
    joint <- as.integer(jointb); #how many joint components users choose to have
    jointcomp <- as.integer(jointcomp); #what joint score to show as the summary table

    fircolnm <- colnames(data1)[1];

    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    
    if(tablenomics1 > dim(fit$W.)[1]){
        return(1);
    }
    if(tablenomics2 > dim(fit$C.)[1]){
        return(1);
    }
  
    if(jointcomp > joint){
      return(1);
    }

    sig_featomics1 <- c();
    sig_featomics2 <- c();

    load1x <- fit$W.[,jointcomp]
    load1xjoint <- melt(sort(abs(load1x), decreasing = T)[1:tablenomics1]) #sum of joint

    load1y <- fit$C.[,jointcomp]
    load1yjoint <-melt(sort(abs(load1y), decreasing = T)[1:tablenomics2])

    sig_featomics1 <- load1xjoint
    sig_featomics2 <- load1yjoint
       
    ## data.frames with filled NA's when the number of rows are different 
    data.frame.na <- function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,
      stringsAsFactors = FALSE){
        data.row.names <- if (check.rows && is.null(row.names))
            function(current, new, i) {
                if (is.character(current))
                    new <- as.character(new)
                if (is.character(new))
                    current <- as.character(current)
                if (anyDuplicated(new))
                    return(current)
                if (is.null(current))
                    return(new)
                if (all(current == new) || all(current == ""))
                    return(new)
                stop(gettextf("mismatch of row names in arguments of 'data.frame', item %d",
                    i), domain = NA)
            }
        else function(current, new, i) {
            if (is.null(current)) {
                if (anyDuplicated(new)) {
                    warning("some row.names duplicated: ", paste(which(duplicated(new)),
                      collapse = ","), " --> row.names NOT used")
                    current
                }
                else new
            }
            else current
        }
        object <- as.list(substitute(list(...)))[-1L]
        mrn <- is.null(row.names)
        x <- list(...)
        n <- length(x)
        if (n < 1L) {
            if (!mrn) {
                if (is.object(row.names) || !is.integer(row.names))
                    row.names <- as.character(row.names)
                if (any(is.na(row.names)))
                    stop("row names contain missing values")
                if (anyDuplicated(row.names))
                    stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                      collapse = ", "))
            }
            else row.names <- integer(0L)
            return(structure(list(), names = character(0L), row.names = row.names,
                class = "data.frame"))
        }
        vnames <- names(x)
        if (length(vnames) != n)
            vnames <- character(n)
        no.vn <- !nzchar(vnames)
        vlist <- vnames <- as.list(vnames)
        nrows <- ncols <- integer(n)
        for (i in seq_len(n)) {
            xi <- if (is.character(x[[i]]) || is.list(x[[i]]))
                as.data.frame(x[[i]], optional = TRUE, stringsAsFactors = stringsAsFactors)
            else as.data.frame(x[[i]], optional = TRUE)
            nrows[i] <- .row_names_info(xi)
            ncols[i] <- length(xi)
            namesi <- names(xi)
            if (ncols[i] > 1L) {
                if (length(namesi) == 0L)
                    namesi <- seq_len(ncols[i])
                if (no.vn[i])
                    vnames[[i]] <- namesi
                else vnames[[i]] <- paste(vnames[[i]], namesi, sep = ".")
            }
            else {
                if (length(namesi))
                    vnames[[i]] <- namesi
                else if (no.vn[[i]]) {
                    tmpname <- deparse(object[[i]])[1L]
                    if (substr(tmpname, 1L, 2L) == "I(") {
                      ntmpn <- nchar(tmpname, "c")
                      if (substr(tmpname, ntmpn, ntmpn) == ")")
                        tmpname <- substr(tmpname, 3L, ntmpn - 1L)
                    }
                    vnames[[i]] <- tmpname
                }
            }
            if (missing(row.names) && nrows[i] > 0L) {
                rowsi <- attr(xi, "row.names")
                nc <- nchar(rowsi, allowNA = FALSE)
                nc <- nc[!is.na(nc)]
                if (length(nc) && any(nc))
                    row.names <- data.row.names(row.names, rowsi,
                      i)
            }
            nrows[i] <- abs(nrows[i])
            vlist[[i]] <- xi
        }
        nr <- max(nrows)
        for (i in seq_len(n)[nrows < nr]) {
            xi <- vlist[[i]]
            if (nrows[i] > 0L) {
                xi <- unclass(xi)
                fixed <- TRUE
                for (j in seq_along(xi)) {
                    ### added NA fill to max length/nrow
                    xi1 <- xi[[j]]
                    if (is.vector(xi1) || is.factor(xi1))
                      xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                    else if (is.character(xi1) && class(xi1) == "AsIs")
                      xi[[j]] <- structure(c(xi1, rep(NA, nr - nrows[i])),
                        class = class(xi1))
                    else if (inherits(xi1, "Date") || inherits(xi1,
                      "POSIXct"))
                      xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                    else {
                      fixed <- FALSE
                      break
                    }
                }
                if (fixed) {
                    vlist[[i]] <- xi
                    next
                }
            }
            stop("arguments imply differing number of rows: ", paste(unique(nrows),
                collapse = ", "))
        }
        value <- unlist(vlist, recursive = FALSE, use.names = FALSE)
        vnames <- unlist(vnames[ncols > 0L])
        noname <- !nzchar(vnames)
        if (any(noname))
            vnames[noname] <- paste("Var", seq_along(vnames), sep = ".")[noname]
        if (check.names)
            vnames <- make.names(vnames, unique = TRUE)
        names(value) <- vnames
        if (!mrn) {
            if (length(row.names) == 1L && nr != 1L) {
                if (is.character(row.names))
                    row.names <- match(row.names, vnames, 0L)
                if (length(row.names) != 1L || row.names < 1L ||
                    row.names > length(vnames))
                    stop("row.names should specify one of the variables")
                i <- row.names
                row.names <- value[[i]]
                value <- value[-i]
            }
            else if (!is.null(row.names) && length(row.names) !=
                nr)
                stop("row names supplied are of the wrong length")
        }
        else if (!is.null(row.names) && length(row.names) != nr) {
            warning("row names were found from a short variable and have been discarded")
            row.names <- NULL
        }
        if (is.null(row.names))
            row.names <- .set_row_names(nr)
        else {
            if (is.object(row.names) || !is.integer(row.names))
                row.names <- as.character(row.names)
            if (any(is.na(row.names)))
                stop("row names contain missing values")
            if (anyDuplicated(row.names))
                stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                    collapse = ", "))
        }
        attr(value, "row.names") <- row.names
        attr(value, "class") <- "data.frame"
        value
    }


      databrick <- data.frame.na(feature_omics1 = rownames(sig_featomics1), 
                                 feature_omics2 = rownames(sig_featomics2)
                                 )
                           
      Cairo(file=imgNm, width=8, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
      grid.table(databrick);
      dev.off();

      return(1);
}



#cross-correlation for significant features from two omics
PlotMainplotsummary2cor <- function(tableName, dpi, format, tableNomics1, tableNomics2, jointcomp, jointb, correlationoption){
    library(pheatmap);
  
    dpi <- as.numeric(dpi)
    imgNm <- paste(tableName, ".", format, sep="");
    options(scipen = 999);
    fit <- dataSet$fit;
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    tablenomics1 <- as.integer(tableNomics1); #how many top variables -> "Loading"
    tablenomics2 <- as.integer(tableNomics2); #how many top variables -> "Loading"
    joint <- as.integer(jointb); #how many joint components users choose to have
    jointcomp <- as.integer(jointcomp); #what joint score to show as the summary table
    correlationoption <- as.integer(correlationoption) #1. pearson 2. spearman

    fircolnm <- colnames(data1)[1];

    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
    
    if(tablenomics1 > dim(fit$W.)[1]){
        return(1);
    }
    if(tablenomics2 > dim(fit$C.)[1]){
        return(1);
    }
  
    if(jointcomp > joint){
      return(1);
    }
    load1x <- fit$W.[,jointcomp]
    load1y <- fit$C.[,jointcomp]

    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
    #add 1 since the first column is sample meta-data
    pheatmap(cor(data1[,order(abs(load1x), decreasing = T)[1:tablenomics1]+1], 
        data2[,order(abs(load1y), decreasing = T)[1:tablenomics2]+1]), 
                   fontsize_row = 5,
                   fontsize_col = 5,
                   cluster_cols = F,
                   cluster_rows = F)
    dev.off();

    return(1);
}


#Main - plot11 (statistics string) update
PerformQualityStatUpdate2<-function(samplenumber, jointchoice, jointb, choice){

    samplenumber <- as.integer(samplenumber); #how many top variables    
    jointchoice <- as.integer(jointchoice); #what joint component to compare
    joint <- as.integer(jointb); #joint component limit
    choice <- as.integer(choice); #1. use all 2. choose samples
    
    if(joint < jointchoice){ #if user wants to compare joint outside of the joint dimensions
        return(1);
    } 
    
    fit <- dataSet$fit;
    if(samplenumber > nrow(fit$U)){ #if user wants to bring out more samples than we have
      return(1);
    }
    
    if(choice == 1){ #all samples
        getting <- lm(fit$U[,jointchoice] ~ fit$Tt[,jointchoice])

        stat.info <- paste("UPDATE!!! Adjusted R^2: ", 
                           signif(summary(getting)[9]$adj.r.squared, 5),
                           " on entire samples: ", dim(fit$U)[1], sep=""); 
    }else{
        getting <- lm(fit$U[,jointchoice] ~ fit$Tt[,jointchoice])
        getting2 <- lm(fit$U[,jointchoice][(order(cooks.distance(getting), decreasing=F)[1:samplenumber])] ~
                        fit$Tt[,jointchoice][(order(cooks.distance(getting), decreasing=F)[1:samplenumber])])
        
        stat.info <- paste("UPDATE!!! Adjusted R^2: ", 
                           signif(summary(getting2)[9]$adj.r.squared, 5),
                           " on newly selected samples: ", samplenumber, sep=""); 
    }

    analSet$quality.stat.info.update2 <- stat.info;
    analSet<<-analSet;

    return(stat.info);
}






################################

#Procrustes analysis

PlotPAplot1 <- function(imgName, dpi, format, sizes, types, colors, arrow){
  library(vegan) #procrustes function from MCMCpack library requires the same number of columns....
  library(RColorBrewer)
  options(scipen = 999);
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgName, ".", format, sep="");
  data1 <- readRDS("norm.data1");
  data2 <- readRDS("norm.data2");
  
  types <- as.integer(types) #1. points (default) 2. text 
  sizes <- as.numeric(sizes) #point, text sizes (0 ~ 1)
  colors <- as.integer(colors)
  arrow <- as.integer(arrow) #1. yes 2. no
  
  if(dim(data2)[2] < dim(data1)[2]){
    #standardize the data to means of 0 and standard deviations of 1
    data1 <- decostand(data1[,-1], method = "standardize")
    data2 <- decostand(data2[,-1], method = "standardize")
    
    #PCA based on Euclidean distance
    data1 <- rda(data1)
    data2 <- rda(data2)
    
    procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
  }else{
    #standardize the data to means of 0 and standard deviations of 1
    data1 <- decostand(data1[,-1], method = "standardize")
    data2 <- decostand(data2[,-1], method = "standardize")
    
    #PCA based on Euclidean distance
    data1 <- rda(data1)
    data2 <- rda(data2)
    
    procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
  }

  if(types == 1){ # point
        if(arrow == 1){ #arrow exists
                if(colors == 1){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "red")
                    dev.off();
                }else if(colors == 2){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 9, name = "YlOrBr")[5], len=0.09, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                    dev.off();
                }else if(colors == 3){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 9, name = "Greens")[5], len=0.09, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                    dev.off();
                }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 11, name = "PRGn")[2], len=0.09, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, 
                           col = brewer.pal(n = 9, name = "YlOrBr")[5])
                    dev.off();
                }
        }else{ #arrow not exists
                if(colors == 1){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "red")
                    dev.off();
                }else if(colors == 2){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 9, name = "YlOrBr")[5], len=0, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                    dev.off();
                }else if(colors == 3){
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 9, name = "Greens")[5], len=0, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                    dev.off();
                }else{
                    Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                    plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                       ar.col = brewer.pal(n = 11, name = "PRGn")[2], len=0, cex = sizes, to.target = T,
                       type = "points", xlab = "PC1", ylab = "PC2") 
                    points(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, 
                           col = brewer.pal(n = 9, name = "YlOrBr")[5])
                    dev.off();
                }
         }       
  }else{ # text
            if(arrow == 1){
              if(colors == 1){
                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                     ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
                     type = "text", xlab = "PC1", ylab = "PC2")
                  text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "red")
                  dev.off();
              }else if(colors == 2){
                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                     ar.col = brewer.pal(n = 9, name = "YlOrBr")[5], len=0.09, cex = sizes, to.target = T,
                     type = "text", xlab = "PC1", ylab = "PC2")
                  text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                  dev.off();
              }else if(colors == 3){
                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                     ar.col = brewer.pal(n = 9, name = "Greens")[5], len=0.09, cex = sizes, to.target = T,
                     type = "text", xlab = "PC1", ylab = "PC2")
                  text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                  dev.off();
              }else{
                  Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                  plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                     ar.col = brewer.pal(n = 11, name = "PRGn")[2], len=0.09, cex = sizes, to.target = T,
                     type = "text", xlab = "PC1", ylab = "PC2")
                  text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, 
                       col = brewer.pal(n = 9, name = "YlOrBr")[5])
                  dev.off();
              }
        }else{
            if(colors == 1){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                   ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0, cex = sizes, to.target = T,
                   type = "text", xlab = "PC1", ylab = "PC2")
                text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "red")
                dev.off();
            }else if(colors == 2){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                   ar.col = brewer.pal(n = 9, name = "YlOrBr")[5], len=0, cex = sizes, to.target = T,
                   type = "text", xlab = "PC1", ylab = "PC2")
                text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                dev.off();
            }else if(colors == 3){
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                   ar.col = brewer.pal(n = 9, name = "Greens")[5], len=0, cex = sizes, to.target = T,
                   type = "text", xlab = "PC1", ylab = "PC2")
                text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, col = "blue")
                dev.off();
            }else{
                Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
                plot(procrus, kind = 1, main = "Interomics Procrustes errors", 
                   ar.col = brewer.pal(n = 11, name = "PRGn")[2], len=0, cex = sizes, to.target = T,
                   type = "text", xlab = "PC1", ylab = "PC2")
                text(procrus, display = c("target", "rotated"), truemean = F, cex = sizes, 
                     col = brewer.pal(n = 9, name = "YlOrBr")[5])
                dev.off();
             }
           }
        }
  }


PlotPAplot2 <- function(imgName, dpi, format, colors){

  options(scipen = 999);
  dpi <- as.numeric(dpi)
  imgNm <- paste(imgName, ".", format, sep="");
  
  data1 <- readRDS("norm.data1");
  data2 <- readRDS("norm.data2");
  
  colors <- as.integer(colors)
  
  if(dim(data2)[2] < dim(data1)[2]){
    #standardize the data to means of 0 and standard deviations of 1
    data1 <- decostand(data1[,-1], method = "standardize")
    data2 <- decostand(data2[,-1], method = "standardize")
    
    #PCA based on Euclidean distance
    data1 <- rda(data1)
    data2 <- rda(data2)
    
    procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
  }else{
    #standardize the data to means of 0 and standard deviations of 1
    data1 <- decostand(data1[,-1], method = "standardize")
    data2 <- decostand(data2[,-1], method = "standardize")
    
    #PCA based on Euclidean distance
    data1 <- rda(data1)
    data2 <- rda(data2)
    
    procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
  }
  
      if(colors == 1){
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off();
      }else if(colors == 2){
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample", 
             col = brewer.pal(n = 9, name = "YlOrBr")[5], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off();
      }else if(colors == 3){
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample", 
             col = brewer.pal(n = 9, name = "Greens")[5], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off();
      }else{
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample", 
             col = brewer.pal(n = 11, name = "PRGn")[2], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off();
      }
} 
  

PerformQualityStatUpdate3<-function(perm){

    options(scipen = 999);
    perm <- as.integer(perm); # how many permutations
    
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    if(dim(data2)[2] < dim(data1)[2]){
      #standardize the data to means of 0 and standard deviations of 1
      data1 <- decostand(data1[,-1], method = "standardize")
      data2 <- decostand(data2[,-1], method = "standardize")
      
      #PCA based on Euclidean distance
      data1 <- rda(data1)
      data2 <- rda(data2)
      
      protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
    }else{
      #standardize the data to means of 0 and standard deviations of 1
      data1 <- decostand(data1[,-1], method = "standardize")
      data2 <- decostand(data2[,-1], method = "standardize")
      
      #PCA based on Euclidean distance
      data1 <- rda(data1)
      data2 <- rda(data2)
      
      protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
    }
    
    stat.info <- paste(
      "UPDATE!!! Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                        signif(protes$t0, 5), 
                       ". The mean and median of Procrustes correlations from permutations is ",
                       round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                        ", respectively, with newly selected permutations of ",
      protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
      round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 

    analSet$quality.stat.info.update3 <- stat.info;
    analSet<<-analSet;

    return(stat.info);
}



PerformQualityStatUpdate4<-function(){
    library(made4);
    options(scipen = 999);
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    xx1 <- prcomp(data1[,-1])
    xx2 <- prcomp(data2[,-1])
    
    coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])

    stat.info <- paste(
      "UPDATE!!! Coinertia RV coefficient is ", 
                      coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 

    analSet$quality.stat.info.update4 <- stat.info;
    analSet<<-analSet;

    return(stat.info);
}



################################
# PA and CA feature contributions

FeatureContributionPAplot1 <- function(imgName, dpi, format, choice, topn1, topn2, adjusting, joints, orth1, orth2, jointcomp){
    library(vegan)
    library(OmicsPLS)
    library(RColorBrewer)
    options(scipen = 999);
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
  
    
    choice <- as.integer(choice); #1. Kruskal 2. Anova 3. o2pls
    topn1 <- as.integer(topn1) #how many significant features to remove for omics1
    topn2 <- as.integer(topn2) #how many significant features to remove for omics2
    
    #<valid if choice == 1 or 2>
    adjusting <- as.integer(adjusting); #1. none 2. FDR -> only for choice == 1 or 2 
    
    #<valid if choice == 3>
    joints <- as.integer(joints) #joint component 
    orth1 <- as.integer(orth1) #omics1 orth
    orth2 <- as.integer(orth2) #omics2 orth
    jointcomp <- as.integer(jointcomp); #what joint score to use to find significant features
    
    
    
    if(choice == 1){
      kruskal <- function(x){
        kruskal.test(x ~ as.factor(data1[,1]))
      }
      kruskals <- apply(data1[,-1], 2, kruskal)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(kruskals)){
        p_val[i] <- kruskals[[i]]$p.value
      }
      
      kruskal2 <- function(x){
        kruskal.test(x ~ as.factor(data2[,1]))
      }
      kruskals2 <- apply(data2[,-1], 2, kruskal2)
      
      p_val2 <- c() 
      for(i in 1:length(kruskals2)){
        p_val2[i] <- kruskals2[[i]]$p.value
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 1, main = "Interomics Procrustes errors - Sig Features removed", 
             ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = 0.35, to.target = T,
             type = "text", xlab = "PC1", ylab = "PC2")
          text(procrus, display = c("target", "rotated"), truemean = F, cex = 0.35, col = "red")
          dev.off(); 
             
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 1, main = "Interomics Procrustes errors - Sig Features removed", 
             ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = 0.35, to.target = T,
             type = "text", xlab = "PC1", ylab = "PC2")
          text(procrus, display = c("target", "rotated"), truemean = F, cex = 0.35, col = "red")
          dev.off(); 
      }
    }else if(choice == 2){
      aof <- function(x){
        anova(aov(x ~ as.factor(data1[,1])))
      }
      anova.res <- apply(data1[,-1], 2, aof)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(anova.res)){
        p_val[i] <- anova.res[i][[1]][[5]][1]
      }
      
      aof2 <- function(x){
        anova(aov(x ~ as.factor(data2[,1])))
      }
      anova.res2 <- apply(data2[,-1], 2, aof2)
      
      p_val2 <- c() 
      for(i in 1:length(anova.res2)){
        p_val2[i] <- anova.res2[i][[1]][[5]][1]
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 1, main = "Interomics Procrustes errors - Sig Features removed", 
             ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = 0.35, to.target = T,
             type = "text", xlab = "PC1", ylab = "PC2")
          text(procrus, display = c("target", "rotated"), truemean = F, cex = 0.35, col = "red")
          dev.off(); 
             
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 1, main = "Interomics Procrustes errors - Sig Features removed", 
             ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = 0.35, to.target = T,
             type = "text", xlab = "PC1", ylab = "PC2")
          text(procrus, display = c("target", "rotated"), truemean = F, cex = 0.35, col = "red")
          dev.off(); 
      }
    }else{ #o2pls loading significant
      if(jointcomp > joints){
        return(1);
      }
      
      fit <- o2m(data1[,-1], data2[,-1], joints, orth1, orth2)
      load1x <- fit$W.[,jointcomp]
      load1y <- fit$C.[,jointcomp]

       if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
             data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 1, main = "Interomics Procrustes errors - Sig Features removed", 
             ar.col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = 0.35, to.target = T,
             type = "text", xlab = "PC1", ylab = "PC2")
          text(procrus, display = c("target", "rotated"), truemean = F, cex = 0.35, col = "red")
          dev.off();   
    }
}


FeatureContributionPAplot2 <- function(imgName, dpi, format, choice, topn1, topn2, adjusting, joints, orth1, orth2, jointcomp){
    options(scipen = 999);
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
  
    
    choice <- as.integer(choice); #1. Kruskal 2. Anova 3. o2pls
    topn1 <- as.integer(topn1) #how many significant features to remove for omics1
    topn2 <- as.integer(topn2) #how many significant features to remove for omics2
    
    #<valid if choice == 1 or 2>
    adjusting <- as.integer(adjusting); #1. none 2. FDR -> only for choice == 1 or 2 
    
    #<valid if choice == 3>
    joints <- as.integer(joints) #joint component 
    orth1 <- as.integer(orth1) #omics1 orth
    orth2 <- as.integer(orth2) #omics2 orth
    jointcomp <- as.integer(jointcomp); #what joint score to use to find significant features
    
    
    
    if(choice == 1){
      kruskal <- function(x){
        kruskal.test(x ~ as.factor(data1[,1]))
      }
      kruskals <- apply(data1[,-1], 2, kruskal)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(kruskals)){
        p_val[i] <- kruskals[[i]]$p.value
      }
      
      kruskal2 <- function(x){
        kruskal.test(x ~ as.factor(data2[,1]))
      }
      kruskals2 <- apply(data2[,-1], 2, kruskal2)
      
      p_val2 <- c() 
      for(i in 1:length(kruskals2)){
        p_val2[i] <- kruskals2[[i]]$p.value
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample - Sig features removed", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off(); 
             
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample - Sig features removed", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off(); 
      }
    }else if(choice == 2){
      aof <- function(x){
        anova(aov(x ~ as.factor(data1[,1])))
      }
      anova.res <- apply(data1[,-1], 2, aof)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(anova.res)){
        p_val[i] <- anova.res[i][[1]][[5]][1]
      }
      
      aof2 <- function(x){
        anova(aov(x ~ as.factor(data2[,1])))
      }
      anova.res2 <- apply(data2[,-1], 2, aof2)
      
      p_val2 <- c() 
      for(i in 1:length(anova.res2)){
        p_val2[i] <- anova.res2[i][[1]][[5]][1]
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample - Sig features removed", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off(); 
             
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample - Sig features removed", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off(); 
      }
    }else{ #o2pls loading significant
      if(jointcomp > joints){
        return(1);
      }
      
      fit <- o2m(data1[,-1], data2[,-1], joints, orth1, orth2)
      load1x <- fit$W.[,jointcomp]
      load1y <- fit$C.[,jointcomp]

       if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data1, data2, scale = T, symmetric = T)
          }else{
            #standardize the data to means of 0 and standard deviations of 1
             data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            procrus <- vegan::procrustes(data2, data1, scale = T, symmetric = T)
          }
          
          Cairo(file=imgNm, width=8, height=6, type=format, bg="white", dpi=dpi, unit="in");
          plot(procrus, kind = 2, main = "Impulse diagram of residuals for each sample - Sig features removed", 
             col = brewer.pal(n = 8, name = "PuBu")[7], len=0.09, cex = sizes, to.target = T,
             xlab = "", xaxt="n") 
          axis(1, at = 1:nrow(procrus$X),labels = rownames(procrus$X), las = 2, 
               tck = 0.01)
          mtext("Samples", side = 1, line = 4)
          dev.off();   
    }
}


FeatureContributionPAplot3 <- function(imgName, dpi, format, choice, topn1, topn2, adjusting, joints, orth1, orth2, jointcomp){
    library(gridExtra);
    library(reshape2);
    options(scipen = 999);
    dpi <- as.numeric(dpi)
    imgNm <- paste(imgName, ".", format, sep="");
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    fircolnm <- colnames(data1)[1];

    correct_colnames <- function(df) {
     delete.columns <- grep("(^X)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
    
      if (length(delete.columns) > 0) {
    
       colnames(df) <- gsub("^X", "",  colnames(df))
       #X might be replaced by different characters, instead of being deleted
      }
    
      return(df)
    }    
    
    data1 <- correct_colnames(data1);
    data2 <- correct_colnames(data2);
    data1 <- cbind(factor(data1[,1]), as.data.frame(data1[,-1]));
    data2 <- cbind(factor(data2[,1]), as.data.frame(data2[,-1]));
    colnames(data1)[1] <- fircolnm;
    colnames(data2)[1] <- fircolnm;
  
    
    choice <- as.integer(choice); #1. Kruskal 2. Anova 3. o2pls
    topn1 <- as.integer(topn1) #how many significant features to remove for omics1
    topn2 <- as.integer(topn2) #how many significant features to remove for omics2
    
    #<valid if choice == 1 or 2>
    adjusting <- as.integer(adjusting); #1. none 2. FDR -> only for choice == 1 or 2 
    
    #<valid if choice == 3>
    joints <- as.integer(joints) #joint component 
    orth1 <- as.integer(orth1) #omics1 orth
    orth2 <- as.integer(orth2) #omics2 orth
    jointcomp <- as.integer(jointcomp); #what joint score to use to find significant features
    
    
    ## data.frames with filled NA's when the number of rows are different 
    data.frame.na <- function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,
      stringsAsFactors = FALSE){
        data.row.names <- if (check.rows && is.null(row.names))
            function(current, new, i) {
                if (is.character(current))
                    new <- as.character(new)
                if (is.character(new))
                    current <- as.character(current)
                if (anyDuplicated(new))
                    return(current)
                if (is.null(current))
                    return(new)
                if (all(current == new) || all(current == ""))
                    return(new)
                stop(gettextf("mismatch of row names in arguments of 'data.frame', item %d",
                    i), domain = NA)
            }
        else function(current, new, i) {
            if (is.null(current)) {
                if (anyDuplicated(new)) {
                    warning("some row.names duplicated: ", paste(which(duplicated(new)),
                      collapse = ","), " --> row.names NOT used")
                    current
                }
                else new
            }
            else current
        }
        object <- as.list(substitute(list(...)))[-1L]
        mrn <- is.null(row.names)
        x <- list(...)
        n <- length(x)
        if (n < 1L) {
            if (!mrn) {
                if (is.object(row.names) || !is.integer(row.names))
                    row.names <- as.character(row.names)
                if (any(is.na(row.names)))
                    stop("row names contain missing values")
                if (anyDuplicated(row.names))
                    stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                      collapse = ", "))
            }
            else row.names <- integer(0L)
            return(structure(list(), names = character(0L), row.names = row.names,
                class = "data.frame"))
        }
        vnames <- names(x)
        if (length(vnames) != n)
            vnames <- character(n)
        no.vn <- !nzchar(vnames)
        vlist <- vnames <- as.list(vnames)
        nrows <- ncols <- integer(n)
        for (i in seq_len(n)) {
            xi <- if (is.character(x[[i]]) || is.list(x[[i]]))
                as.data.frame(x[[i]], optional = TRUE, stringsAsFactors = stringsAsFactors)
            else as.data.frame(x[[i]], optional = TRUE)
            nrows[i] <- .row_names_info(xi)
            ncols[i] <- length(xi)
            namesi <- names(xi)
            if (ncols[i] > 1L) {
                if (length(namesi) == 0L)
                    namesi <- seq_len(ncols[i])
                if (no.vn[i])
                    vnames[[i]] <- namesi
                else vnames[[i]] <- paste(vnames[[i]], namesi, sep = ".")
            }
            else {
                if (length(namesi))
                    vnames[[i]] <- namesi
                else if (no.vn[[i]]) {
                    tmpname <- deparse(object[[i]])[1L]
                    if (substr(tmpname, 1L, 2L) == "I(") {
                      ntmpn <- nchar(tmpname, "c")
                      if (substr(tmpname, ntmpn, ntmpn) == ")")
                        tmpname <- substr(tmpname, 3L, ntmpn - 1L)
                    }
                    vnames[[i]] <- tmpname
                }
            }
            if (missing(row.names) && nrows[i] > 0L) {
                rowsi <- attr(xi, "row.names")
                nc <- nchar(rowsi, allowNA = FALSE)
                nc <- nc[!is.na(nc)]
                if (length(nc) && any(nc))
                    row.names <- data.row.names(row.names, rowsi,
                      i)
            }
            nrows[i] <- abs(nrows[i])
            vlist[[i]] <- xi
        }
        nr <- max(nrows)
        for (i in seq_len(n)[nrows < nr]) {
            xi <- vlist[[i]]
            if (nrows[i] > 0L) {
                xi <- unclass(xi)
                fixed <- TRUE
                for (j in seq_along(xi)) {
                    ### added NA fill to max length/nrow
                    xi1 <- xi[[j]]
                    if (is.vector(xi1) || is.factor(xi1))
                      xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                    else if (is.character(xi1) && class(xi1) == "AsIs")
                      xi[[j]] <- structure(c(xi1, rep(NA, nr - nrows[i])),
                        class = class(xi1))
                    else if (inherits(xi1, "Date") || inherits(xi1,
                      "POSIXct"))
                      xi[[j]] <- c(xi1, rep(NA, nr - nrows[i]))
                    else {
                      fixed <- FALSE
                      break
                    }
                }
                if (fixed) {
                    vlist[[i]] <- xi
                    next
                }
            }
            stop("arguments imply differing number of rows: ", paste(unique(nrows),
                collapse = ", "))
        }
        value <- unlist(vlist, recursive = FALSE, use.names = FALSE)
        vnames <- unlist(vnames[ncols > 0L])
        noname <- !nzchar(vnames)
        if (any(noname))
            vnames[noname] <- paste("Var", seq_along(vnames), sep = ".")[noname]
        if (check.names)
            vnames <- make.names(vnames, unique = TRUE)
        names(value) <- vnames
        if (!mrn) {
            if (length(row.names) == 1L && nr != 1L) {
                if (is.character(row.names))
                    row.names <- match(row.names, vnames, 0L)
                if (length(row.names) != 1L || row.names < 1L ||
                    row.names > length(vnames))
                    stop("row.names should specify one of the variables")
                i <- row.names
                row.names <- value[[i]]
                value <- value[-i]
            }
            else if (!is.null(row.names) && length(row.names) !=
                nr)
                stop("row names supplied are of the wrong length")
        }
        else if (!is.null(row.names) && length(row.names) != nr) {
            warning("row names were found from a short variable and have been discarded")
            row.names <- NULL
        }
        if (is.null(row.names))
            row.names <- .set_row_names(nr)
        else {
            if (is.object(row.names) || !is.integer(row.names))
                row.names <- as.character(row.names)
            if (any(is.na(row.names)))
                stop("row names contain missing values")
            if (anyDuplicated(row.names))
                stop("duplicate row.names: ", paste(unique(row.names[duplicated(row.names)]),
                    collapse = ", "))
        }
        attr(value, "row.names") <- row.names
        attr(value, "class") <- "data.frame"
        value
    }
    
    if(choice == 1){
      kruskal <- function(x){
        kruskal.test(x ~ as.factor(data1[,1]))
      }
      kruskals <- apply(data1[,-1], 2, kruskal)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(kruskals)){
        p_val[i] <- kruskals[[i]]$p.value
      }
      
      kruskal2 <- function(x){
        kruskal.test(x ~ as.factor(data2[,1]))
      }
      kruskals2 <- apply(data2[,-1], 2, kruskal2)
      
      p_val2 <- c() 
      for(i in 1:length(kruskals2)){
        p_val2[i] <- kruskals2[[i]]$p.value
      }
      
      if(adjusting == 1){
          databrick <- data.frame.na(
              feature_omics1_removed = colnames(data1[,(order(p_val, decreasing = F) + 1)[1:topn1]]),
              feature_omics1_p.val = round(sort(p_val, decreasing = F)[1:topn1], 5),
              feature_omics2_removed = colnames(data2[,(order(p_val2, decreasing = F) + 1)[1:topn2]]),
              feature_omics2_p.val = round(sort(p_val2, decreasing = F)[1:topn2], 5)
          ) 
          
          Cairo(file=imgNm, width=13, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
          grid.table(databrick);
          dev.off(); 
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         databrick <- data.frame.na(
              feature_omics1_removed = colnames(data1[,(order(fdr, decreasing = F) + 1)[1:topn1]]),
              feature_omics1_FDRp.val = round(sort(fdr, decreasing = F)[1:topn1], 5),
              feature_omics2_removed = colnames(data2[,(order(fdr2, decreasing = F) + 1)[1:topn2]]),
              feature_omics2_FDRp.val = round(sort(fdr2, decreasing = F)[1:topn2], 5)
          ) 
          
          Cairo(file=imgNm, width=13, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
          grid.table(databrick);
          dev.off(); 
      }
    }else if(choice == 2){
      aof <- function(x){
        anova(aov(x ~ as.factor(data1[,1])))
      }
      anova.res <- apply(data1[,-1], 2, aof)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(anova.res)){
        p_val[i] <- anova.res[i][[1]][[5]][1]
      }
      
      aof2 <- function(x){
        anova(aov(x ~ as.factor(data2[,1])))
      }
      anova.res2 <- apply(data2[,-1], 2, aof2)
      
      p_val2 <- c() 
      for(i in 1:length(anova.res2)){
        p_val2[i] <- anova.res2[i][[1]][[5]][1]
      }
      
      if(adjusting == 1){
          databrick <- data.frame.na(
              feature_omics1_removed = colnames(data1[,(order(p_val, decreasing = F) + 1)[1:topn1]]),
              feature_omics1_p.val = round(sort(p_val, decreasing = F)[1:topn1], 5),
              feature_omics2_removed = colnames(data2[,(order(p_val2, decreasing = F) + 1)[1:topn2]]),
              feature_omics2_p.val = round(sort(p_val2, decreasing = F)[1:topn2], 5)
          ) 
          
          Cairo(file=imgNm, width=11, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
          grid.table(databrick);
          dev.off(); 
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         databrick <- data.frame.na(
              feature_omics1_removed = colnames(data1[,(order(fdr, decreasing = F) + 1)[1:topn1]]),
              feature_omics1_FDRp.val = round(sort(fdr, decreasing = F)[1:topn1], 5),
              feature_omics2_removed = colnames(data2[,(order(fdr2, decreasing = F) + 1)[1:topn2]]),
              feature_omics2_FDRp.val = round(sort(fdr2, decreasing = F)[1:topn2], 5)
          ) 
          
          Cairo(file=imgNm, width=13, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
          grid.table(databrick);
          dev.off(); 
      }
    }else{ #o2pls loading significant
        if(jointcomp > joints){
          return(1);
        }
        
        fit <- o2m(data1[,-1], data2[,-1], joints, orth1, orth2)
        load1x <- fit$W.[,jointcomp]
        load1y <- fit$C.[,jointcomp]
        load1xjoint <- melt(sort(abs(load1x), decreasing = T)[1:topn1]) #sum of joint
        #colnames(data1[,(order(abs(load1x), decreasing = T)[1:topn1]+1)])
        load1yjoint <- melt(sort(abs(load1y), decreasing = T)[1:topn2])

        databrick <- data.frame.na(
              feature_omics1_removed = rownames(load1xjoint),
              feature_omics1_loading = load1xjoint$value,
              feature_omics2_removed = rownames(load1yjoint),
              feature_omics2_loading = load1yjoint$value
          ) 
          
        Cairo(file=imgNm, width=13, height=nrow(databrick)/2.6, type=format, bg="white", dpi=dpi, unit="in");
        grid.table(databrick);
        dev.off(); 
    }
}



PerformQualityStatUpdate5 <- function(choice, topn1, topn2, adjusting, joints, orth1, orth2, jointcomp, perm){
    options(scipen = 999);
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    choice <- as.integer(choice); #1. Kruskal 2. Anova 3. o2pls
    topn1 <- as.integer(topn1) #how many significant features to remove for omics1
    topn2 <- as.integer(topn2) #how many significant features to remove for omics2
    
    #<valid if choice == 1 or 2>
    adjusting <- as.integer(adjusting); #1. none 2. FDR -> only for choice == 1 or 2 
    
    #<valid if choice == 3>
    joints <- as.integer(joints) #joint component 
    orth1 <- as.integer(orth1) #omics1 orth
    orth2 <- as.integer(orth2) #omics2 orth
    jointcomp <- as.integer(jointcomp); #what joint score to use to find significant features
    
    perm <- as.integer(perm); # how many permutations
    
    if(choice == 1){
      kruskal <- function(x){
        kruskal.test(x ~ as.factor(data1[,1]))
      }
      kruskals <- apply(data1[,-1], 2, kruskal)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(kruskals)){
        p_val[i] <- kruskals[[i]]$p.value
      }
      
      kruskal2 <- function(x){
        kruskal.test(x ~ as.factor(data2[,1]))
      }
      kruskals2 <- apply(data2[,-1], 2, kruskal2)
      
      p_val2 <- c() 
      for(i in 1:length(kruskals2)){
        p_val2[i] <- kruskals2[[i]]$p.value
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
          }
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with Kruskal-Wallis with p-value, Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                              signif(protes$t0, 5), 
                             ". The mean and median of Procrustes correlations from permutations is ",
                             round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                              ", respectively, with newly selected permutations of ",
            protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
            round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 
      
          analSet$quality.stat.info.update5 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
          }
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with Kruskal-Wallis with FDR adjusted p-value, Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                              signif(protes$t0, 5), 
                             ". The mean and median of Procrustes correlations from permutations is ",
                             round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                              ", respectively, with newly selected permutations of ",
            protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
            round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 
      
          analSet$quality.stat.info.update5 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }
    }else if(choice == 2){
      aof <- function(x){
        anova(aov(x ~ as.factor(data1[,1])))
      }
      anova.res <- apply(data1[,-1], 2, aof)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(anova.res)){
        p_val[i] <- anova.res[i][[1]][[5]][1]
      }
      
      aof2 <- function(x){
        anova(aov(x ~ as.factor(data2[,1])))
      }
      anova.res2 <- apply(data2[,-1], 2, aof2)
      
      p_val2 <- c() 
      for(i in 1:length(anova.res2)){
        p_val2[i] <- anova.res2[i][[1]][[5]][1]
      }
      
      if(adjusting == 1){
          if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
          }
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with ANOVA with p-value, Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                              signif(protes$t0, 5), 
                             ". The mean and median of Procrustes correlations from permutations is ",
                             round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                              ", respectively, with newly selected permutations of ",
            protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
            round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 
      
          analSet$quality.stat.info.update5 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
          }else{
            #standardize the data to means of 0 and standard deviations of 1
            data1 <- decostand(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
          }
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with ANOVA with FDR adjusted p-value, Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                              signif(protes$t0, 5), 
                             ". The mean and median of Procrustes correlations from permutations is ",
                             round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                              ", respectively, with newly selected permutations of ",
            protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
            round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 
      
          analSet$quality.stat.info.update5 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }
    }else{ #o2pls loading significant
      if(jointcomp > joints){
        return(1);
      }
      
      fit <- o2m(data1[,-1], data2[,-1], joints, orth1, orth2)
      load1x <- fit$W.[,jointcomp]
      load1y <- fit$C.[,jointcomp]

       if((dim(data2)[2] - topn2) < (dim(data1)[2] - topn1)){
            #standardize the data to means of 0 and standard deviations of 1
            # +1 as the first column is a sample metadata
            data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data1, data2, scale = T, symmetric = T, permutations = how(nperm = perm))
          }else{
            #standardize the data to means of 0 and standard deviations of 1
             data1 <- decostand(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))], method = "standardize")
            data2 <- decostand(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))], method = "standardize")
            
            #PCA based on Euclidean distance
            data1 <- rda(data1)
            data2 <- rda(data2)
            
            protes <- vegan::protest(data2, data1, scale = T, symmetric = T, permutations = how(nperm = perm))
          }
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with O2PLS absolute joint loading values, Procrustes correlation from non-permuted value in a symmetric Procrustes rotation is ", 
                              signif(protes$t0, 5), 
                             ". The mean and median of Procrustes correlations from permutations is ",
                             round(mean(protes$t), 5), " and ", round(unname(quantile(protes$t)[3]), 5), 
                              ", respectively, with newly selected permutations of ",
            protes$permutations, " times. ", "The significance of Procrustes correlations is P-val <= ",
            round(protes$signif, 5), ". (the minimum significance is ", round((1/(1+protes$permutations)), 5), ")", sep=""); 
      
          analSet$quality.stat.info.update5 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);  
    }
}


PerformQualityStatUpdate6 <- function(choice, topn1, topn2, adjusting, joints, orth1, orth2, jointcomp){
    library(made4);  
    options(scipen = 999);
    
    data1 <- readRDS("norm.data1");
    data2 <- readRDS("norm.data2");
    
    choice <- as.integer(choice); #1. Kruskal 2. Anova 3. o2pls
    topn1 <- as.integer(topn1) #how many significant features to remove for omics1
    topn2 <- as.integer(topn2) #how many significant features to remove for omics2
    
    #<valid if choice == 1 or 2>
    adjusting <- as.integer(adjusting); #1. none 2. FDR -> only for choice == 1 or 2 
    
    #<valid if choice == 3>
    joints <- as.integer(joints) #joint component 
    orth1 <- as.integer(orth1) #omics1 orth
    orth2 <- as.integer(orth2) #omics2 orth
    jointcomp <- as.integer(jointcomp); #what joint score to use to find significant features
    
    if(choice == 1){
      kruskal <- function(x){
        kruskal.test(x ~ as.factor(data1[,1]))
      }
      kruskals <- apply(data1[,-1], 2, kruskal)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(kruskals)){
        p_val[i] <- kruskals[[i]]$p.value
      }
      
      kruskal2 <- function(x){
        kruskal.test(x ~ as.factor(data2[,1]))
      }
      kruskals2 <- apply(data2[,-1], 2, kruskal2)
      
      p_val2 <- c() 
      for(i in 1:length(kruskals2)){
        p_val2[i] <- kruskals2[[i]]$p.value
      }
      
      if(adjusting == 1){
          xx1 <- prcomp(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])])
          xx2 <- prcomp(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])])
        
          coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with Kruskal-Wallis with p-value, Coinertia RV coefficient is ", 
                            coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 
      
          analSet$quality.stat.info.update6 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         xx1 <- prcomp(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])])
         xx2 <- prcomp(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])])

         coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with Kruskal-Wallis with FDR adjusted p-value, Coinertia RV coefficient is ", 
                            coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 
      
          analSet$quality.stat.info.update6 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }
    }else if(choice == 2){
      aof <- function(x){
        anova(aov(x ~ as.factor(data1[,1])))
      }
      anova.res <- apply(data1[,-1], 2, aof)
      
      p_val <- c() #p.val small --> corresponding feature is different along the sample meta-data
      for(i in 1:length(anova.res)){
        p_val[i] <- anova.res[i][[1]][[5]][1]
      }
      
      aof2 <- function(x){
        anova(aov(x ~ as.factor(data2[,1])))
      }
      anova.res2 <- apply(data2[,-1], 2, aof2)
      
      p_val2 <- c() 
      for(i in 1:length(anova.res2)){
        p_val2[i] <- anova.res2[i][[1]][[5]][1]
      }
      
      if(adjusting == 1){
          xx1 <- prcomp(data1[,-c(1, (order(p_val, decreasing = F) + 1)[1:topn1])])
          xx2 <- prcomp(data2[,-c(1, (order(p_val2, decreasing = F) + 1)[1:topn2])])
        
          coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with ANOVA with p-value, Coinertia RV coefficient is ", 
                            coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 
      
          analSet$quality.stat.info.update6 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }else{
         fdr <- p.adjust(p_val, method = "BH")
         fdr2 <- p.adjust(p_val2, method = "BH")
         
         xx1 <- prcomp(data1[,-c(1, (order(fdr, decreasing = F) + 1)[1:topn1])])
         xx2 <- prcomp(data2[,-c(1, (order(fdr2, decreasing = F) + 1)[1:topn2])])

         coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])
          
          stat.info <- paste(
            "UPDATE!!! After removing significant features with ANOVA with FDR adjusted p-value, Coinertia RV coefficient is ", 
                            coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 
      
          analSet$quality.stat.info.update6 <- stat.info;
          analSet<<-analSet;
      
          return(stat.info);
      }
    }else{ #o2pls loading significant
      if(jointcomp > joints){
        return(1);
      }
      
      fit <- o2m(data1[,-1], data2[,-1], joints, orth1, orth2)
      load1x <- fit$W.[,jointcomp]
      load1y <- fit$C.[,jointcomp]

      xx1 <- prcomp(data1[,-c(1, (order(abs(load1x), decreasing = T)[1:topn1]+1))])
      xx2 <- prcomp(data2[,-c(1, (order(abs(load1y), decreasing = T)[1:topn2]+1))])
         
      coinertias <- cia(xx1$x[,1:5], xx2$x[,1:5])
          
      stat.info <- paste(
        "UPDATE!!! After removing significant features with O2PLS absolute joint loading values, Coinertia RV coefficient is ", 
                        coinertias$coinertia$RV, ", based on five PCA components.", sep=""); 
  
      analSet$quality.stat.info.update6 <- stat.info;
      analSet<<-analSet;
  
      return(stat.info);
    }
}
