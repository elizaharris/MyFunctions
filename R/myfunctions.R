
#' Function to convert time from a string to a Posix object
#' 
#' @param timechar A string or array of strings containing the time/date
#' @param t_form The format of the time (see strptime for format details) - default t_form="%Y-%m-%d %H:%M:%S"
#' @examples time.c("2019-01-01 01:01:01", t_form="%Y-%m-%d %H:%M:%S")
#' @export
time.c <- function(timechar,t_form="%Y-%m-%d %H:%M:%S"){
  as.POSIXct(strptime(timechar,t_form), tz = "CEST") } 

#' Function to fill in times at a certain freq
#' 
#' Useful to create a continuous time vec for interpolation
#' @param startt Start datetime (Posix)
#' @param endt End datetime (Posix)
#' @param freq Spacing of times to fill in (must be in seconds!)
#' @param updown If updown == "up" times are rounded up to the nearest second; otherwise down
#' @export
filltime = function(startt,endt,freq,updown){ 
  if (missing(updown)){updown="up"} # default values
  if (missing(freq)){freq=1}
  a = as.numeric(difftime(endt,startt,units="secs"))/freq
  if (updown == "up"){a = ceiling(a)} 
  if (updown == "down"){a = floor(a)}
  x = rep(startt,a)
  for (n in 1:a){ x[n] = startt + (n-1)*freq }
  return(x) }

#' Function for weighted mean and standard deviation
#' 
#' @param x Vector of data
#' @param w Vector of weights or standard deviations
#' @param weightsaresd If TRUE (default) then the weights are assumed to be the standard deviation
#' in the data, eg. a high value means a low weighting. If weightsaresd = FALSE the weights are
#' used directly, eg. a high weight is a high weighting. If weightsaresd = FALSE the returned 
#' standard deviation in the weighted mean may be incorrect.
#' @return A vector of three numbers
#' \itemize{
#'   \item The weighted mean
#'   \item The weighted standard deviation
#'   \item The number of non-nan value-weight pairs
#' }
#' @export
weighted.m.sd = function(x,w,weightsaresd = TRUE){
  if (length(x)!=length(w)){stop("Weights and data are not the same length")}
  tmp = which(!is.na(x) & !is.na(w) & is.finite(x) & is.finite(w))
  x = x[tmp]; w = w[tmp]
  if (length(x)>1){ 
    if (weightsaresd){
      b = weighted.mean(x,1/w)
    } else {b = weighted.mean(x,w)}
    top = sum(w*((x-b)^2),na.rm=TRUE) # numerator of weighted sd
    bottom = (length(w)-1)*sum(w)/length(w)
    wsd = sqrt(top/bottom)} 
  if (length(x)==1){ print("Only one non-NaN point; returning input data")
    b=x; wsd=w }
  if (length(x)==0){ print("No non-NaN points; returning NaNs")
    b=NaN; wsd=NaN }
  return(c(b,wsd,length(x)))
}

# find flattest portion of peak
flatpeak = function(data,sdthresh){
  # first remove big conc changes
  diffs = abs(data[11:length(data)]-data[1:(length(data)-10)])
  if (length(unique(diffs))>5){ # check for measurement artefact; if instrument briefly crashed data increases/decreases monotonically
    if (2*sd(diffs)>1){diffthresh = 2*sd(diffs)} else {diffthresh=1}
    b = c(which(diffs>(mean(diffs)+diffthresh)),length(data))
    if (length(b)>1){b = b[which(c(1,diff(b))>2)]}
    if (b[1]!=1){b=c(1,b)}
    sdb = b[1:(length(b)-1)]
    for (n in 1:(length(b)-1)){ # sd between each "large" change normalised to length
      sdb[n] = sd(data[b[n]:b[n+1]])/sqrt(b[n+1]-b[n])} 
    b = c((b[which.min(sdb)]),(b[which.min(sdb)+1]))
    if (max(b)>length(data)){b[2]=length(data)}
    if (length(b)==0 | diff(b)==0){ b = c(1,length(data)) }
    # then filter for sd threshold
    istart = b[1]
    iend = b[2]
    finished = 0
    while (finished!=2){
      low = mean(data[istart:iend])-sdthresh*sd(data[istart:iend])
      high = mean(data[istart:iend])+sdthresh*sd(data[istart:iend])
      if (data[istart]<low | data[istart]>high){istart = istart+5} else {finished=finished+1}
      if (data[iend]<low | data[iend]>high){iend = iend-5} else {finished=finished+1}
      if (finished!=2){finished=0}
      if (iend < istart+10){finished=2} 
    }
    } else { istart=0; iend=0 }
  return(c(istart,iend))
}

# linear comparison with plot
par(mfrow=c(2,2))
compare = function(x,y,np=TRUE,cols="black",wts=NULL,xlimits=c(min(x,na.rm=TRUE),max(x,na.rm=TRUE)),ylimits=c(min(y,na.rm=TRUE),max(y,na.rm=TRUE)),ylabel="y"){
  if (np==TRUE){ plot(x,y,col=cols,xlim=xlimits,ylim=ylimits) } else { points(x,y,col=cols) }
  res = lm(y~x,weights=wts)
  print(summary(res)$coefficients[2,4])
  if (summary(res)$coefficients[2,4] < 0.1){
    lines(x,x*summary(res)$coefficients[2,1]+
            summary(res)$coefficients[1,1],col="red") }
  return(summary(res)$coefficients[2,4])
}

# wilcox test that checks if both sets of observations have enough non-NA values and returns NA if not, otherwise p value
wilcox.check = function(x,y,paired=FALSE){ 
  x = as.numeric(x)
  y = as.numeric(y)
  if (sum(!is.na(x))>1 & sum(!is.na(y))>1){ res = wilcox.test(x,y,paired=paired)$p.value } else {res = NA}
  return(res) }

# t-test that checks if both sets of observations have enough non-NA values and returns NA if not, otherwise p value
ttest.check = function(x,y,paired=FALSE){ 
  x = as.numeric(x)
  y = as.numeric(y)
  if (sum(!is.na(x))>1 & sum(!is.na(y))>1){ res = t.test(x,y,paired=paired)$p.value } else {res = NA}
  return(res) }

# RMSE
RMSE = function(x,y){
  diffs = x-y
  diffs = diffs[!is.na(diffs)]
  res = (sum(diffs^2)/length(diffs))^0.5
  return(res)
}

# plot PCA results
pca.biplot = function(res.pca,components=2,obs=TRUE){ # PCA results, number of components to plot must be given (component 1 always plotted!)
  PCs = summary(res.pca)$importance # summary of the components
  var.wts = res.pca$rotation # weightings of the variables
  obs.wts = res.pca$x # weightings of the observations
  par(mfrow=c(2,1),mar=c(4,4,1,1)) # set margins
  for (n in 1:components){ # plot as many components as are selected
    plot(var.wts[,1],var.wts[,n+1],col="cyan",pch=16,xlab=paste0("PC1 (",round(PCs[2,1]*100),"%)"),ylab=paste0("PC",n+1," (",round(PCs[2,n+1]*100),"%)"))
    text(var.wts[,1],var.wts[,n+1],col="blue",labels=rownames(var.wts),cex=0.5,pos=1)
    lines(c(-1,1),c(0,0),col="grey") # axis lines
    lines(c(0,0),c(-1,1),col="grey") # axis lines
    for (i in 1:length(var.wts[,1])){ # a line to each variable on the PCA plot
      lines(c(0,var.wts[i,1]),c(0,var.wts[i,n+1]),col="cyan") }
    if (obs){
      plot(obs.wts[,1],obs.wts[,n+1],col="red",pch=3,xlab=paste0("PC1 (",round(PCs[2,1]*100),"%)"),ylab=paste0("PC",n+1," (",round(PCs[2,n+1]*100),"%)"))
      text(obs.wts[,1],obs.wts[,n+1],labels=rownames(obs.wts),col="red",cex=0.5,pos=1) }
  }
}

# combine files
combineFiles = function(path,filepattern,fileformat,sheetchoice="Daten",save="n"){
    filenames = list.files(path)
    filenames = filenames[grep(filepattern,filenames)]
    if (fileformat=="xlsx"){
        for (n in seq_along(filenames)){
            if (n==1){
                names = colnames(read_excel(paste0(path,filenames[n]),skip=1,sheet=sheetchoice))
                data = read_excel(paste0(path,filenames[n]),col_names=names,skip=4,sheet=sheetchoice)
            } else {
                newnames = colnames(read_excel(paste0(path,filenames[n]),skip=1,sheet=sheetchoice))
                newdata = read_excel(paste0(path,filenames[n]),col_names=names,skip=4,sheet=sheetchoice)
                data = bind_rows(data,newdata)
            }
        } }
    if (fileformat=="dat"){
        for (n in seq_along(filenames)){
            print(n)
            if (n==1){
                names = colnames(read_delim(paste0(path,filenames[n]),skip=1,delim=","))
                data = read_delim(paste0(path,filenames[n]),col_names=names,skip=4,delim=",")
            } else {
                names = colnames(read_delim(paste0(path,filenames[n]),skip=1,delim=","))
                newdata = read_delim(paste0(path,filenames[n]),col_names=names,skip=4,delim=",")
                data = bind_rows(data,newdata)
            }
        } }
    if (save!="n"){
        write.csv(data,paste0(path,save,".csv"))
    }
    return(data)
}


