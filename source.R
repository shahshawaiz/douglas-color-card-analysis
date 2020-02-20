## ----echo=FALSE, message=FALSE, warning=FALSE, results=FALSE-------------
# 
  setwd(".")

# 
  library(raster)
  library(colorspace)
  library(scatterplot3d)
  library(plyr)
  library(rgl)
  library(colordistance)
  library(ggplot2)
  library(grid)

  options(rgl.useNULL = TRUE)
  options(rgl.printRglwidget = TRUE)
  
# 
measrmnts <- read.csv2("LabMeasurements-Color-Card.csv")
master <- read.csv2("MasterColorCard.csv")

# add sheet column
  measrmnts['Sheet'] <- rep(1:13, 42)



## ----echo=FALSE, message=FALSE, warning=FALSE, results=FALSE-------------

# calculate delta E
getColorDifferenceDf <- function(df1, df2){
  
  # res <- sqrt(
  #   abs(spot1[1]-spot2[1])^2 - abs(spot1[2]-spot2[2])^2 - abs(spot1[3]-spot2[3])^2
  # )

  res <- df1 - df2
  
  return(res)
}



## ----echo=FALSE, message=FALSE, warning=FALSE, results=FALSE-------------

# 
groupBySheet <- function(dat){
  groupColumns = c('sheet')
  dataColumns = names(dat)[1:3]
  measrmnts_grouped_by_sheet <- ddply(
    dat, groupColumns, function(x) colMeans(x[dataColumns])
  )
  return(measrmnts_grouped_by_sheet)
}

# 
groupByColorSpot <- function(dat){
  groupColumns = c('spotRow', 'spotCol')
  dataColumns = names(dat)[1:3]
  measrmnts_grouped_by_colorspot <- ddply(
    dat, groupColumns, function(x) colMeans(x[dataColumns])
  )
  return(measrmnts_grouped_by_colorspot)
}

groupByTarget <- function(dat){
  groupColumns = c('targetRow', 'targetCol')
  dataColumns = names(dat)[1:3]
  measrmnts_grouped_by_colorspot <- ddply(
    dat, groupColumns, function(x) colMeans(x[dataColumns])
  )
  return(measrmnts_grouped_by_colorspot)
}

# prase master card
parseMaster <- function(df){

  res <- df[
    c('L','a','b')
  ]

  return(res)
}

# parse measurement cards
parseMeasurment <- function(df){

  res <- data.frame()

  # iterate measurments
  for(i in 1:8) {

    for(j in 1:8) {
      cardSlice <- df[
        c(
          paste0("L",i,j),
          paste0("a",i,j),
          paste0("b",i,j)
        )
      ]
      
      names(cardSlice) <- c("L","a", "b")
      
      cardSlice$spotRow <- i
      cardSlice$spotCol <- j
      cardSlice$targetRow <- df$Row
      cardSlice$targetCol <- df$Column
      cardSlice$sheet <- df$Sheet
      
      res <- rbind(
        res, data.frame(cardSlice)
      )

    }
  }
  return(res)
}

# parse measurement cards
parseMeasurmentBySheetNo <- function(dfSheet, sheetNo){
  return(dfSheet[sheetNo,])
}

datMaster <- parseMaster(master)
datMeasurements <- parseMeasurment(measrmnts)



## ----echo=FALSE, message=FALSE, warning=FALSE, results=FALSE-------------

# group data by sheet, target, colorspot
measrmnts_grouped_by_sheet <- groupBySheet(datMeasurements)
measrmnts_grouped_by_colorspot <- groupByColorSpot(datMeasurements)
measrmnts_grouped_by_target <- groupByTarget(datMeasurements)



## ----echo=FALSE, message=FALSE, warning=FALSE, results=FALSE-------------

library(hexbin)

LabMeasurements_color_code <- read.csv2("LabMeasurements-Color-Card.csv", header = TRUE, sep = ";", dec = ",")
MasterColorCard <- read.csv2("MasterColorCard.csv", header = TRUE, sep = ";", dec = ",")

### calculate delta E
getColorDifference <- function(spot1, spot2, command){
  res<-0
  if (command=="all"){
    res <- sqrt( 
      (spot1$L-spot2$L)^2 + (spot1$a-spot2$a)^2 + (spot1$b-spot2$b)^2
    ) 
  }
  else if(command=="L"){
    res <- sqrt((spot1$L-spot2$L)^2)
  }
  else if(command=="a"){
    res <- sqrt((spot1$a-spot2$a)^2)
  }
  else{
    res <- sqrt((spot1$b-spot2$b)^2)
  }
  
  return(res)
}


pos_avg13_sheets<-function(pos_lab){
  pos_color_sheet <- data.frame(Crow=double(),
                                Ccol=double(),
                                L=double(),
                                a=double(),
                                b=double(),
                                stringsAsFactors=FALSE)
  Ccol<-1
  Crow<-1
  lower<-3
  for (i in 1:64){
    
    upper<-lower+2
    temp<-pos_lab[,lower:upper]
    names(temp)<-c("L","a","b")
    temp["Crow"] <- Crow
    temp["Ccol"] <- Ccol
    pos_color_sheet <- rbind(pos_color_sheet, temp)
    lower<-lower+3
    if (Ccol%%8==0){
      Crow<-Crow+1 
      Ccol<-1
    }
    else{
      Ccol<-Ccol+1 
    }
  }
  return (pos_color_sheet)
}

avg_lab_13sheets<-function(LabMeasurements_color_code){
  avg_color_sheet <- data.frame(Crow=double(),
                                Ccol=double(),
                                L=double(),
                                a=double(),
                                b=double(),
                                stringsAsFactors=FALSE)
  
  avg_lab<-LabMeasurements_color_code[,3:194]
  avg_lab<-as.data.frame(t(apply((avg_lab),MARGIN=2,FUN=mean)))
  Ccol<-1
  Crow<-1
  lower<-1
  for (i in 1:64){
    upper<-lower+2
    #print("NEW")
    #print (i)
    #print(lower)
    #print(upper)
    #print(avg_lab[lower:upper])
    temp<-avg_lab[lower:upper]
    names(temp)<-c("L","a","b")
    temp["Crow"] <- Crow
    temp["Ccol"] <- Ccol
    avg_color_sheet <- rbind(avg_color_sheet, temp)
    lower<-lower+3
    if (Ccol%%8==0){
      Crow<-Crow+1 
      Ccol<-1
    }
    else{
      Ccol<-Ccol+1 
    }
  }
  return (avg_color_sheet)
}

# Extract average L,a,b values of 13 color sheets
avg_color_sheet<-avg_lab_13sheets(LabMeasurements_color_code)

get_lab_byPosition<-function(LabMeasurements_color_code){

  pos<-matrix(list(), nrow=7, ncol=6)
  delta_e_pos<-matrix(, nrow=7, ncol=6)
  delta_L_pos<-matrix(, nrow=7, ncol=6)
  delta_a_pos<-matrix(, nrow=7, ncol=6)
  delta_b_pos<-matrix(, nrow=7, ncol=6)
  for (i in 1:7){
    for (j in 1:6){
      #print(c(i,j))
      temp_store<-LabMeasurements_color_code[which(LabMeasurements_color_code$Row == i & LabMeasurements_color_code$Column == j), ]
      #pos[[i,j]]<-as.data.frame(t(apply((temp_store),MARGIN=2,FUN=mean)))
      temp<-as.data.frame(t(apply((temp_store),MARGIN=2,FUN=mean)))
      #print(temp)
      #print (c(temp$Row,temp$Column))
      pos_sheets<-pos_avg13_sheets(temp)
      delta_l<-mean(getColorDifference(MasterColorCard,pos_sheets,"L"))
      delta_a<-mean(getColorDifference(MasterColorCard,pos_sheets,"a"))
      delta_b<-mean(getColorDifference(MasterColorCard,pos_sheets,"b"))
      delta_e<-mean(getColorDifference(MasterColorCard,pos_sheets,"all"))

      delta_e_pos[i,j]<-delta_e
      delta_L_pos[i,j]<-delta_l
      delta_a_pos[i,j]<-delta_a
      delta_b_pos[i,j]<-delta_b
    }
  }
  return(new("ByIndex",delta_l=delta_L_pos,delta_a=delta_a_pos,delta_b=delta_b_pos,delta_e=delta_e_pos))
}

#Extract values with respects to the color sheets position

setClass(Class="ByIndex",
           representation(
             delta_l="matrix",
             delta_a="matrix",
             delta_b="matrix",
             delta_e="matrix"
           )
  )
pos<-get_lab_byPosition(LabMeasurements_color_code)



## ----echo=FALSE, message=FALSE, warning=FALSE, fig1, fig.align="center", fig.height=8, out.width = '100%'----

# dimension comparision
  deltaE <- getColorDifferenceDf(datMaster, measrmnts_grouped_by_colorspot[-1:-2])

  par(mfrow=c(3,1))

  # L 
    plot(
      datMaster$L, type = "o", col = 10, pch=17, xlab="Colorspot", ylab = "Color Level", main = "Channel Comparision, Master vs Measurement: L")
    lines(measrmnts_grouped_by_colorspot$L, type = "o", col = 100, pch=19)
    legend("bottomleft", 
      legend = c("Master", "Measurement"), 
      col = c(10,100),
      pch = c(17,19)
    )

  # A
  plot(
    datMaster$a, type = "o", col = 10, pch=17, xlab="Colorspot", ylab = "Color Level", main = "Channel Comparision, Master vs Measurement: a")
  lines(measrmnts_grouped_by_colorspot$a, type = "o", col = 100, pch=19)
  legend("bottomleft", 
    legend = c("Master", "Measurement"), 
    col = c(10,100),
    pch = c(17,19)
  )  
  
  # B
  plot(
    datMaster$b, type = "o", col = 10, pch=17, xlab="Colorspot", ylab = "Color Level", main = "Channel Comparision, Master vs Measurement: b")
  lines(measrmnts_grouped_by_colorspot$b, type = "o", col = 100, pch=19)  
  legend("bottomleft", 
    legend = c("Master", "Measurement"), 
    col = c(10,100),
    pch = c(17,19)
  )




## ----echo=FALSE, message=FALSE, warning=FALSE, results='asis', fig.align="center", out.width = '100%'----

  plot3d( MasterColorCard[,9:11][,1], MasterColorCard[,9:11][,2], MasterColorCard[,9:11][,3], col = "red" , xlab="Luminance", ylab="a(green-red)", zlab="b(blue-yellow)", size = 5)
  
  plot3d( avg_color_sheet[,1:3][,1], avg_color_sheet[,1:3][,2], avg_color_sheet[,1:3][,3], add=TRUE,size=5)
  


## ----echo=FALSE, message=FALSE, warning=FALSE, results='asis', fig.align="center", out.width = '100%'----

  abs_diff<-avg_color_sheet[,1:3]-MasterColorCard[,9:11]

  plot3d( MasterColorCard[,9:11][,1], MasterColorCard[,9:11][,2],
          MasterColorCard[,9:11][,3], xlab="Luminance", ylab="a", zlab="b", type="s", radius = apply((abs_diff),MARGIN=1,FUN=mean),lit=TRUE)
  


## ----echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '80%'----


###How do color chanels behave relative to the print master?
  
# par(mfrow=c(1,3))

plot(hexbin(avg_color_sheet$L,MasterColorCard$L),  xlab="Avg Measurement L", ylab="Master L", main="Colors Behaving w.r.t Master")
# plot(hexbin(avg_color_sheet$a,MasterColorCard$a),  xlab="Avg Measurement a", ylab="Master a", main="Colors Behaving w.r.t Master")
# plot(hexbin(avg_color_sheet$b,MasterColorCard$b),  xlab="Avg Measurement b", ylab="Master b", main="Colors Behaving w.r.t Master")



## ----echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '100%'----

###How does color dispersion behave?
delta_l_byTargetRow<-apply(pos@delta_l, 1, mean)
delta_a_byTargetRow<-apply(pos@delta_a, 1, mean)
delta_b_byTargetRow<-apply(pos@delta_b, 1, mean)
delta_e_byTargetRow<-apply(pos@delta_e, 1, mean)



index<-seq(1:7)

par(mfrow=c(2,2))

symbols(index, 2*delta_l_byTargetRow/max(delta_l_byTargetRow), circles=sqrt(delta_l_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Row", ylab="Error L", main="Row Cards vs Delta L")
symbols(index, 2*delta_a_byTargetRow/max(delta_a_byTargetRow), circles=sqrt(delta_a_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Row", ylab="Error a", main="Row Cards vs Delta a")
symbols(index, 2*delta_b_byTargetRow/max(delta_b_byTargetRow), circles=sqrt(delta_b_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Row", ylab="Error b",main="Row Cards vs Delta b")
symbols(index, 2*delta_e_byTargetRow/max(delta_e_byTargetRow), circles=sqrt(delta_e_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Row", ylab="Delta E", main="Row Cards vs Delta E")




## ----echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '100%'----
###How does color dispersion behave?
delta_l_byTargetRow<-apply(pos@delta_l, 2, mean)
delta_a_byTargetRow<-apply(pos@delta_a, 2, mean)
delta_b_byTargetRow<-apply(pos@delta_b, 2, mean)
delta_e_byTargetRow<-apply(pos@delta_e, 2, mean)



index<-seq(1:6)

par(mfrow=c(2,2))

symbols(index, 2*delta_l_byTargetRow/max(delta_l_byTargetRow), circles=sqrt(delta_l_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Column", ylab="Error L", main="Column Cards vs Delta L")
symbols(index, 2*delta_a_byTargetRow/max(delta_a_byTargetRow), circles=sqrt(delta_a_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Column", ylab="Error a", main="Column Cards vs Delta a")
symbols(index, 2*delta_b_byTargetRow/max(delta_b_byTargetRow), circles=sqrt(delta_b_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Column", ylab="Error b", main="Column Cards vs Delta b")
symbols(index, 2*delta_e_byTargetRow/max(delta_e_byTargetRow), circles=sqrt(delta_e_byTargetRow/100), bg="grey50", inches=FALSE, xlab="Column", ylab="Delta E", main="Column Cards vs Delta E")



## ----cars, echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '80%'----

generate_nxn_matrix<-function(delta)
{
  delta_e_matrix<-matrix(, nrow = 8, ncol = 8)
  lower<-1
  for (i in 1:8){
    upper<-lower+7
    temp<-delta[lower:upper]
    if (i<2){
      delta_e_matrix<-rbind(temp)
    }
    else{
      #print(temp)
      delta_e_matrix=rbind(delta_e_matrix,temp)
      #print(delta_e_matrix)
    }
    lower<-lower+8
  }
  return (delta_e_matrix)
}

generate_color_palette<-function(delta, delta_title)
{
  library(plotrix)
  
  #Build the matrix data to look like a correlation matrix
  n <- 8
  x <- delta/max(delta)
  xmin <- 0
  xmax <- 1
  #for (i in 1:n) x[i, i] <- 1.0 #Make the diagonal all 1's
  
  #Generate the palette for the matrix and the legend.  Generate labels for the legend
  palmat <- color.scale(x, c(1, 0.4), c(1, 0.4), c(0.86, 1))
  palleg <- color.gradient(c(1, 0.4), c(1, 0.4), c(0.86, 1), nslices=100)
  lableg <- c(formatC(xmin, format="f", digits=2), formatC(1*(xmax-xmin)/4, format="f", digits=2), formatC(2*(xmax-xmin)/4, format="f", digits=2), formatC(3*(xmax-xmin)/4, format="f", digits=2), formatC(xmax, format="f", digits=2))
  
  #Set up the plot area and plot the matrix
  par(mar=c(5, 5, 5, 8))
  color2D.matplot(x, cellcolors=palmat, main=paste(delta_title, " on ", n, " X ", n, " Color Spots", sep=""), show.values=2, vcol=rgb(0,0,0), axes=FALSE, vcex=0.7)
  axis(2, at=seq(1, n, 1)-0.5, labels=seq(n, 1, -1), tck=-0.01, padj=-1)
  axis(1, at=seq(1, n, 1)-0.5, labels=seq(1, n, 1), tck=-0.01, padj=0.7)
  
  #Plot the legend
  pardat <- par()
  color.legend(pardat$usr[2]+0.5, 0, pardat$usr[2]+1, pardat$usr[2], paste(" ", lableg, sep=""), palleg, align="rb", gradient="y", cex=0.7)
  
}

###How does color dispersion behave?
delta_l<-getColorDifference(MasterColorCard,avg_color_sheet,"L")
delta_a<-getColorDifference(MasterColorCard,avg_color_sheet,"a")
delta_b<-getColorDifference(MasterColorCard,avg_color_sheet,"b")
delta_e<-getColorDifference(MasterColorCard,avg_color_sheet,"all")

# par(mfrow=c(2,2))

delta_e_matrix<-generate_nxn_matrix(delta_e)
generate_color_palette(delta_e_matrix,"Delta E (Avg Lab)")

delta_L_matrix<-generate_nxn_matrix(delta_l)
generate_color_palette(delta_L_matrix,"Delta L (Avg L)")

delta_a_matrix<-generate_nxn_matrix(delta_a)
generate_color_palette(delta_a_matrix,"Delta a (Avg a)")

delta_b_matrix<-generate_nxn_matrix(delta_b)
generate_color_palette(delta_b_matrix,"Delta b (Avg b) ")




## ----echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '100%'----

###Which color chanel contributes more to error?
impact_on_error<-function(delta_e, delta, x_title, ..){
  
  # We create 2 vectors x and y. It is a polynomial function.
  x <- delta
  y <-  delta_e
  
  # Basic plot of x and y :
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3 , xlab=x_title , ylab="Delta E", main=paste("Impact of chanels on Delta E")) 
  
  # Can we find a polynome that fit this function ?
  model=lm(y ~ x + I(x^2) + I(x^3))
  
  # I can get the features of this model :
  summary(model)
  model$coefficients
  summary(model)$adj.r.squared
  
 
  suppressWarnings({ 
    
      #For each value of x, I can get the value of y estimated by the model, and the confidence interval around this value.
      myPredict <- predict( model , interval="predict" )
      
      #Finally, I can add it to the plot using the line and the polygon function with transparency.
      ix <- sort(x,index.return=T)$ix
      #lines(x[ix], myPredict[ix , 1], col=2, lwd=2 )  
      polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
    
    })


}

delta_e_flatten<-as.vector(t(pos@delta_e))
delta_l_flatten<-as.vector(t(pos@delta_l))
delta_a_flatten<-as.vector(t(pos@delta_a))
delta_b_flatten<-as.vector(t(pos@delta_b))

par(mfrow=c(1,3))

impact_on_error(delta_e_flatten, delta_l_flatten, "Delta L")
impact_on_error(delta_e_flatten, delta_a_flatten, "Delta a")
impact_on_error(delta_e_flatten, delta_b_flatten, "Delta b")
  


## ----echo=FALSE, message=FALSE, warning=FALSE, fig.align="center", out.width = '100%'----

###Which color chanel contributes more to error?
impact_on_error<-function(delta_e, delta, x_title, ..){
  
# We create 2 vectors x and y. It is a polynomial function.
x <- delta
y <-  delta_e

# Basic plot of x and y :
plot(x,y,col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3 , xlab=x_title , ylab="Delta E", main=paste("Impact of color chanels on Error")) 

# Can we find a polynome that fit this function ?
model=lm(y ~ x + I(x^2) + I(x^3))

# I can get the features of this model :
summary(model)
model$coefficients
summary(model)$adj.r.squared

 suppressWarnings({
    #For each value of x, I can get the value of y estimated by the model, and the confidence interval around this value.
    myPredict <- predict( model , interval="predict" )
    
    #Finally, I can add it to the plot using the line and the polygon function with transparency.
    ix <- sort(x,index.return=T)$ix
    #lines(x[ix], myPredict[ix , 1], col=2, lwd=2 )  
    polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
 })

}

par(mfrow=c(1,3))

impact_on_error(delta_e, delta_l, "Delta L")
impact_on_error(delta_e, delta_a, "Delta a")
impact_on_error(delta_e, delta_b, "Delta b")
  


## ----echo=FALSE, message=FALSE, warning=FALSE, fig4, fig.align="center", fig.height=3, out.width = '100%'----
# distribution comparision
# iterate measurments

  deltaE <- getColorDifferenceDf(datMaster, measrmnts_grouped_by_colorspot[-1:-2])
  
  par(mfrow=c(1,3))
  
  # L
  plot(density(datMaster$L), col="3", main="Distribution Comparision L")
  lines(density(measrmnts_grouped_by_colorspot$L), col="4")

  # A
  plot(density(datMaster$a), col="3", main="Distribution Comparision a")
  lines(density(measrmnts_grouped_by_colorspot$a), col="4")

  # B
  plot(density(datMaster$b), col="3", main="Distribution Comparision b")
  lines(density(measrmnts_grouped_by_colorspot$b), col="4")  

  

