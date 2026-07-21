# 1 SETUP ENVRIONMENT ----
## 1.1 Required packages ----
library(tidyverse)    #for data manipulation
library(sjPlot)       #for assumptions checks
library(car)          #for type III sums of squares ANOVA
library(MuMIn)        #for R-squared estimates from linear mixed effects models
library(effectsize)   #for effect size estimates
library(nlme)         #for statistical models
library(extrafont)    #for Arial font
library(cowplot)      #for plotting multi-panel figures
library(emmeans)      #for least-squares means
library(reghelper)    #for simple slope tests

## 1.2 Constants ----
setClrs=c("#32006E",       #FHL colour
          "#FFBF00",       #Bodega colour
          "#002664")       #Villanova colour
setShapes=c(21,            #FHL shape
            22,            #Bodega shape
            24)            #Villanova shape
setLineTypes=c("solid",    #FHL shape
               "dashed",   #Bodega shape
               "dotted")   #Villanova shape

## 1.3 Utility functions ----
cohensd <- function(mod,modType){
  # Computes Cohen's d for Strain Rate (1/s) from a t-value.
  #
  # INPUT
  # mod: output from summary() of a model object.
  # modType: designation of either summary of an lme object from nlme ("lme"),
  #          or an lm object from base ("lm").
  #
  # OUTPUT
  # cohensd: an "effectsize_table" with the cohen's d value for your input 
  #          t-value and the associated 95% confidence intervals.
  #
  # EXAMPLE CALL
  # mod.ten.d=cohensd(mod.ten.lmO,"lme")
  #
  # Written by: Keegan Lutek
  
  if (modType=="lme"){
    t_to_d(mod$tTable["rate","t-value"],          #estimate effect size for rate from t-value
           mod$tTable["rate","DF"]) %>%           #pass degrees of freedom for t-value
      mutate(CI_diff=d-CI_low)                    #calculate DI size
  } else if (modType=="lm"){
    t_to_d(mod$coefficients["rate","t value"],
           mod$df[[2]]) %>%
      mutate(CI_diff=d-CI_low)
  }
}
genNewData <- function(rawData,mod,newVar,plotGlobal=FALSE){
  # Generate dataset with higher resolution along rate for linear model predictions.
  #
  # INPUT
  # rawData: datset used to generate mod. Must have variables "rate" (Strain Rate) 
  #          and "set" (Animal Set).
  # mod: linear mixed effects model from nlme with "rate" (Strain Rate) and "set"
  #      (Animal Set) as predictors.
  # newVar: name of your predictions in the output. Character variable.
  #
  # OPTIONAL INPUT
  # plotGlobal: indicator whether the predictions should be by Animal Set or global.
  #             Default=TRUE.
  #
  # OUTPUT
  # genNewData: a dataframe with your predictions (as named in the input), their
  #             standard error (when available, SE), the Strain Rate at which they
  #             were calculated (rate) and, if plotGlobal==FALSE, the Animal Set 
  #             for each prediction.
  #
  # EXAMPLE CALL
  # pred.ten=genNewData(tfDPull,mod.ten,"ten.p",plotGlobal=TRUE)
  #
  # Written by: Keegan Lutek
  
  if (nlevels(rawData$set)==3){                     #i.e. if you have all Animal Sets in the dataset
    
    if (!plotGlobal) {                              #if you're plotting predictions for each Animal Set
      dataout=rawData %>%
        reframe(
          rate=seq(min(rate),
                   max(rate),
                   length.out=100),
          .by=set) %>%                              #generate fine-resolution x across the range for each Animal Set
        mutate(!!sym(newVar):=predict(mod,
                                      newdata=.,
                                      level=0),     #generate predictions at those x values
               set=fct_relevel(set,
                               "WA2","CA8")) %>%    #relevel Animal Set so it matches raw data
        select(set,rate,!!sym(newVar))           #select only necessary variables
      
    } else {                                        #if you're plotting global predictions
      
      newX=seq(min(rawData$rate),
               max(rawData$rate),
               length.out=100)
      
      dataout=emmip(mod, ~rate, 
                    at=list(rate=newX),
                    plotit=FALSE) %>%
        mutate(!!sym(newVar):= yvar) %>%
        select(rate,!!sym(newVar),SE)
    }
    
  } else {                                         #if there are only two Animal Sets in the dataset
    
    if (!plotGlobal){
      dataout=rawData %>%
        reframe(
          rate=seq(min(rate),
                   max(rate),
                   length.out=100),
          .by=set) %>%
        mutate(!!sym(newVar):=predict(mod,
                                      newdata=.,
                                      level=0),
               set=fct_relevel(set,
                               "WA2")) %>%
        select(set,rate,!!sym(newVar))
      
    } else {
      
      newX=seq(min(rawData$rate),
               max(rawData$rate),
               length.out=100)
      
      dataout=emmip(mod,~rate, 
                    at=list(rate=newX),
                    plotit=FALSE) %>%
        mutate(!!sym(newVar):=yvar) %>%
        select(rate,!!sym(newVar),SE)
      
    }
  }
  
  return(dataout)
  
}
genNewData.glm <- function(rawData,mod,newVar){
  # Generate dataset with higher resolution along rate for generalized linear 
  # model predictions.
  #
  # INPUT
  # rawData: datset used to generate mod. Must have variables "rate" (Strain Rate)
  #          and "set" (Animal Set).
  # mod: linear mixed effects model from nlme with "rate" (Strain Rate) and "set"
  #      (Animal Set) as predictors.
  # newVar: name of your predictions in the output. Character variable.
  #
  # OUTPUT
  # genNewData: a dataframe with your predictions (as named in the input),and the 
  #             Strain Rate at which they were calculated (rate).
  #
  # EXAMPLE CALL
  # pred.adh=genNewData.glm(waPull,mod.adh,"adh.p")
  #
  # Written by: Keegan Lutek
  
  dataout=rawData %>%
    reframe(
      rate=seq(min(rate),
               max(rate),
               length.out=100),
      .by=set) %>%                                  #generate fine-resolution x across the range for each Animal Set
    mutate(!!sym(newVar):=predict(mod,
                                  newdata=.,
                                  level=0,            #generate predictions at those x values
                                  type="response"),
           set=fct_relevel(set,
                           "WA2","CA8")) %>%          #relevel Animal Set so it matches raw data
    select(set,rate,!!sym(newVar))
}
get_legend2 <- function(plot, legend=NULL) {
  # A replacement for cowplot's get_legend that is compatible with
  # the latest version of ggplot.
  #
  # INPUT
  # plot: a ggplot or grob object
  #
  # OUTPUT
  # the legend from your input
  #
  # EXAMPLE CALL
  # legend <- get_legend2(plot.adh)
  #
  # Written by: teunbrand and Tim Bainbridge (from https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0)
  
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what="zeroGrob", 
    FUN.VALUE=logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}
makeTable <- function(varName,modType,N,predictors,mod.aO,mod.lmO,mod.R2){
  # Make the properly formatted row(s) of Table 2 and Table S1 from the 
  # required statistical output.
  # 
  # INPUT
  # varName: the name of your variable of interest. Character variable.
  # modType: indicator for the type of model input later. Options are "lme",
  #          "glm", or "lm". Character variable.
  # N: the sample size for your variable of interest.
  # predictors: the name(s) of your predictors in the model. Input as a 
  #             concatenated set of character variables.
  # mod.aO: anova-type output for your model from anova.lme, if linear
  #         mixed effects model, or from Anova, if generalized linear model or
  #         linear model.
  # mod.lmO: summary of your model.
  # mod.R2: r-squared value for your model. From MuMIn::r.squaredGLMM if linear
  #         mixed effects model.
  #
  # OUTPUT
  # a table with columns for the response variable, form of the model, type of model
  # predictors in the model, F-value, degrees of freedom, and p-value for each predictor
  # and the R-squared value(s).
  #
  # EXAMPLE CALL
  #table.ten=makeTable("ten","lme",
  #                    nlevels(tfDPull$ind),
  #                    c("rate","set"),
  #                    mod.ten.aO,mod.ten.lmO,
  #                    mod.ten.R2)
  #
  # Written by: Keegan Lutek
  
  # Extract model-type specific information
  if (modType=="lme"){
    
    # Fixed effects model formula
    fitForm=as.character(mod.lmO$call[["fixed"]])
    # R-squared value(s)
    R2=paste0(round(mod.R2[1],2),"[",round(mod.R2[2],2),"]")
    # Name for the F-value in the model object
    Fstring="F-value"
    # Name for the p-values in the model object
    pstring="p-value"
    
  } else if (modType=="glm"){
    
    fitForm=as.character(mod.lmO$call[["formula"]])
    R2=round(mod.R2,2)
    Fstring="F values"
    pstring="Pr(>F)"
    
  } else if (modType=="lm"){
    
    fitForm=as.character(mod.lmO$call[["formula"]])
    R2=round(mod.R2,2)
    Fstring="F value"
    pstring="Pr(>F)"
    
  } #end if statement from line 247
  
  # Create variable for output.
  Predictors=predictors
  
  # Preset variables with only one filled row
  Response=rep("",length(predictors))
  Form=rep("",length(predictors))
  Type=rep("",length(predictors))
  nInd=rep("",length(predictors))
  Fval=rep("",length(predictors))
  df=rep("",length(predictors))
  p=rep("",length(predictors))
  R2mR2c=rep("",length(predictors))
  
  # Iterate through the predictors and extract the remaining information.
  count=1
  for (i in predictors){
    if (count==1){ #For all columns
      # Response variable
      Response[count]=varName
      # Formula for the fixed effects model
      Form[count]=fitForm[[3]]
      # Type of model
      Type[count]=modType
      # Sample size
      nInd[count]=N
      # F-value for each predictor
      Fval[count]=round(mod.aO[predictors[count],Fstring],2)
      # Degrees of freedom for each predictor
      if (modType=="lme"){
        df[count]=paste0(mod.aO[predictors[count],"numDF"],",",mod.aO[predictors[count],"denDF"])
      } else if (modType %in% c("glm","lm")){
        df[count]=paste0(mod.aO[predictors[count],"Df"],",",mod.aO["Residuals","Df"])
      }
      # P-value for each predictor
      p[count]=round_p(mod.aO[predictors[count],pstring])
      # R-squared value
      R2mR2c[count]=R2
    } else { #For columns that require more than one row
      Response[count]=""
      Form[count]=""
      Type[count]=""
      nInd[count]=""
      Fval[count]=round(mod.aO[predictors[count],Fstring],2)
      if (modType=="lme"){
        df[count]=paste0(mod.aO[predictors[count],"numDF"],",",mod.aO[predictors[count],"denDF"])
      } else if (modType %in% c("glm","lm")){
        df[count]=paste0(mod.aO[predictors[count],"Df"],",",mod.aO["Residuals","Df"])
      }
      p[count]=round_p(mod.aO[predictors[count],pstring])
      R2mR2c[count]=""
    }
    count=count+1
  }
  
  # Combine into single dataframe
  out=data.frame(Response,
                 Form,
                 Type,
                 nInd,
                 Predictors,
                 Fval,
                 df,
                 p,
                 R2mR2c)
}
plotAnnot <- function(mod.aO,mod.R2,mod.d){
  # Generate annotations for linear mixed effects plots.
  #
  # INPUT
  # mod.aO: ANOVA-type output from anova.lme.
  # mod.R2: Estimated R-squared values from r.squaredGLMM.
  # mod.d: Cohen's d estimate data from cohensd.
  #
  # OUTPUT
  # an expression with the above information formatted for plot overlay.
  #
  # EXAMPLE CALL
  # plot.annot=plotAnnot(mod.ten.aO,mod.ten.R2,mod.ten.d)
  #
  # Written by: Keegan Lutek
  
  #p-value for Strain Rate
  p.rate=as.expression(bquote(p[rate]==.(round_p(mod.aO["rate","p-value"]))))
  #p-value for Animal Set
  p.set=as.expression(bquote(p[set]==.(round_p(mod.aO["set","p-value"]))))
  #R-squared for whole model
  R2=as.expression(bquote(R^2 == .(round(mod.R2[1],2))~"["*.(round(mod.R2[2],2))*"]"))
  #effect size for rate
  d=as.expression(bquote(d[rate] == .(round(mod.d$d,2))*"±"*.(round(mod.d$CI_diff,2))))
  
  out=c(p.rate,p.set,R2,d)
  
  return(out)
}
plotAnnot.glm <- function(mod.aO,mod.R2,mod.d){
  # Generate annotations for generalized linear model plots.
  #
  # INPUT
  # mod.aO: ANOVA-type output from anova.lme.
  # mod.R2: Estimated R-squared values from r.squaredGLMM.
  # mod.d: Cohen's d estimate data from cohensd.
  #
  # OUTPUT
  # an expression with the above information formatted for plot overlay.
  #
  # EXAMPLE CALL
  # plot.annot=plotAnnot.glm(mod.adh.aO,mod.adh.R2,mod.adh.d)
  #
  # Written by: Keegan Lutek
  
  #p-value for Strain Rate
  p.rate=as.expression(bquote(p[rate]==.(round_p(mod.aO["rate","Pr(>F)"]))))
  #p-value for Animal Set
  p.set=as.expression(bquote(p[set]==.(round_p(mod.aO["set","Pr(>F)"]))))
  #R-squared for whole model
  R2=as.expression(bquote(R^2 == .(round(mod.R2[1],2))))
  #effect size for rate
  d=as.expression(bquote(d[rate] == .(round(mod.d$d,2))*"±"*.(round(mod.d$CI_diff,2))))
  
  out=c(p.rate,p.set,R2,d)
  
  return(out)
}
plotLME <- function(rawData,yVar,preds,yPred,
                    xCoords,yCoords,statAnnot,
                    yTitle,clrs,shapeType,lineType,
                    ptSize=2,plotGlobal=FALSE,plotLines=TRUE){
  # Generate a plot for a linear mixed effects model.
  #
  # INPUT
  # rawData: data frame that contains the raw data to be plotted. Must include the columns
  #          "rate" (Strain Rate) and "set" (Animal Set).
  # yVar: a character variable containing the name of your y-variable.
  # preds: a dataframe with your predictions. Must contain the variable "rate" (Strain Rate).
  #        If plotGlobal==FALSE (see below), then this dataframe must also contain the variable
  #        "set" (Animal Set).
  # yPred: a character variable containing the name of your y-variable predictions.
  # xCoords: the x-coordinate for the position of your annotations (from plotAnnot(.glm)).
  # yCoords: the y-coordinate for the positions of your annotations (from plotAnnot(.glm)).
  # statAnnot: the output from plotAnnot(.glm).
  # yTitle: a character variable with the title for your y-axis.
  # clrs: a character variable with the colours to be used for plotting.
  # shapeType: a vector with the number codes for the shapes to be used for plotting. See
  #            ggplot documentation for options.
  # lineType: a character variable with the line types to be used for plotting. See ggplot
  #           documentation for options.
  #
  # OPTIONAL INPUT
  # ptSize: a one-element vector with the size of points you would like to use for your plot.
  #         Default=2.
  # plotGlobal: indicator whether the predictions should be by Animal Set or global.
  #             Default=TRUE.
  #
  # EXAMPLE CALL
  #plot.ten=plotLME(tfDPull,"ten",pred.ten,"ten.p",
  #                 4.2,c(0.57,0.53,0.49,0.45),plot.annot,
  #                 "Tube Foot Disc Tenacity (MPa)",setClrs,setShapes,
  #                 setLineTypes,plotGlobal=TRUE)
  #
  # Written by: Keegan Lutek
  
  plot.out=ggplot()+     
    geom_point(data=rawData,                                        #dataset from above
               aes(x=rate,y=!!sym(yVar),shape=set,color=set),       #set aesthetics
               size=ptSize,alpha=0.4,                               #set size & transparency of point
               stroke=1)                                            #set stroke width
  if (plotLines) {
    if (!plotGlobal) {
      plot.out=plot.out+
        geom_line(data=preds,                                       #prediction dataset
                  aes(x=rate,y=!!sym(yPred),color=set,linetype=set),  #set aesthetics
                  linewidth=1.3,lineend="round")                      #set line width & style
    } else {
      plot.out=plot.out+
        geom_line(data=preds,
                  aes(x=rate,y=!!sym(yPred)),
                  color="black",linewidth=1.3,lineend="round")      
    }
  }
  
  plot.out=plot.out+
    scale_shape_manual(values=shapeType,                            #set shape of points
                       name="Set:")+                                #name for legend
    scale_color_manual(values=clrs,                                 #set colour values
                       name="Set:")                                 #name for legend
  
  if (plotLines) {
    plot.out=plot.out+
      scale_linetype_manual(values=lineType,                        #set linetype
                            name="Set:")                            #name for legend
  }
  
  plot.out=plot.out+
    annotate("text",                                                #text annotations
             x=xCoords, y=yCoords,                                  #annotation position
             label=statAnnot,                                       #labels from above
             hjust=1,
             size=8/.pt)+                                           #align right edge
    ylab(yTitle)+
    xlab(expression("Strain Rate (s"^{-1}*")"))+
    theme_JEB()+                                                    #adjust aesthetics
    theme(legend.position="top",                                    #move legend to top
          axis.title.y=element_text(margin=margin(r=15,unit="pt")), #y-axis title position
          plot.margin=margin(r=3,l=8,unit="pt"))                    #adjust margins around plot
  
  plot.out
}
plotGeneral <- function(rawData,xVar,yVar,preds,yPred,
                        xCoords,yCoords,statAnnot,
                        xTitle,yTitle,hjustVal,
                        clrs,shapeType,lineType,ptSize=2){
  # Generate a plot for a (generalized) linear model. This is effectively the same as for
  # plotLME, but modifies where relevant information is pulled from.
  #
  # INPUT
  # rawData: data frame that contains the raw data to be plotted. Must include the columns
  #          "rate" (Strain Rate) and "set" (Animal Set).
  # yVar: a character variable containing the name of your y-variable.
  # preds: a dataframe with your predictions. Must contain the variable "rate" (Strain Rate).
  #        If plotGlobal==FALSE (see below), then this dataframe must also contain the variable
  #        "set" (Animal Set).
  # yPred: a character variable containing the name of your y-variable predictions.
  # xCoords: the x-coordinate for the position of your annotations (from plotAnnot(.glm)).
  # yCoords: the y-coordinate for the positions of your annotations (from plotAnnot(.glm)).
  # statAnnot: the output from plotAnnot(.glm).
  # yTitle: a character variable with the title for your y-axis.
  # clrs: a character variable with the colours to be used for plotting.
  # shapeType: a vector with the number codes for the shapes to be used for plotting. See
  #            ggplot documentation for options.
  # lineType: a character variable with the line types to be used for plotting. See ggplot
  #           documentation for options.
  #
  # OPTIONAL INPUT
  # ptSize: a one-element vector with the size of points you would like to use for your plot.
  #         Default=2.
  #
  # OUTPUT
  # A ggplot object with your plotted data and predictions.
  #
  # EXAMPLE CALL
  #plot.fp=plotGeneral(waPull,"fp","adh",pred.fp,"fp.p",
  #                    0.01,seq(2.5,1.7,length=2),plot.annot,
  #                    expression("Footprints (cm"^{-2}*")"),
  #                    expression("Whole Animal Adhesion (N cm"^{-1}*")"),
  #                    0,setClrs,setShapes,setLineTypes)
  #
  # Written by: Keegan Lutek
  
  plot.out=ggplot()+     
    geom_point(data=rawData,                                         #dataset from above
               aes(x=!!sym(xVar),y=!!sym(yVar),
                   shape=set,color=set),                             #set aesthetics
               size=ptSize,alpha=0.4,                                #set size & transparency of point
               stroke=1)+                                            #set stroke width
    geom_line(data=preds,                                            #prediction dataset
              aes(x=!!sym(xVar),y=!!sym(yPred),
                  color=set,linetype=set),                           #set aesthetics
              linewidth=1.3,lineend="round")+                        #set line width & style
    scale_shape_manual(values=shapeType,                             #set shape of points
                       name="Set:")+                                 #name for legend
    scale_color_manual(values=clrs,                                  #set colour values
                       name="Set:")+                                 #name for legend
    scale_linetype_manual(values=lineType, #set linetype
                          name="Set:")+                              #name for legend
    annotate("text",                                                 #text annotations
             x=xCoords, y=yCoords,                               #annotation position
             label=statAnnot,                                      #labels from above
             hjust=hjustVal,
             size=8/.pt)+                                            #align right edge
    ylab(yTitle)+
    xlab(xTitle)+
    theme_JEB()+                                                     #adjust aesthetics
    theme(legend.position="top",                                     #move legend to top
          axis.title.y=element_text(margin=margin(r=15,unit="pt")),  #y-axis title position
          plot.margin=margin(r=3,l=8,unit='pt'))                     #adjust margins around plot
  
  plot.out                                                           #display the generated plot
}
plotInset <- function(emmFrame,expClrs){
  
  nColours=length(unique(emmFrame$set))
  
  plot.out=ggplot()+
    geom_point(data=emmFrame,
               aes(x=set,y=emmean,colour=set),
               size=3)+
    geom_errorbar(data=emmFrame,
                  aes(x=set,ymin=emmean-SE,ymax=emmean+SE,colour=set),
                  width=0.2)
    
    if (nColours==2) {
    plot.out=plot.out+
      scale_colour_manual(values=expClrs[c(1,3)])
    } else {
      plot.out=plot.out+
        scale_colour_manual(values=expClrs)
    }
    
  plot.out=plot.out+
    theme_JEB()+                                                    #adjust aesthetics
    theme(legend.position="top",                                    #move legend to top
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=margin(r=3,l=8,unit="pt"),
          plot.background=element_rect(fill='transparent',
                                       colour='transparent'))
}
round_p <- function(number){
  # A function to round your p-value to an appropriate representation for reporting in 
  # scientific journals.
  #
  # INPUT
  # number: your number to be rounded.
  # 
  # OUTPUT
  # The input number rounded (or changed to a relevant inequality). Returned as a character.
  #
  # EXAMPLE CALL
  # roundedNum=round_p(0.234)
  #
  # Written by: Keegan Lutek
  
  # Start by rounding your number to two significant digits
  r.number=signif(number,digits=2)
  
  # If the value is less than 0.0001, replace with a relevant inequality
  if (r.number<0.0001){
    r.number="<0.0001"
    
    # If the value is less than 0.001, replace with a relevant inequality
  } else if (r.number<0.001){
    r.number="<0.001"
  }
  
  # Change the number to a character variable.
  r.number=as.character(r.number)
  
  return(r.number)
}
theme_JEB <- function(){ 
  # An alternate theme for ggplot2 that aligns with JEB's style guide. Use at end of
  # ggplot calls to adjust the format accordingly.
  
  font="Arial"                                   #assign font family up front
  
  theme_classic() %+replace%                     #replace elements we want to change
    
    theme(plot.title=element_blank(),            #get rid of the plot title
          plot.subtitle=element_blank(),         #get rid of the plot subtitle
          plot.caption=element_blank(),          #get rid of the caption
          axis.title=element_text(family=font,   #modify the axis title formatting
                                  size=8),
          axis.text=element_text(family=font,    #modify the axis text formatting
                                 size=8),
          legend.text=element_text(family=font,  #modify the legend text formatting
                                   size=8),
          legend.title=element_text(family=font, #modify the legend title text.
                                    size=8))
}

## 1.4 Setup for Saving PDF Plots ----
# We have used GhostScript to save pdf plots reliably. The path below needs to 
#  be adjusted for your computer, and GhostScript can be downloaded here:
#  https://ghostscript.com/releases/gsdnld.html
Sys.setenv(R_GSCMD="C:/Program Files/gs/gs10.04.0/bin/gswin64c.exe")
loadfonts()              #load in the fonts from the extrafonts package

# 2 DATA ANALYSIS ----
## 2.1 Disc Tenacity (MPa) ----
### 2.1.1 Load in dataset ----
tfDPull=read_csv("tfDPull.csv",col_types='ffnn') %>%
  rename(ind=Individual,                           #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         ten=`Tube Foot Disc Tenacity (MPa)`) %>%
  mutate(set=fct_relevel(set,                      #order levels by time removed from ocean
                         "WA2")) %>%
  drop_na()                                        #drop rows with missing NaNs
### 2.1.2 Statistical analysis ----
mod.ten=lme(data=tfDPull,                                    #dataset from above
            ten~rate+set,                                    #model form: fixed effects
            random=~1|ind)                                   #model form: random effects
mod.ten.assumptions=plot_model(mod.ten,type="diag")          #generate plots for assumptions checks
mod.ten.assumptions[[1]]                                     #check normality (qq plot)
mod.ten.assumptions[[2]]                                     #check normality (histogram)
mod.ten.assumptions[[3]]                                     #check homoskedasticity
mod.ten.lmO=summary(mod.ten)                                 #print model summary
mod.ten.aO=anova.lme(update(mod.ten,                         #model from above
                            contrasts=list(set=contr.sum)),  #update with contrasts for type III sums of squares
                     type="marginal")                        #request type III sums of squares
mod.ten.R2=r.squaredGLMM(mod.ten)                            #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.ten.d=cohensd(mod.ten.lmO,"lme")                         #calculate Cohen's d
pred.ten=genNewData(tfDPull,mod.ten,"ten.p",plotGlobal=TRUE) #generate model predictions


table.ten=makeTable("ten","lme",                             #create the row(s) of the table
                    nlevels(tfDPull$ind),
                    c("rate","set"),
                    mod.ten.aO,mod.ten.lmO,
                    mod.ten.R2)

### 2.1.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.ten.aO,mod.ten.R2,mod.ten.d)

#### b. plot ----
plot.ten=plotLME(tfDPull,"ten",pred.ten,"ten.p",
                 4.2,c(0.57,0.53,0.49,0.45),plot.annot,
                 "Tube Foot Disc Tenacity (MPa)",setClrs,setShapes,
                 setLineTypes,plotGlobal=TRUE)+
  annotate("text",
           x=4.25,
           y=0.57,
           label="*",
           hjust=0,
           size=8/.pt)

#### c. clean up environment ----
rm(tfDPull,mod.ten,mod.ten.assumptions,
   mod.ten.lmO,mod.ten.aO,mod.ten.R2,
   mod.ten.d,pred.ten,plot.annot)

## 2.2 Tube Foot Stem Breaking Force (N) ----
### 2.2.1 Load in dataset ----
tfSPull=read_csv("tfSPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                      #select variables necessary
         `Strain Rate (1/s)`,
         `Tube Foot Stem Breaking Force (N)`) %>%
  rename(ind=Individual,                               #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         bF=`Tube Foot Stem Breaking Force (N)`) %>%
  mutate(set=fct_relevel(set,                          #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                            #drop rows with NaNs
### 2.2.2 Statistical analysis ----
mod.bF=lme(data=tfSPull,                                     #dataset from above
           bF~rate+set,                                      #model form: fixed effects
           random=~1|ind)                                    #model form: random effects
mod.bF.assumptions=plot_model(mod.bF,type="diag")            #generate plots for assumptions checks
mod.bF.assumptions[[1]]                                      #check normality (qq plot)
mod.bF.assumptions[[2]]                                      #check normality (histogram)
mod.bF.assumptions[[3]]                                      #check homoskedasticity
mod.bF.lmO=summary(mod.bF)                                   #print model summary
mod.bF.aO=anova.lme(update(mod.bF,                           #model from above
                           contrasts=list(set=contr.sum)),   #update with contrasts for type III sums of squares
                    type="marginal")                         #request type III sums of squares
mod.bF.R2=r.squaredGLMM(mod.bF)                              #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.bF.d=cohensd(mod.bF.lmO,"lme")                           #calculate Cohen's d for Strain Rate (1/s)

pred.bF=genNewData(tfSPull,mod.bF,"bF.p")                    #generate model predictions

emmean.bF=data.frame(emmeans(mod.bF,
                  ~set))

table.bF=makeTable("bF","lme",                               #create the row(s) of the table
                   nlevels(tfSPull$ind),
                   c("rate","set"),
                   mod.bF.aO,mod.bF.lmO,
                   mod.bF.R2)
### 2.2.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.bF.aO,mod.bF.R2,mod.bF.d)

#### b. plot ----
plot.bF=plotLME(tfSPull,"bF",pred.bF,"bF.p",
                4.2,c(0.86,0.823,0.785,0.74),plot.annot,
                "Tube Foot Stem Breaking Force (N)",setClrs[c(1,3)],setShapes[c(1,3)],
                setLineTypes[c(1,3)],plotLines=FALSE)+
  annotate("text",
           x=4.25,
           y=0.823,
           label="*",
           hjust=0,
           size=8/.pt)

plot.bF.inset=plotInset(emmean.bF,expClrs)+
  annotate("text",
           x=c(1,2),
           y=c(0.58,0.45),
           label=c("a","b"),
           size=8/.pt)

#### c. clean up environment ----
rm(tfSPull,mod.bF,mod.bF.assumptions,
   mod.bF.lmO,mod.bF.aO,mod.bF.R2,
   mod.bF.d,pred.bF,plot.annot)

## 2.3 Tube Foot Stem Extension (m) ----
### 2.3.1 Load in dataset ----
tfSPull=read_csv("tfSPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                    #select variables necessary
         `Strain Rate (1/s)`,
         `Tube Foot Stem Extension (m)`) %>%
  rename(ind=Individual,                             #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         xT=`Tube Foot Stem Extension (m)`) %>%
  mutate(set=fct_relevel(set,                        #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                          #drop rows with NaNs
### 2.3.2 Statistical analysis ----
mod.xT=lme(data=tfSPull,                                     #dataset from above
           xT~rate+set,                                      #model form: fixed effects
           random=~1|ind)                                    #model form: random effects
mod.xT.assumptions=plot_model(mod.xT,type="diag")            #generate plots for assumptions checks
mod.xT.assumptions[[1]]                                      #check normality (qq plot)
mod.xT.assumptions[[2]]                                      #check normality (histogram)
mod.xT.assumptions[[3]]                                      #check homoskedasticity
mod.xT.lmO=summary(mod.xT)                                   #print model summary
mod.xT.aO=anova.lme(update(mod.xT,                           #model from above
                           contrasts=list(set=contr.sum)),   #update with contrasts for type III sums of squares
                    type="marginal")                         #request type III sums of squares
mod.xT.R2=r.squaredGLMM(mod.xT)                              #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.xT.d=cohensd(mod.xT.lmO,"lme")                           #calculate Cohen's d for Strain Rate (1/s)

pred.xT=genNewData(tfSPull,mod.xT,"xT.p")                    #generate model predictions

emmean.xT=data.frame(emmeans(mod.xT,
                             ~set))

table.xT=makeTable("xT","lme",                               #create the row(s) of the table
                   nlevels(tfSPull$ind),
                   c("rate","set"),
                   mod.xT.aO,mod.xT.lmO,
                   mod.xT.R2)
### 2.3.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.xT.aO,mod.xT.R2,mod.xT.d)

#### b. plot ----
plot.xT=plotLME(tfSPull,"xT",pred.xT,"xT.p",
                4.2,seq(0.12,0.09,length=4),plot.annot,
                "Tube Foot Stem Extension (m)",setClrs[c(1,3)],setShapes[c(1,3)],
                setLineTypes[c(1,3)],
                plotLines=FALSE)+
  theme(axis.title.y=element_text(margin=margin(r=25,unit="pt")))+   #adjust position of y-axis title
  annotate("text",
           x=4.25,
           y=0.11,
           label="*",
           hjust=0,
           size=8/.pt)

plot.xT.inset=plotInset(emmean.xT,expClrs)+
  annotate("text",
           x=c(1,2),
           y=c(0.05,0.042),
           label=c("a","b"),
           size=8/.pt)

#### c. clean up environment ----
rm(tfSPull,mod.xT,mod.xT.assumptions,
   mod.xT.lmO,mod.xT.aO,mod.xT.R2,
   mod.xT.d,pred.xT,plot.annot)

## 2.4 Tube Foot Stem Initial Spring Constant (N/m) ----
### 2.4.1 Load in dataset ----
tfSPull=read_csv("tfSPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                    #select variables necessary
         `Strain Rate (1/s)`,
         `Tube Foot Stem Initial Spring Constant (N/m)`) %>%
  rename(ind=Individual,                             #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         iSC=`Tube Foot Stem Initial Spring Constant (N/m)`) %>%
  mutate(set=fct_relevel(set,                        #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                          #drop rows with NaNs

### 2.4.2 Statistical analysis ----
mod.iSC=lme(data=tfSPull,                                    #dataset from above
            iSC~rate+set,                                    #model form: fixed effects
            random=~1|ind)                                   #model form: random effects
mod.iSC.assumptions=plot_model(mod.iSC,type="diag")          #generate plots for assumptions checks
mod.iSC.assumptions[[1]]                                     #check normality (qq plot)
mod.iSC.assumptions[[2]]                                     #check normality (histogram)
mod.iSC.assumptions[[3]]                                     #check homoskedasticity
mod.iSC.lmO=summary(mod.iSC)                                 #print model summary
mod.iSC.aO=anova.lme(update(mod.iSC,                         #model from above
                            contrasts=list(set=contr.sum)),  #update with contrasts for type III sums of squares
                     type="marginal")                        #request type III sums of squares
mod.iSC.R2=r.squaredGLMM(mod.iSC)                            #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.iSC.d=cohensd(mod.iSC.lmO,"lme")                         #calculate Cohen's d for Strain Rate (1/s)

pred.iSC=genNewData(tfSPull,mod.iSC,"iSC.p",plotGlobal=TRUE) #generate model predictions

table.iSC=makeTable("iSC","lme",                             #create the row(s) of the table
                    nlevels(tfSPull$ind),
                    c("rate","set"),
                    mod.iSC.aO,mod.iSC.lmO,
                    mod.iSC.R2)
### 2.4.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.iSC.aO,mod.iSC.R2,mod.iSC.d)

#### b. plot ----
plot.iSC=plotLME(tfSPull,"iSC",pred.iSC,"iSC.p",
                 4.2,seq(3,2.2,length=4),plot.annot,
                 "Tube Foot Stem Initial Spring Constant (N/m)",setClrs[c(1,3)],setShapes[c(1,3)],
                 setLineTypes[c(1,3)],
                 plotLines=FALSE)

#### c. clean up environment ----
rm(tfSPull,mod.iSC,mod.iSC.assumptions,
   mod.iSC.lmO,mod.iSC.aO,mod.iSC.R2,
   mod.iSC.d,pred.iSC,plot.annot)

## 2.5 Tube Foot Stem Final Spring Constant (N/m) ----
### 2.5.1 Load in dataset ----
tfSPull=read_csv("tfSPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                    #select variables necessary
         `Strain Rate (1/s)`,
         `Tube Foot Stem Final Spring Constant (N/m)`) %>%
  rename(ind=Individual,                             #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         fSC=`Tube Foot Stem Final Spring Constant (N/m)`) %>%
  mutate(set=fct_relevel(set,                        #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                          #drop rows with NaNs

### 2.5.2 Statistical analysis ----
mod.fSC=lme(data=tfSPull,                                     #dataset from above
            fSC~rate+set,                                     #model form: fixed effects
            random=~1|ind)                                    #model form: random effects
mod.fSC.assumptions=plot_model(mod.fSC,type="diag")           #generate plots for assumptions checks
mod.fSC.assumptions[[1]]                                      #check normality (qq plot)
mod.fSC.assumptions[[2]]                                      #check normality (histogram)
mod.fSC.assumptions[[3]]                                      #check homoskedasticity
mod.fSC.lmO=summary(mod.fSC)                                  #print model summary
mod.fSC.aO=anova.lme(update(mod.fSC,                          #model from above
                            contrasts=list(set=contr.sum)),   #update with contrasts for type III sums of squares
                     type="marginal")                         #request type III sums of squares
mod.fSC.R2=r.squaredGLMM(mod.fSC)                             #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.fSC.d=cohensd(mod.fSC.lmO,"lme")                          #calculate Cohen's d for Strain Rate (1/s)

pred.fSC=genNewData(tfSPull,mod.fSC,"fSC.p",plotGlobal=TRUE)  #generate model predictions

table.fSC=makeTable("fSC","lme",                              #create the row(s) of the table
                    nlevels(tfSPull$ind),
                    c("rate","set"),
                    mod.fSC.aO,mod.fSC.lmO,
                    mod.fSC.R2)
### 2.5.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.fSC.aO,mod.fSC.R2,mod.fSC.d)

#### b. plot ----
plot.fSC=plotLME(tfSPull,"fSC",pred.fSC,"fSC.p",
                 4.2,seq(80,58,length=4),plot.annot,
                 "Tube Foot Stem Final Spring Constant (N/m)",setClrs[c(1,3)],setShapes[c(1,3)],
                 setLineTypes[c(1,3)],
                 plotLines=FALSE)

#### c. clean up environment ----
rm(tfSPull,mod.fSC,mod.fSC.assumptions,
   mod.fSC.lmO,mod.fSC.aO,mod.fSC.R2,
   mod.fSC.d,pred.fSC,plot.annot)

## 2.6 Tube Foot Stem Work to Failure (J) ----
### 2.6.1 Load in dataset ----
tfSPull=read_csv("tfSPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                    #select variables necessary
         `Strain Rate (1/s)`,
         `Tube Foot Stem Work to Failure (J)`) %>%
  rename(ind=Individual,                             #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         wF=`Tube Foot Stem Work to Failure (J)`) %>%
  mutate(set=fct_relevel(set,                        #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                          #drop rows with NaNs

### 2.6.2 Statistical analysis ----
mod.wF=lme(data=tfSPull,                                    #dataset from above
           wF~rate+set,                                     #model form: fixed effects
           random=~1|ind)                                   #model form: random effects
mod.wF.assumptions=plot_model(mod.wF,type="diag")           #generate plots for assumptions checks
mod.wF.assumptions[[1]]                                     #check normality (qq plot)
mod.wF.assumptions[[2]]                                     #check normality (histogram)
mod.wF.assumptions[[3]]                                     #check homoskedasticity
mod.wF.lmO=summary(mod.wF)                                  #print model summary
mod.wF.aO=anova.lme(update(mod.wF,                          #model from above
                           contrasts=list(set=contr.sum)),  #update with contrasts for type III sums of squares
                    type="marginal")                        #request type III sums of squares
mod.wF.R2=r.squaredGLMM(mod.wF)                             #estimate R-squared values (conditional - fixed effects; marginal - full model)
mod.wF.d=cohensd(mod.wF.lmO,"lme")                          #calculate Cohen's d for Strain Rate (1/s)

pred.wF=genNewData(tfSPull,mod.wF,"wF.p")                   #generate model predictions

emmean.wF=data.frame(emmeans(mod.wF,
                             ~set))

table.wF=makeTable("wF","lme",                              #create the row(s) of the table
                   nlevels(tfSPull$ind),
                   c("rate","set"),
                   mod.wF.aO,mod.wF.lmO,
                   mod.wF.R2)
### 2.6.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot(mod.wF.aO,mod.wF.R2,mod.wF.d)

#### b. plot ----
plot.wF=plotLME(tfSPull,"wF",pred.wF,"wF.p",
                4.2,seq(0.07,0.051,length=4),plot.annot,
                "Tube Foot Stem Work to Failure (J)",setClrs[c(1,3)],setShapes[c(1,3)],
                setLineTypes[c(1,3)],
                plotLines=FALSE)+
  ylim(0,0.07)+
  annotate("text",
           x=4.25,
           y=0.06366667,
           label="*",
           hjust=0,
           size=8/.pt)

plot.wF.inset=plotInset(emmean.wF,expClrs)+
  scale_y_continuous(breaks=c(0.01,0.015,0.02))+
  annotate("text",
           x=c(1,2),
           y=c(0.022,0.015),
           label=c("a","b"),
           size=8/.pt)

#### c. clean up environment ----
rm(tfSPull,mod.wF,mod.wF.assumptions,
   mod.wF.lmO,mod.wF.aO,mod.wF.R2,
   mod.wF.d,pred.wF,plot.annot)

## 2.7 Whole Animal Adhesion (N/cm2)| Strain Rate (1/s) ----
### 2.7.1 Load in dataset ----
waPull=read_csv("waPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                   #select variables necessary
         `Strain Rate (1/s)`,
         `Whole Animal Adhesion (N/cm2)`) %>%
  rename(ind=Individual,                            #rename variables for ease of use
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         adh=`Whole Animal Adhesion (N/cm2)`) %>%
  mutate(set=fct_relevel(set,                       #order levels by time out of ocean
                         "WA2")) %>%
  drop_na()                                         #drop rows with NaNs
### 2.7.2 Statistical analysis ----
mod.adh=glm(data=waPull,adh~rate+set,family=Gamma(link="log"))
residualPlot(mod.adh)
mod.adh.lmO=summary(mod.adh)                                     #print model summary
mod.adh.aO=Anova(update(mod.adh,                                 #model from above
                        contrasts=list(set=contr.sum)),          #update with contrasts for type III sums of squares
                 type="III",                                     #request type III sums of squares
                 test.statistic="F")
mod.adh.R2=1-(mod.adh.lmO$deviance/mod.adh.lmO$null.deviance)
mod.adh.d=t_to_d(mod.adh.lmO$coefficients["rate","t value"],     #estimate effect size for rate from t-value
                 mod.adh.lmO$df.residual) %>%                    #pass degrees of freedom for t-value
  mutate(CI_diff=d-CI_low)                                       #calculate DI size

pred.adh=genNewData.glm(waPull,mod.adh,"adh.p")                  #generate model predictions

emmean.adh=data.frame(emmeans(mod.adh,
                             ~set,
                             type='response')) %>%
  mutate(emmean=response)

table.adh=makeTable("adh","glm",                                 #create the row(s) of the table
                    nlevels(waPull$ind),
                    c("rate","set"),
                    mod.adh.aO,mod.adh.lmO,
                    mod.adh.R2)

pairs(emmeans(mod.adh,~set,type="response"))                     #print post-hoc comparisons

### 2.7.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot.glm(mod.adh.aO,mod.adh.R2,mod.adh.d)

#### b. plot ----
plot.adh=plotLME(waPull,"adh",pred.adh,"adh.p",
                 10,seq(0.78,0.63,length=4),plot.annot,
                 expression("Whole Animal Adhesion (N cm"^{-1}*")"),
                 setClrs,setShapes,setLineTypes,
                 plotLines=FALSE)+
  scale_x_continuous(trans='log10')+                           #present x-axis on log scale because so many low strain rate tests.
  annotate("text",
           x=10.1190476,
           y=0.73,
           label="*",
           hjust=0,
           size=8/.pt)

plot.adh.inset=plotInset(emmean.adh,expClrs)+
  ylim(0,0.41)+
  annotate("text",
           x=c(1,2,3),
           y=c(0.41,0.22,0.16),
           label=c("a","a","b"),
           size=8/.pt)

#### c. clean up environment ----
rm(waPull,mod.adh,
   mod.adh.lmO,mod.adh.aO,mod.adh.R2,
   mod.adh.d,pred.adh,plot.annot)

## 2.8 Tube Foot Attachments (1/cm2) | Strain Rate (1/s) ----
### 2.8.1 Load in dataset ----
waPull=read_csv("waPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                   #select variables necessary
         `Footprints (1/cm2)`,
         `Broken Tube Feet (1/cm2)`,
         `Strain Rate (1/s)`) %>%
  rename(ind=Individual,                            #rename variables for ease of use
         set=`Animal Set`,
         fp=`Footprints (1/cm2)`,
         btf=`Broken Tube Feet (1/cm2)`,
         rate=`Strain Rate (1/s)`) %>%
  mutate(set=fct_relevel(set,                       #order levels by time out of ocean
                         "WA2"), 
         atch=fp+btf) %>%                           #calculate number of attachments
  drop_na()                                         #drop rows with NaNs

### 2.8.2 Statistical analysis ----
mod.atch=lm(data=waPull,
            atch~rate+set)
mod.atch.assumptions=plot_model(mod.atch,type="diag")          #generate plots for assumptions checks
mod.atch.assumptions[[2]]                                      #check normality (qq plot)
mod.atch.assumptions[[3]]                                      #check normality (histogram)
mod.atch.assumptions[[4]]                                      #check homoskedasticity
mod.atch.lmO=summary(mod.atch)                                 #print model summary
mod.atch.aO=Anova(update(mod.atch,                             #model from above
                         contrasts=list(set=contr.sum)),       #update with contrasts for type III sums of squares
                  type=3)                                      #request type III sums of squares
mod.atch.d=cohensd(mod.atch.lmO,"lm")                          #calculate Cohen's d for Strain Rate (1/s)

pred.atch=genNewData(waPull,mod.atch,"atch.p")                 #generate model predictions

emmean.atch=data.frame(emmeans(mod.atch,
                              ~set))

table.atch=makeTable("atch","lm",                              #create the row(s) of the table
                     nlevels(waPull$ind),
                     c("rate","set"),
                     mod.atch.aO,mod.atch.lmO,
                     mod.atch.lmO$r.squared)

pairs(emmeans(mod.atch,~set,type="response"))                  #print post-hoc comparisons

### 2.8.3 Generate plot ----
#### a. generate annotations for plot ----
plot.annot=plotAnnot.glm(mod.atch.aO,mod.atch.lmO$r.squared,mod.atch.d)

#### b. plot ----
plot.atch=plotLME(waPull,"atch",pred.atch,"atch.p",
                  10,seq(5.9,4.8,length=4),plot.annot,
                  expression("Number of Attachments (cm"^{-2}*")"),
                  setClrs,setShapes,setLineTypes,
                  plotLines=FALSE)+
  scale_x_continuous(trans='log10')+                                 #present x-axis on log scale because so many low strain rate tests.
  annotate("text",
           x=10.1190476,
           y=5.533333,
           label="*",
           hjust=0,
           size=8/.pt)

plot.atch.inset=plotInset(emmean.atch,expClrs)+
  annotate("text",
           x=c(1,2,3),
           y=c(4.1,2.1,2),
           label=c("a","b","b"),
           size=8/.pt)

#### c. clean up environment ----
rm(waPull,mod.atch,mod.atch.assumptions,
   mod.atch.lmO,mod.atch.aO,
   mod.atch.d,pred.atch,plot.annot)

## 2.9 Whole Animal Adhesion (N/cm2) | Footprints (1/cm2) ----
### 2.9.1 Load in dataset ----
waPull=read_csv("waPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                     #select variables necessary
         `Footprints (1/cm2)`,
         `Whole Animal Adhesion (N/cm2)`) %>%
  rename(ind=Individual,                              #rename variables for ease of use
         set=`Animal Set`,
         fp=`Footprints (1/cm2)`,
         adh=`Whole Animal Adhesion (N/cm2)`) %>%
  mutate(set=fct_relevel(set,                         #order levels by time out of ocean
                         "WA2","CA8")) %>%
  drop_na()                                           #drop rows with NaNs

### 2.9.2 Statistical analysis ----
mod.fp=glm(data=waPull,
           adh~fp*set,
           family=Gamma(link="log"))
residualPlot(mod.fp)
mod.fp.lmO=summary(mod.fp)                                    #print model summary
mod.fp.aO=Anova(update(mod.fp,                                #model from above
                       contrasts=list(set=contr.sum)),        #update with contrasts for type III sums of squares
                type="III",                                   #request type III sums of squares
                test.statistic="F")
mod.fp.R2=1-(mod.fp.lmO$deviance/mod.fp.lmO$null.deviance)    #calculate the R-squared value
pred.fp=waPull %>%
  reframe(
    fp=seq(min(fp),
           max(fp),
           length.out=100),
    .by=set) %>%                                              #generate fine-resolution x across the range for each Animal Set
  mutate(fp.p= predict(mod.fp,
                       newdata=.,
                       level=0,                               #generate predictions at those x values
                       type="response"),
         set=fct_relevel(set,
                         "WA2","CA8"))                        #relevel Animal Set so it matches raw data

simple_slopes(mod.fp)                                         #run simple slopes test - the results for the simple slopes are those rows with sstest in column 1

table.fp=makeTable("adh","glm",                               #create the row(s) of the table
                   nlevels(waPull$ind),
                   c("fp","set","fp:set"),
                   mod.fp.aO,mod.fp.lmO,
                   mod.fp.R2)
### 2.9.3 Generate plot ----
#### a. generate annotations for plot ----
#p-value for interaction
p.inter=as.expression(bquote(p[fp:set]==.(round_p(mod.fp.aO["fp:set","Pr(>F)"]))))
#R-squared value
R2=as.expression(bquote(R^2 == .(round(mod.fp.R2,2))))
#Combine into plot annotation
plot.annot=c(p.inter,R2)

#### b. plot ----
plot.fp=plotGeneral(waPull,"fp","adh",pred.fp,"fp.p",
                    0.01,seq(2.5,1.7,length=2),plot.annot,
                    expression("Footprints (cm"^{-2}*")"),
                    expression("Whole Animal Adhesion (N cm"^{-1}*")"),
                    0,setClrs,setShapes,setLineTypes)+
  scale_y_continuous(trans='log10')+                                      #present x-axis on log scale because so many low strain rate tests.
  annotate("text",                                                        #add annotations for differences between Animal Sets
           x=c(2.5),
           y=c(0.7),
           label=c("‡"),
           family="Arial",
           size=8/.pt)

#### c. clean up environment ----
rm(waPull,mod.fp,mod.fp.lmO,
   mod.fp.aO,mod.fp.R2,
   pred.fp,plot.annot,
   p.inter,R2)

## 2.10 Whole Animal Adhesion (N/cm2) | Broken Tubefeet (1/cm2) ----
### 2.10.1 Load in dataset ----
waPull=read_csv("waPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,                    #select variables necessary
         `Broken Tube Feet (1/cm2)`,
         `Whole Animal Adhesion (N/cm2)`) %>%
  rename(ind=Individual,                             #rename variables for ease of use
         set=`Animal Set`,
         btf=`Broken Tube Feet (1/cm2)`,
         adh=`Whole Animal Adhesion (N/cm2)`) %>%
  mutate(set=fct_relevel(set,                        #order levels by time out of ocean
                         "WA2","CA8")) %>%
  drop_na()                                          #drop rows with NaNs

### 2.10.2 Statistical analysis ----
mod.btf=glm(data=waPull,
            adh~btf*set,
            family=Gamma(link="log"))
residualPlot(mod.btf)
mod.btf.lmO=summary(mod.btf)                                    #print model summary
mod.btf.aO=Anova(update(mod.btf,                                #model from above
                        contrasts=list(set=contr.sum)),         #update with contrasts for type III sums of squares
                 type="III",                                    #request type III sums of squares
                 test.statistic="F")
mod.btf.R2=1-(mod.btf.lmO$deviance/mod.btf.lmO$null.deviance)   #calculate the R-squared value
pred.btf=waPull %>%
  reframe(
    btf=seq(min(btf),
            max(btf),
            length.out=100),
    .by=set) %>%                                                #generate fine-resolution x across the range for each Animal Set
  mutate(btf.p= predict(mod.btf,
                        newdata=.,
                        level=0,                             #generate predictions at those x values
                        type="response"),
         set=fct_relevel(set,
                         "WA2","CA8"))                          #relevel Animal Set so it matches raw data

simple_slopes(mod.btf)                                         #run simple slopes test - the results for the simple slopes are those rows with sstest in column 1

table.btf=makeTable("adh","glm",                                #create the row(s) of the table
                    nlevels(waPull$ind),
                    c("btf","set","btf:set"),
                    mod.btf.aO,mod.btf.lmO,
                    mod.btf.R2)
### 2.10.3 Generate plot ----
#### a. generate annotations for plot ----
#p-value for interaction
p.inter=as.expression(bquote(p[btf:set]==.(round_p(mod.btf.aO["btf:set","Pr(>F)"]))))
#R-squared value
R2=as.expression(bquote(R^2 == .(round(mod.btf.R2,2))))
#Combine into plot annotation
plot.annot=c(p.inter,R2)

#### b. plot ----
plot.btf=plotGeneral(waPull,"btf","adh",pred.btf,"btf.p",
                     0.01,seq(6,4,length=2),plot.annot,
                     expression("Broken Tube Feet (cm"^{-2}*")"),
                     expression("Whole Animal Adhesion (N cm"^{-1}*")"),
                     0,setClrs,setShapes,setLineTypes)+
  scale_y_continuous(trans='log10')+                                      #present x-axis on log scale because so many low strain rate tests.
  annotate("text",                                                        #add annotations for differences between Animal Sets
           x=c(1.07,2.6),
           y=c(0.18,0.75),
           label=c("‡","†"),
           family="Arial",
           size=8/.pt)

#### c. clean up environment ----
rm(waPull,mod.btf,mod.btf.lmO,
   mod.btf.aO,mod.btf.R2,
   pred.btf,plot.annot,
   p.inter,R2)

# 3 FIGURES ----
## 3.1 Figure 2 | across strain rate ----
### 3.1.1 Get legend ----
legend <- get_legend2(plot.adh+
                        theme(legend.box.margin=margin(0, 0, 0, 0)))

### 3.1.2 Combine plots ----
plots.2=cowplot::plot_grid(plot.ten+theme(legend.position="none"),
                           plot.bF+theme(legend.position="none"),
                           plot.adh+theme(legend.position="none"),
                           plot.atch+theme(legend.position="none"),
                           labels=c("A","B","C","D"),label_size=12,           #labels for the panels
                           nrow=2,ncol=2,
                           hjust=-0.3,                                        #adjust horizontal alignment
                           vjust=0.8,                                         #adjust vertical alignment
                           align="v")                                         #align plots vertically.

plots.2=ggdraw(plots.2)+
  draw_plot(plot.bF.inset+theme(legend.position="none"),
            0.6,0.86,
            0.12,0.14)+
  draw_plot(plot.adh.inset+theme(legend.position="none"),
            0.093,0.36,
            0.15,0.14)+
  draw_plot(plot.atch.inset+theme(legend.position="none"),
            0.58,0.36,
            0.15,0.14)

### 3.1.3 Combine plots and legend ----
fig2=cowplot::plot_grid(legend,
                        plots.2,
                        labels=c("",""),label_size=12,
                        nrow=2,
                        ncol=1,
                        rel_heights=c(0.04,1))
### 3.1.4 Save figure ----
ggsave("figures/Fig2.pdf",
       width=180,height=145,
       units="mm",bg="white",device=cairo_pdf)
fig2
dev.off()
embed_fonts("figures/Fig2.pdf")

## 3.2 Figure 3 | across attachment points ----
### 3.2.1 Get legend ----
legend <- get_legend2(plot.fp+
                        theme(legend.box.margin=margin(0, 0, 0, 0)))

### 3.2.2 Combine plots ----
plots.3=cowplot::plot_grid(plot.fp+theme(legend.position="none"),
                           plot.btf+theme(legend.position="none"),
                           labels="AUTO",label_size=12,
                           nrow=2,
                           rel_heights=c(1,1),
                           hjust=-0.3,
                           vjust=1.2)

### 3.2.3 Combine plots and legend ----
fig3=cowplot::plot_grid(legend,
                        plots.3,
                        labels=c("",""),
                        nrow=2,
                        ncol=1,
                        rel_heights=c(0.05,1))

### 3.2.4 Save figure ----
ggsave("figures/Fig3.pdf",
       width=90/25.4,height=145/25.4,
       units="in",bg="white",device=cairo_pdf)
fig3
dev.off()
embed_fonts("figures/Fig3.pdf")

## 3.3 Figure S1 | material properties ----
### 3.3.1 Get legend ----
legend <- get_legend2(plot.xT+
                        theme(legend.box.margin=margin(0, 0, 0, 0)))

### 3.3.2 Combine plots ----
plots.S1=cowplot::plot_grid(plot.xT+theme(legend.position="none"),
                            plot.iSC+theme(legend.position="none"),
                            plot.fSC+theme(legend.position="none"),
                            plot.wF+theme(legend.position="none"),
                            labels=c("A","B","C","D"),label_size=12,
                            nrow=2,ncol=2,
                            hjust=-0.2,
                            vjust=0.4,
                            align="v")

plots.S1=ggdraw(plots.S1)+
  draw_plot(plot.xT.inset+theme(legend.position="none"),
            0.115,0.86,
            0.12,0.14)+
  draw_plot(plot.wF.inset+theme(legend.position="none"),
            0.62,0.36,
            0.12,0.14)

### 3.3.3 Combine plots and legend ----
figS1=cowplot::plot_grid(legend,
                         plots.S1,
                         labels=c("",""),
                         nrow=2,
                         ncol=1,
                         rel_heights=c(0.05,1))

### 3.3.4 Save figure ----
ggsave("figures/FigS1.pdf",
       width=180,height=145,
       units="mm",bg="white",device=cairo_pdf)
figS1
dev.off()
embed_fonts("figures/FigS1.pdf")

# 4 TABLES ----
## 4.1 Table 1 ----
#Get Tube Foot Disc Tenacity (MPa) summary statistics
tfD.descStats=read_csv("tfDPull.csv",col_types='ffnn') %>%
  rename(ind=Individual,                             #rename variables
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         ten=`Tube Foot Disc Tenacity (MPa)`) %>%
  mutate(set=fct_relevel(set,                        #order levels in time removed from ocean
                         "WA2")) %>% 
  group_by(set) %>%
  summarise(mean.ten=mean(ten,na.rm=TRUE),
            min.ten=min(ten,na.rm=TRUE),
            max.ten=max(ten,na.rm=TRUE))

#Get Tube Foot Stem Material Properties summary statistics
tfS.descStats=read_csv("tfSPull.csv",col_types='ff') %>%
  rename(ind=Individual,
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         bF=`Tube Foot Stem Breaking Force (N)`,
         xT=`Tube Foot Stem Extension (m)`,
         iSC=`Tube Foot Stem Initial Spring Constant (N/m)`,
         fSC=`Tube Foot Stem Final Spring Constant (N/m)`,
         wF=`Tube Foot Stem Work to Failure (J)`) %>%
  mutate(set=fct_relevel(set,
                         "WA2"),
         set=fct_expand(set,
                        "Bodega",
                        after=1)) %>%
  group_by(set) %>%
  summarise(mean.bF=mean(bF,na.rm=TRUE),
            min.bF=min(bF,na.rm=TRUE),
            max.bF=max(bF,na.rm=TRUE),
            mean.xT=mean(xT,na.rm=TRUE),
            min.xT=min(xT,na.rm=TRUE),
            max.xT=max(xT,na.rm=TRUE),
            mean.iSC=mean(iSC,na.rm=TRUE),
            min.iSC=min(iSC,na.rm=TRUE),
            max.iSC=max(iSC,na.rm=TRUE),
            mean.fSC=mean(fSC,na.rm=TRUE),
            min.fSC=min(fSC,na.rm=TRUE),
            max.fSC=max(fSC,na.rm=TRUE),
            mean.wF=mean(wF,na.rm=TRUE),
            min.wF=min(wF,na.rm=TRUE),
            max.wF=max(wF,na.rm=TRUE)) %>%
  add_row(set="Bodega",.before=2)

#Get Whole Animal variables summary statistics
wa.descStats=read_csv("waPull.csv",col_types='ff') %>%
  select(Individual,`Animal Set`,
         `Strain Rate (1/s)`,
         `Maximum Adhesive Force (N)`,`Whole Animal Adhesion (N/cm2)`,
         `Footprints (1/cm2)`,`Broken Tube Feet (1/cm2)`) %>%
  rename(ind=Individual,
         set=`Animal Set`,
         rate=`Strain Rate (1/s)`,
         mxf=`Maximum Adhesive Force (N)`,
         adh=`Whole Animal Adhesion (N/cm2)`,
         fp=`Footprints (1/cm2)`,
         btf=`Broken Tube Feet (1/cm2)`) %>%
  mutate(set=fct_relevel(set,
                         "WA2"),
         atch=fp+btf) %>%
  group_by(set) %>%
  summarise(mean.mxf=mean(mxf,na.rm=TRUE),
            min.mxf=min(mxf,na.rm=TRUE),
            max.mxf=max(mxf,na.rm=TRUE),
            mean.adh=mean(adh,na.rm=TRUE),
            min.adh=min(adh,na.rm=TRUE),
            max.adh=max(adh,na.rm=TRUE),
            mean.atch=mean(atch,na.rm=TRUE),
            min.atch=min(atch,na.rm=TRUE),
            max.atch=max(atch,na.rm=TRUE))

# Format the above statistics so they fit into the table.
ten.Formatted=paste0(sprintf('%.3f',tfD.descStats$mean.ten),
                     "[",
                     sprintf('%.3f',tfD.descStats$min.ten),",",
                     sprintf('%.3f',tfD.descStats$max.ten),
                     "]")
bF.Formatted=paste0(sprintf('%.2f',tfS.descStats$mean.bF),
                    "[",
                    sprintf('%.2f',tfS.descStats$min.bF),",",
                    sprintf('%.2f',tfS.descStats$max.bF),
                    "]")
xT.Formatted=paste0(sprintf('%.1f',tfS.descStats$mean.xT*100),
                    "[",
                    sprintf('%.1f',tfS.descStats$min.xT*100),",",
                    sprintf('%.1f',tfS.descStats$max.xT*100),
                    "]")
iSC.Formatted=paste0(sprintf('%.2f',tfS.descStats$mean.iSC),
                     "[",
                     sprintf('%.2f',tfS.descStats$min.iSC),",",
                     sprintf('%.2f',tfS.descStats$max.iSC),
                     "]")
fSC.Formatted=paste0(sprintf('%.1f',tfS.descStats$mean.fSC),
                     "[",
                     sprintf('%.1f',tfS.descStats$min.fSC),",",
                     sprintf('%.1f',tfS.descStats$max.fSC),
                     "]")
wF.Formatted=paste0(sprintf('%.3f',tfS.descStats$mean.wF),
                    "[",
                    sprintf('%.3f',tfS.descStats$min.wF),",",
                    sprintf('%.3f',tfS.descStats$max.wF),
                    "]")
mxf.Formatted=paste0(sprintf('%.1f',wa.descStats$mean.mxf),
                     "[",
                     sprintf('%.1f',wa.descStats$min.mxf),",",
                     sprintf('%.1f',wa.descStats$max.mxf),
                     "]")
adh.Formatted=paste0(sprintf('%.2f',wa.descStats$mean.adh),
                     "[",
                     sprintf('%.2f',wa.descStats$min.adh),",",
                     sprintf('%.2f',wa.descStats$max.adh),
                     "]")
atch.Formatted=paste0(sprintf('%.2f',wa.descStats$mean.atch),
                      "[",
                      sprintf('%.2f',wa.descStats$min.atch),",",
                      sprintf('%.2f',wa.descStats$max.atch),
                      "]")

# Combine into a table
table1=data.frame(tfD.descStats$set,
                  ten.Formatted,
                  bF.Formatted,xT.Formatted,
                  iSC.Formatted,fSC.Formatted,
                  wF.Formatted,
                  atch.Formatted,
                  mxf.Formatted,adh.Formatted) %>%
  rename(Set=tfD.descStats.set,
         `Tube Foot Disc Tenacity (MPa)`=ten.Formatted,
         `Tube Foot Stem Breaking Force (N)`=bF.Formatted,
         `Tube Foot Stem Extension (mm)`=xT.Formatted,
         `Tube Foot Stem Initial Spring Constant (N/m)`=iSC.Formatted,
         `Tube Foot Stem Final Spring Constant (N/m)`=fSC.Formatted,
         `Tube Foot Stem Work to Failure (J)`=wF.Formatted,
         `Number of Attachments (1/cm2)`=atch.Formatted,
         `Maximum Adhesive Force (N)`=mxf.Formatted,
         `Whole Animal Adhesion (N/cm2)`=adh.Formatted)

# Save table
write.csv(table1,"tables/table1.csv")

## 4.2 Table 2 ----
# Make table
table2=rbind(table.ten,
             table.bF,
             table.atch,
             table.adh,
             table.fp,
             table.btf)

# Save table
write.csv(table2,"tables/table2.csv")
## 4.3 Table S1 ----
# Make table
tableS1=rbind(table.xT,
              table.iSC,
              table.fSC,
              table.wF)

# Save table
write.csv(tableS1,"tables/tableS1.csv")
