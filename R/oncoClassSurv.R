#' oncoClassSurv
#'
#'Classify tumors into several molecular sub-types using "random forest" or
#'"support vector machine" according to RNA-Seq data.
#'Visualization for survival risk based on Cox regression.
#' @param exp.type A character for the normalized format of RNA-Seq expression matrix. The default
#' is "fpkm", and the other option is "tpm". Please use the consistent
#' "exp.type" in both train expression matrix and input expression matrix.
#' @param train.exp.path The path of the RNA-Seq expression profile of the training cohort.
#' The default path is 'system.file("extdata", paste0("train.tumor.exp.", exp.type,
#' ".txt"),package = "oncoClassSurv")' depending on the "exp.type".
#' This path can be any other customized path as per the user's demands.
#' @param train.clin.path The path of a table file containing clinical information. When performing the
#' classification task, the information on clusters in the training cohort is
#' essential for machine learning. When performing the prediction for survival
#' risk, survival data of the training cohort is essential for fitting the Cox
#' regression model.
#' @param train_cluster.feature.path This option can be customized. The path of a table file containing the marker
#' genes of clusters. When performing the classification task, the cluster-
#' specific marker genes is essential for precise training. Theoretically,
#' users can use all genes in the expression profile, but this is not recommended.
#' @param train_survival.feature.path Optional item. The path of a table file (.txt, csv, xlsx, xls, etc. or .rds)
#' including all covariables for training in the Cox regression model, when
#' prediction for survival risk needs to be performed. The default value is null.
#' @param input.exp.path The path of the RNA-Seq expression profile of the user's cohort. Please use
#' the consistent normalized format per the expression matrix in the training
#' cohort.
#' @param input.clin.path Optional item. The path of a table file (.txt, csv, xlsx, xls, etc. or .rds)
#' including additional clinical variables (e.g. tumor stage, age, gender, etc.),
#' when prediction for the survival risk of input data needs to be performed.
#' The default value is null.
#' @param miss_go.on Logical value. The function will automatically check whether missing genes
#' exist in the user's input expression matrix. If there are missing genes,
#' selecting "TRUE" means using the common genes (the training expression
#' matrix and the user's input expression matrix) and continuing to perform
#' the current task. This may cause imprecise prediction.
#' @param rm.batch.effect Logical value. Whether need to remove the batch effects between
#' the expression matrix of the training cohort and the expression matrix input by the user.
#' @param task Number value. Perform prediction for classifications when the task is 1.
#' Perform prediction for survival risk when the task is 2. Perform prediction for classifications and survival risk when the task is 3. The default value is 3.
#' @param plot.combatch Logical value. After removing the batch effects, visualize and return the
#' comparison results. Using the principal component analysis. The default value is FALSE.
#' @param print.combat.plots Logical value. After removing the batch effects, visualize the comparison
#' results in the current graph panel. This item works when the item of "plot.combatch" is TRUE.
#' @param cluster.method A character indicating the algorithm used in the prediction for classification.
#' Optional values include "RF" (random forest) and "SVM" (support vector machine).
#' @param nodesize A parameter in the random forest. Details can be seen in
#' the \code{randomForest::\link{randomForest}}.
#' @param ntree A parameter in the random forest. Details can be seen in the
#' \code{randomForest::\link{randomForest}}.
#' @param mtry A parameter in the random forest. Details can be seen in the
#' \code{randomForest::\link{randomForest}}.
#' @param importance A parameter in the random forest. Details can be seen in the
#' \code{randomForest::\link{randomForest}}.
#' @param kernel A parameter in the support vector machine. Details can be seen in the
#' \code{e1071::\link{svm}}.
#' @param cost A parameter in the support vector machine. Details can be seen in the
#' \code{e1071::\link{svm}}.
#' @param surv.t.custom When predicting the survival risk of patients, users can specify
#' concerned time points. The default value is null, which means predicting for
#' all time points.
#' @param time A character indicating the categories of survival time, such as overall
#' survival (OS), recurrence-free survival (RFS), etc. The default value is "OS".
#' @param event A character or number indicating the categories of survival event, such as
#' death or recurrence customized by users. The default value is 1.
#' @param plot.surv.curve A logical value. Whether return curves of survival risk over time
#' for individual patients. The default value is TRUE.
#' @param plot.samples A number. When there are too many patients, plotting multiple survival curves
#' on one graph panel will appear too crowded. Therefore, it is recommended to
#' select only the patients concerned by the user (default is numbered 1, 2, 3,
#' 4, 5) to show their survival curves.
#' @param survtime.unit A character. The unit of the survival time, such as "days", "months", "year",
#' etc. as per the user's demands.
#' @param survcurve.break.x.by A number. The breakpoints of the time axis for the survival curves.
#' The default value is 12. Details can be seen in the
#' \code{survminer::\link{ggsurvplot}}.
#' @param print.survplot A logical value. Whether visualize curves of survival risk over time for
#' individual patients in the current graph panel. This item works when the item
#' of "plot.surv.curve" is TRUE. The default value is FALSE.
#' @param show.message Whether to display redundant detection information. The default value is TRUE.
#' @author Yang Li
#'
#' @return A list including classifications or survival risk results or both of them. The returned results depend on the task performed.
#' @export
#' @importFrom data.table fread
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom tibble column_to_rownames
#' @importFrom stats prcomp
#' @importFrom stats as.formula
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @import e1071
#' @import ggplot2
#' @import ggfortify
#' @import patchwork
#' @import survival
#' @import survminer
#' @import BiocManager
#' @import limma
#' @import sva
#' @examples
#' \dontrun{
#' #perform task 1 (classify by SVM)
#' results<-oncoClassSurv(input.exp.path = system.file("extdata",
#' "icgc.tumor.exp.fpkm.txt",package = "oncoClassSurv"),
#' miss_go.on=T,task=1,rm.batch.effect=TRUE,
#' plot.combatch=TRUE,print.combat.plots=TRUE,cluster.method="SVM",
#' show.message=FALSE)
#' }
#'
#' \dontrun{
#' #perform task 1 (classify by RF)
#' results<-oncoClassSurv(input.exp.path = system.file("extdata",
#' "icgc.tumor.exp.fpkm.txt",package = "oncoClassSurv"),
#' task=1,rm.batch.effect=TRUE,
#' plot.combatch=TRUE,print.combat.plots=TRUE,cluster.method="RF",
#' show.message=FALSE)
#' }
#'
#' \dontrun{
#' #perform task 2
#' results<-oncoClassSurv(train_survival.feature.path=system.file("extdata",
#'  "train_survival.features.rds",package = "oncoClassSurv"),
#'  input.exp.path = system.file("extdata", "icgc.tumor.exp.fpkm.txt",
#'  package = "oncoClassSurv"),input.clin.path = system.file("extdata",
#'  "input_clinsurv.txt",package = "oncoClassSurv"),
#'  task=2,rm.batch.effect=TRUE,plot.combatch=TRUE,
#'  print.combat.plots=TRUE,surv.t.custom=NULL,
#'  plot.surv.curve=TRUE,survcurve.break.x.by = 12,
#'  print.survplot = TRUE,plot.samples=c(1:5),
#'  show.message=FALSE)
#' }
#'
#' \dontrun{
#' #perform task 3
#' results<-oncoClassSurv(train_survival.feature.path=system.file("extdata",
#' "train_survival.features.rds",package = "oncoClassSurv"),
#' input.exp.path = system.file("extdata", "icgc.tumor.exp.fpkm.txt",
#' package = "oncoClassSurv"),input.clin.path = system.file("extdata",
#' "input_clinsurv.txt",package = "oncoClassSurv"),task=3,
#' rm.batch.effect=TRUE,plot.combatch=TRUE,print.combat.plots=TRUE,
#' cluster.method="SVM",surv.t.custom=NULL,plot.surv.curve=TRUE,
#' survcurve.break.x.by = 12,print.survplot = TRUE,plot.samples=c(1:5),
#' show.message=FALSE)
#' }
##' @seealso
##' \code{sva::\link{ComBat}};
##' \code{randomForest::\link{randomForest}};
##' \code{e1071::\link{svm}};
##' \code{survminer::\link{ggsurvplot}}

oncoClassSurv<-function(exp.type="fpkm",
                        train.exp.path=system.file("extdata", paste0("train.tumor.exp.",exp.type,".txt"),
                                                   package = "oncoClassSurv"),
                        train.clin.path=system.file("extdata", "train.cluster.surv.rds",
                                                    package = "oncoClassSurv"),
                        train_cluster.feature.path=system.file("extdata", "train_cluster.features.rds",
                                                               package = "oncoClassSurv"),
                        train_survival.feature.path=NULL,
                        input.exp.path,input.clin.path=NULL,
                        miss_go.on=TRUE,rm.batch.effect=TRUE,
                        task=3,plot.combatch=FALSE,
                        print.combat.plots=FALSE,
                        cluster.method="RF",nodesize = 3,
                        ntree = 3000,mtry = 5,importance = T,
                        kernel="radial",cost=4,
                        surv.t.custom=NULL,time="OS",event=1,
                        plot.surv.curve=TRUE,plot.samples=c(1:5),
                        survtime.unit="months",
                        survcurve.break.x.by=12,print.survplot=FALSE,show.message=TRUE){

  #收集分析结果
  return.list<-list()

  #1.上传数据：表达谱####
  train.tumor.exp<-data.table::fread(file = train.exp.path,data.table = F,showProgress = T)%>%
    tibble::column_to_rownames(var = "Features")
  input.tumor.exp<-data.table::fread(file = input.exp.path,data.table = F,showProgress = T)%>%
    tibble::column_to_rownames(var = "Features")

  rownames(train.tumor.exp)<-gsub(rownames(train.tumor.exp),pattern="-",replacement="_")
  rownames(input.tumor.exp)<-gsub(rownames(input.tumor.exp),pattern="-",replacement="_")

  #2.统计基因与样本数目####
  if(show.message){
    message(paste0("Input: ",dim(input.tumor.exp)[1]," Features",", ",dim(input.tumor.exp)[2]," Samples.\n"))
  }

  #3.缺失值检查check.NA####
  check.NA<-function (data) {
    data_class <- class(data)
    #message(paste0("Data class of your input: ",data_class,"."))
    #require(data.table)
    is_data_table <- data.table::is.data.table(data)
    if (!is_data_table) {
      data <- data.table::data.table(data)
      missing_summary <-
        data.table::data.table(Features = names(data),
                               Number_missing = base::sapply(data, function(x) {sum(is.na(x))}),
                               Locate_missing = base::sapply(data, function(x) {which(is.na(x))}))
    }
    return(missing_summary)
  }
  check.NA.results<-check.NA(data = input.tumor.exp)
  if(check.NA.results%>%as.data.frame()%>%.$Number_missing%>%sum()>0){stop("Stop for missing data.Please check the NA.results in the missing analysis.")}
  return.list[["check.NA.results"]]<-check.NA.results
  #4.基因名称转化####
  #Features名称标准化:"-"转化为:"_"。
  if(task==1){
    if(grepl(x=train_cluster.feature.path,pattern = "\\.rds$")){
      cluster_markergenes<-readRDS(file = train_cluster.feature.path)
    }else{
      cluster_markergenes<-data.table::fread(train_cluster.feature.path,data.table = F)[,1]
    }
    cluster_markergenes<-gsub(cluster_markergenes,pattern="-",replacement="_")

    #根据要执行的操作，check是否包含classifier模型所需要的基因？
    None.feature<-cluster_markergenes[!cluster_markergenes%in%rownames(input.tumor.exp)]
    if(length(None.feature)>0){
      if(miss_go.on){
        if(show.message){
          message(paste0("\nPlease check your input features!\nNotice: Default option is to continue when missing features exist in the input data, which may cause reduced accuracy.",
                         "\nNumber of Marker features in the train data: ",length(cluster_markergenes),
                         ", but ",nrow(input.tumor.exp)," in your input data."))
          message(paste0("\nMissing features: ",paste0(None.feature,collapse = ",")))
        }
        #最后，选择cluster_markergenes与input的共同基因用于后续分型
        cluster_markergenes<-base::intersect(cluster_markergenes,rownames(input.tumor.exp))
        train.tumor.exp<-train.tumor.exp[cluster_markergenes,]
        input.tumor.exp<-input.tumor.exp[cluster_markergenes,]
      }else{
        stop("Need complete input file including all essential features.")
      }
    }
  }
  if(task==2){
    if(grepl(x=train_survival.feature.path,pattern = "\\.rds$")){
      prog.signif.features<-readRDS(file = train_survival.feature.path)
    }else{
      prog.signif.features<-data.table::fread(train_survival.feature.path,data.table = F)[,1]
    }
    prog.signif.features<-gsub(prog.signif.features,pattern="-",replacement="_")
    #保留表达谱文件中的共同基因
    co_genes<-base::intersect(rownames(train.tumor.exp),rownames(input.tumor.exp))
    train.tumor.exp<-train.tumor.exp[co_genes,]
    input.tumor.exp<-input.tumor.exp[co_genes,]
  }
  if(task==3){
    if(grepl(x=train_cluster.feature.path,pattern = "\\.rds$")){
      cluster_markergenes<-readRDS(file = train_cluster.feature.path)
    }else{
      cluster_markergenes<-data.table::fread(train_cluster.feature.path,data.table = F)[,1]
    }

    if(grepl(x=train_survival.feature.path,pattern = "\\.rds$")){
      prog.signif.features<-readRDS(file = train_survival.feature.path)
    }else{
      prog.signif.features<-data.table::fread(train_survival.feature.path,data.table = F)[,1]
    }

    cluster_markergenes<-gsub(cluster_markergenes,pattern="-",replacement="_")
    prog.signif.features<-gsub(prog.signif.features,pattern="-",replacement="_")

    #根据要执行的操作，check是否包含classifier模型所需要的基因？
    #此处仅检查cluster_markergenes是否存在缺失基因。由于预后prog.signif.features通常所涉及的基因数量有限，且可能存在临床变量，不便自动检查。因此，应由用户自行评估是否存在缺失基因。
    None.feature<-cluster_markergenes[!cluster_markergenes%in%rownames(input.tumor.exp)]
    if(length(None.feature)>0){
      if(miss_go.on){
        message(paste0("\nPlease check your input features!\nNotice: Default option is to continue when missing features exist in the input data, which may cause reduced accuracy.",
                       "\nNumber of Marker features in the train data: ",length(cluster_markergenes),
                       ", but ",nrow(input.tumor.exp)," in your input data."))
        message(paste0("\nMissing features: ",paste0(None.feature,collapse = ",")))
        #最后，选择cluster_markergenes与input的共同基因用于后续分型
        cluster_markergenes<-base::intersect(cluster_markergenes,rownames(input.tumor.exp))
        train.tumor.exp<-train.tumor.exp[cluster_markergenes,]
        input.tumor.exp<-input.tumor.exp[cluster_markergenes,]
      }else{
        stop("Need complete input file including all essential features.")
      }
    }
  }

  #5.标准化####
  #5.1使用limma去重基因，平均，标准化
  input.tumor.exp=as.matrix(input.tumor.exp)
  input.tumor.exp=limma::avereps(input.tumor.exp)
  #limma包标准化函数normalizeBetweenArrays()
  input.tumor.exp <- limma::normalizeBetweenArrays(input.tumor.exp)
  #5.2是否已经进行log2(FPKM+1)？若无，需要对表达谱进行log2(FPKM+1)转化
  #如果max(x)>50,警告，需要FPKM数据，请user检查数据！！
  #log2(fpkm+1)
  train.tumor.exp<-log2(train.tumor.exp+1)

  if(min(input.tumor.exp)>=0){
    if(max(input.tumor.exp)>50){
      message("Input data is performing log2(expression+1)...")
      input.tumor.exp<-log2(input.tumor.exp+1)
      message("log2(expression+1) finished.")
    }
  }else{
    stop("Check your input data for negative expression values.")
  }

  #6.选择是否进行combat()####
  if(rm.batch.effect){
    exp.join<-cbind(train.tumor.exp,input.tumor.exp[rownames(train.tumor.exp),])
    batch<-base::rep(c("Train","Input"),c(dim(train.tumor.exp)[2],dim(input.tumor.exp)[2]))
    batch<-factor(batch,levels = c("Train","Input"))

    #library(sva)
    batchremove_combat <- sva::ComBat(dat = as.matrix(exp.join), batch = batch)

    if(plot.combatch){
      #library(ggfortify)
      exp.join.dt <- as.data.frame(t(exp.join))
      exp.join.dt$batch <- batch
      exp.join.pca<-ggplot2::autoplot(stats::prcomp(exp.join.dt[,1:(ncol(exp.join.dt)-1)] ),
                                      data=exp.join.dt,colour = 'batch',
                                      frame.type = 'norm')+
        ggplot2::theme_bw()+
        ggplot2::theme(line =ggplot2::element_line(linewidth = 0.2),
                       legend.position = "top",
                       legend.title =ggplot2::element_text(size = 6,color = "black") ,
                       legend.text =ggplot2::element_text(size = 6),
                       legend.key.size =ggplot2::unit(3,units = "mm"),
                       axis.text = ggplot2::element_text(size = 6,color = "black"),
                       axis.text.x =ggplot2::element_blank(),
                       axis.title =ggplot2::element_text(size = 8,color = "black"))+
        ggplot2::labs(title ="Before ComBat")

      combat <- as.data.frame(t(batchremove_combat))
      combat$batch <- batch
      combat.pca<-ggplot2::autoplot(stats::prcomp(combat[,1:(ncol(combat)-1)] ),
                                    data=combat,colour = 'batch',
                                    frame.type = 'norm')+
        ggplot2::theme_bw()+
        ggplot2::theme(line =ggplot2::element_line(linewidth = 0.2),
                       legend.position = "top",
                       legend.title =ggplot2::element_text(size = 6,color = "black") ,
                       legend.text =ggplot2::element_text(size = 6),
                       legend.key.size =ggplot2::unit(3,units = "mm"),
                       axis.text = ggplot2::element_text(size = 6,color = "black"),
                       axis.text.x =ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(size = 8,color = "black"))+
        ggplot2::labs(title ="After ComBat")

      original_combat.p<-exp.join.pca+combat.pca
      if(print.combat.plots){
        print(original_combat.p)
      }
      return.list[["original_combat.plots"]]<-original_combat.p
    }
  }


  #7.根据第6步的选项，采用原模型，或使用combat后的数据拟合模型（分型模型(支持向量机/随机森林)/预后模型）
  #根据task选项，进行预测（分型；随时间变化的生存概率）

  #合并clinsurv数据####
  if(grepl(x=train.clin.path,pattern = "\\.rds$")){
    train.tumor.clin.surv<-readRDS(file = train.clin.path)
  }else{
    train.tumor.clin.surv<-data.table::fread(train.clin.path,data.table = F)
    rownames(train.tumor.clin.surv)<-train.tumor.clin.surv[,1]
    train.tumor.clin.surv[,1]=NULL
  }

  if(task%in%c(2,3)&!is.null(input.clin.path)){
    input.clinsurv<-data.table::fread(file = input.clin.path,data.table = F)
    rownames(input.clinsurv)<-input.clinsurv[,"sample_name"]
  }


  if(rm.batch.effect){
    train.tumor.exp<-batchremove_combat[,1:dim(train.tumor.exp)[2]]
    train.tumor.exp<-train.tumor.exp%>%t()%>%as.data.frame()
    train.tumor.dt<-cbind(train.tumor.clin.surv,train.tumor.exp[rownames(train.tumor.clin.surv),])

    input.tumor.exp=batchremove_combat[,(dim(train.tumor.exp)[1]+1):dim(batchremove_combat)[2]]
    input.tumor.exp<-input.tumor.exp%>%t()%>%as.data.frame()

    if(task%in%c(2,3)&!is.null(input.clin.path)){
      input.tumor.dt<-cbind(input.clinsurv[rownames(input.tumor.exp),],input.tumor.exp)
    }else{if(task%in%c(2,3)){
      input.tumor.dt<-input.tumor.exp
    }
    }
  }else{
    train.tumor.exp<-train.tumor.exp%>%t()%>%as.data.frame()
    train.tumor.dt<-cbind(train.tumor.clin.surv,train.tumor.exp[rownames(train.tumor.clin.surv),])

    input.tumor.exp=input.tumor.exp
    input.tumor.exp<-input.tumor.exp%>%t()%>%as.data.frame()
    if(task%in%c(2,3)&!is.null(input.clin.path)){
      input.tumor.dt<-cbind(input.clinsurv[rownames(input.tumor.exp),],input.tumor.exp)
    }else{if(task%in%c(2,3)){
      input.tumor.dt<-input.tumor.exp
    }
    }
  }



  #choose task and perform that#
  if(task==1){
    #构建公式
    cluster.fmla <- stats::as.formula(
      paste0("Cluster ~ ",paste0("`", cluster_markergenes ,"`", collapse = " + ")))

    if(cluster.method=="RF"){
      #library(randomForest)
      # 构建RF模型
      set.seed(1234) # 保证结果的可重复性
      RF.fit<- randomForest::randomForest(
        cluster.fmla,
        nodesize = nodesize,
        data = train.tumor.dt,
        ntree = ntree, # 决策树棵数
        mtry = mtry, # 每个节点可供选择的变量数目
        importance = importance # 输出变量重要性
      )
      #预测
      rf.cluster.pred <- stats::predict(RF.fit,newdata = input.tumor.exp,type = "class")
      rf.cluster.pred.dt<-data.frame(ID=names(rf.cluster.pred),rf.cluster.pred)
      return.list[["rf.cluster"]]<-list(cluster.method=cluster.method,
                                        rf.cluster.pred=rf.cluster.pred.dt)
    }else{
      if(cluster.method=="SVM"){
        #library(e1071)
        #构建SVM模型
        SVM.fit.r <- e1071::svm(cluster.fmla, data=train.tumor.dt,kernel=kernel,cost=cost)
        #预测
        svm.cluster.pred <- stats::predict(SVM.fit.r,newdata = input.tumor.exp,type = "class")
        svm.cluster.pred.dt<-data.frame(ID=names(svm.cluster.pred),svm.cluster.pred)
        return.list[["svm.cluster"]]<-list(cluster.method=cluster.method,
                                           svm.cluster.pred=svm.cluster.pred.dt)
      }
    }
  }
  if(task==2){
    #构建surv公式
    surv.fmla <- stats::as.formula(
      paste0("Surv(",time,",","status==", event,") ~",paste0("`",prog.signif.features ,"`", collapse = " + ")))
    #构建coxphfit模型
    train.coxfit<-survival::coxph(surv.fmla,data = train.tumor.dt)

    #预测生存曲线
    input.surv.curve<-survival::survfit(train.coxfit,newdata = input.tumor.dt)
    #根据生存曲线计算生存概率
    input.surv.probablity=data.frame(Time=summary(input.surv.curve)$time,
                                     SurvivalProbablity=summary(input.surv.curve)$surv)

    return.list[["surv.probablity"]]<-input.surv.probablity

    #自定义时间点####
    if(!is.null(surv.t.custom)){
      surv.probablity.custom<-summary(input.surv.curve,times = surv.t.custom)$surv
      return.list[["surv.probablity.custom"]]<-list(surv.probablity.custom)
    }

    #是否绘制生存曲线####
    if(plot.surv.curve==T){
      #library(survminer)
      if(length(plot.samples)<=25){
        gg.ncol=5;
        gg.size=5;
        gg.key=1;
        gg.legend.position="top"
      }else{
        gg.ncol=ceiling(length(plot.samples)/5);
        gg.size=5*5/gg.ncol;
        gg.key=1*5/gg.ncol;
        gg.legend.position="top"
      }
      if(length(plot.samples)>50){
        gg.legend.position="none"
      }

      if(!is.null(survtime.unit)){survtime.axis=
        paste0("Time ","(",survtime.unit,")")}else{
          survtime.axis="Time"
        }
      input.ggsurv.curve<-survminer::ggsurvplot(input.surv.curve[plot.samples],data =input.tumor.dt,
                                                fun = "pct",conf.int = F,surv.median.line = "hv",
                                                legend.title="Sample",legend=gg.legend.position,
                                                break.x.by = survcurve.break.x.by,
                                                xlab=survtime.axis,
                                                ggtheme = survminer::theme_survminer()+
                                                  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0),
                                                                 panel.border = ggplot2::element_rect(color="black",fill ="transparent"),
                                                                 legend.title = ggplot2::element_text(size =gg.size+1),
                                                                 legend.text = ggplot2::element_text(size=gg.size)))+
        ggplot2::guides(col=ggplot2::guide_legend(ncol = gg.ncol,
                                                  keywidth = ggplot2::unit(gg.key,"mm"),
                                                  keyheight = ggplot2::unit(gg.key,"mm")))
      if(print.survplot==T){
        print(input.ggsurv.curve)
      }
      return.list[["ggsurv.curve"]]<-list( plot.samples=plot.samples,
                                           ggsurv.curve=input.ggsurv.curve)}}
  if(task==3){
    #cluster分型预测
    #构建公式
    cluster.fmla <- stats::as.formula(
      paste0("Cluster ~ ",paste0("`", cluster_markergenes ,"`", collapse = " + ")))

    if(cluster.method=="RF"){
      #library(randomForest)
      # 构建RF模型
      set.seed(1234) # 保证结果的可重复性
      RF.fit<- randomForest::randomForest(
        cluster.fmla,
        nodesize = nodesize,
        data = train.tumor.dt,
        ntree = ntree, # 决策树棵数
        mtry = mtry, # 每个节点可供选择的变量数目
        importance = importance # 输出变量重要性
      )
      #预测
      rf.cluster.pred <- stats::predict(RF.fit,newdata = input.tumor.exp,type = "class")
      rf.cluster.pred.dt<-data.frame(ID=names(rf.cluster.pred),rf.cluster.pred)
      return.list[["rf.cluster"]]<-list(cluster.method=cluster.method,
                                        rf.cluster.pred=rf.cluster.pred.dt)
    }else{
      if(cluster.method=="SVM"){
        #library(e1071)
        #构建SVM模型
        SVM.fit.r <- e1071::svm(cluster.fmla, data=train.tumor.dt,kernel=kernel,cost=cost)
        #预测
        svm.cluster.pred <- stats::predict(SVM.fit.r,newdata = input.tumor.exp,type = "class")
        svm.cluster.pred.dt<-data.frame(ID=names(svm.cluster.pred),svm.cluster.pred)
        return.list[["svm.cluster"]]<-list(cluster.method=cluster.method,
                                           svm.cluster.pred=svm.cluster.pred.dt)
      }
    }
    #生存预测
    #构建surv公式
    surv.fmla <- stats::as.formula(
      paste0("Surv(",time,",","status==", event,") ~",paste0("`",prog.signif.features ,"`", collapse = " + ")))
    #构建coxphfit模型
    train.coxfit<-survival::coxph(surv.fmla,data = train.tumor.dt)
    #预测生存曲线
    input.surv.curve<-survival::survfit(train.coxfit,newdata = input.tumor.dt)
    #根据生存曲线计算生存概率
    input.surv.probablity=data.frame(Time=summary(input.surv.curve)$time,
                                     SurvivalProbablity=summary(input.surv.curve)$surv)
    return.list[["surv.probablity"]]<-input.surv.probablity

    #自定义时间点####
    if(!is.null(surv.t.custom)){
      surv.probablity.custom<-summary(input.surv.curve,times = surv.t.custom)$surv
      return.list[["surv.probablity.custom"]]<-list(surv.probablity.custom)
    }

    #是否绘制生存曲线####
    if(plot.surv.curve==T){
      #library(survminer)
      if(length(plot.samples)<=25){
        gg.ncol=5;
        gg.size=5;
        gg.key=1;
        gg.legend.position="top"
      }else{
        gg.ncol=ceiling(length(plot.samples)/5);
        gg.size=5*5/gg.ncol;
        gg.key=1*5/gg.ncol;
        gg.legend.position="top"
      }
      if(length(plot.samples)>50){
        gg.legend.position="none"
      }

      if(!is.null(survtime.unit)){survtime.axis=
        paste0("Time ","(",survtime.unit,")")}else{
          survtime.axis="Time"
        }
      input.ggsurv.curve<-survminer::ggsurvplot(input.surv.curve[plot.samples],data =input.tumor.dt,
                                                fun = "pct",conf.int = F,surv.median.line = "hv",
                                                break.x.by = survcurve.break.x.by,
                                                legend.title="Sample",legend=gg.legend.position,
                                                xlab=survtime.axis,
                                                ggtheme = survminer::theme_survminer()+
                                                  ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0),
                                                                 panel.border = ggplot2::element_rect(color="black",fill ="transparent"),
                                                                 legend.title = ggplot2::element_text(size =gg.size+1),
                                                                 legend.text = ggplot2::element_text(size=gg.size)))+
        ggplot2::guides(col=ggplot2::guide_legend(ncol = gg.ncol,
                                                  keywidth = ggplot2::unit(gg.key,"mm"),
                                                  keyheight = ggplot2::unit(gg.key,"mm")))
      if(print.survplot==T){
        print(input.ggsurv.curve)
      }
      return.list[["ggsurv.curve"]]<-list( plot.samples=plot.samples,
                                           ggsurv.curve=input.ggsurv.curve)
    }
  }

  #8.生成预测结果，输出表格
  return(return.list)
}
