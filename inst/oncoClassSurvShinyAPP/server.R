suppressMessages(library(shiny))
suppressMessages(library(oncoClassSurv))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(ggplotify))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

#设置文件上传大小限制
options(shiny.maxRequestSize=1024^6)

mycol<-c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#E41A1C",
         "#4DAF4A","#F0027F","#FDB462","#BC806D","#377EB8","#FF00FF",
         "#A65628","#4253ff","#984EA3","#386c80","#F781BF","#c7475b",
         "#D8D155","#64495D","#B3DE69","#BEAED4","#80B1D3")

shinyServer(function(input, output,session) {
  output$value.gene <- renderPrint({ input$text.input.gene.explore })#打印输入文本

  #reactive Data####
  #0. exp.type Expression normalization method: fpkm/tpm （Optional. only for default.）
  exp.type.active<-eventReactive(input$run,{input$exp.type})

  #1.file.train.exp
  train.exp.path1<-eventReactive(input$run,{input$file.train.exp$datapath})

  #2.file.input.exp
  input.exp.path1<-eventReactive(input$run,{input$file.input.exp$datapath})

  #3.file.train.clin
  train.clin.path1<-eventReactive(input$run,{input$file.train.clin$datapath})

  #4.file.input.clin
  input.clin.path1<-eventReactive(input$run,{input$file.input.clin$datapath})

  #5.file.train.marker
  train_cluster.feature.path1<-eventReactive(input$run,{input$file.train.marker$datapath})

  #6.file.train.Prog.features
  train_survival.feature.path1<-eventReactive(input$run,{input$file.train.Prog.features$datapath})

  #reactive Parameter####
  #1. s1.task
  s1.task.active<-eventReactive(input$run,{input$s1.task})

  #2. s2.rm.batch.effect
  s2.rm.batch.effect.active<-eventReactive(input$run,{input$s2.rm.batch.effect})

  #3. s3.plot.combatch
  s3.plot.combatch.active<-eventReactive(input$run,{input$s3.plot.combatch})

  #4. s4.cluster.method
  s4.cluster.method.active<-eventReactive(input$run,{input$s4.cluster.method})

  #5. s5.kernel
  s5.kernel.active<-eventReactive(input$run,{input$s5.kernel})

  #6. s6.plot.surv.curve
  s6.plot.surv.curve.active<-eventReactive(input$run,{input$s6.plot.surv.curve})

  #7. ntree
  ntree.active<-eventReactive(input$run,{input$ntree})

  #8. nodesize
  nodesize.active<-eventReactive(input$run,{input$nodesize})

  #9. mtry
  mtry.active<-eventReactive(input$run,{input$mtry})

  #10. cost
  cost.active<-eventReactive(input$run,{input$cost})

  #11. survcurve.break.x.by
  survcurve.break.x.by.active<-eventReactive(input$run,{input$survcurve.break.x.by})

  #12. file.PartSamples.Pred Whether to select parts of samples to perform prediction analysis
  file.PartSamples.Pred.path1<-eventReactive(input$run,{input$file.PartSamples.Pred$datapath})

  #13. 根据update.surv, 刷新plot.samples 获取选择器输入框的值
  #plot.samples.active=NULL
  plot.samples.active <- eventReactive(input$run,{input$plot.samples})

  #14. legend.size
  leg.size.active <- eventReactive(input$run,{input$leg.size})

  #15. sel.geneset #勾选基因集
  sel.geneset.active <- eventReactive(input$run,{input$sel.geneset})

  #16. selgene.uplod #fileinput基因集
  selgene.uplod.active <- eventReactive(input$run,{
    data.table::fread(input$selgene.uplod$datapath,
                      data.table = F,header = F)[,1]})

  #17. text.input.gene.explore #输入text兴趣基因
  text.input.gene.explore.active <- eventReactive(input$run,{
    text.genes<-strsplit(x=input$text.input.gene.explore,split = ",|，|;| ")[[1]]
    text.genes<-text.genes[which(text.genes!="")]
    text.genes
    })

  #18. select.hp.gene #选择哪个输入组的基因呢？
  select.hp.gene.active <- eventReactive(input$run,{
    input$select.hp.gene
  })

  #19. survtime.unit #输入生存曲线时间单位
  survtime.unit.active <- eventReactive(input$run,{input$survtime.unit})

  #20. hp.sample #热图是否选择sample进行绘制，是否显示sample名称
  hp.sample.active<-eventReactive(input$run,{input$hp.sample})

  #21. limmaNorm.train Whether to perform limma::normalizeBetweenArrays()in the train cohort.
  limmaNorm.train.active<-eventReactive(input$run,{input$limmaNorm.train})

  #22. limmaNorm.input Whether to perform limma::normalizeBetweenArrays()in the input cohort.
  limmaNorm.input.active<-eventReactive(input$run,{input$limmaNorm.input})

  #23. RandomSamTrain Whether to use random samples to fit training models in the training cohort.
  RandomSamTrain.active<-eventReactive(input$run,{input$RandomSamTrain})

  #24. TrainSeed If RandomSamTrain is TRUE, use a random seed to control the results of random sampling.
  TrainSeed.active<-eventReactive(input$run,{input$TrainSeed})

  #25. SamplingProb If RandomSamTrain is TRUE, provide a sampling probability. From 0.01 to 1
  SamplingProb.active<-eventReactive(input$run,{input$SamplingProb})

  #save reactive values: Use reactiveValues to store the reactive value
  actvalue <- reactiveValues(path = NULL)

  runornot<-eventReactive(input$run,{TRUE})
  showNotification("Ready to run!",
                   action = a(href = "javascript:location.reload();",
                              "Reload"),id="ready")

  observeEvent(input$run,{

    if(runornot()){message("Running...")}

    showPageSpinner()


    #Get default Data path#####
    #0. exp.type Expression normalization method: fpkm/tpm
    if(!is.null(exp.type.active())){
      actvalue$exp.type<-ifelse(exp.type.active()=="1","fpkm","tpm")
    }else{
      actvalue$exp.type<-"fpkm"
    }

    #1.file.train.exp
    if(!is.null(input$file.train.exp)){
      actvalue$train.exp.path<-train.exp.path1()
      }else{
      actvalue$train.exp.path<-system.file("extdata", paste0("train.tumor.exp.",actvalue$exp.type,".txt"),
                                           package = "oncoClassSurv")
      }


    #2.file.input.exp
    if(!is.null(input$file.input.exp)){
      actvalue$input.exp.path<-input.exp.path1()
      }else{
      actvalue$input.exp.path<-system.file("extdata", paste0("icgc.tumor.exp.",actvalue$exp.type,".txt"),
                                           package = "oncoClassSurv")
      }

    #3.file.train.clin
    if(!is.null(input$file.train.clin)){
      actvalue$train.clin.path<-train.clin.path1()
      }else{
      actvalue$train.clin.path<-system.file("extdata", "train.cluster.surv.rds",
                                            package = "oncoClassSurv")
      }

    #4.file.input.clin
    if(!is.null(input$file.input.clin)){
      actvalue$input.clin.path<-input.clin.path1()
      }else{
      actvalue$input.clin.path<-system.file("extdata", "input_clinsurv.txt",
                                            package = "oncoClassSurv")
      }

    #5.file.train.marker
    if(!is.null(input$file.train.marker)){
      actvalue$train_cluster.feature.path<-train_cluster.feature.path1()
      }else{
      actvalue$train_cluster.feature.path<-system.file("extdata", "train_cluster.features.rds",
                                                       package = "oncoClassSurv")
      }

    #6.file.train.Prog.features
    if(!is.null(input$file.train.Prog.features)){
      actvalue$train_survival.feature.path<-train_survival.feature.path1()
      }else{
      actvalue$train_survival.feature.path<-system.file("extdata", "train_survival.features.rds",
                                                        package = "oncoClassSurv")
      }


    #Get default Parameter####
      #1. s1.task
      if(!is.null(s1.task.active())){
        actvalue$task<-as.numeric(s1.task.active())
      }else{
        actvalue$task<-1
      }


      #2. s2.rm.batch.effect
      if(!is.null(s2.rm.batch.effect.active())){
        actvalue$rm.batch.effect<-s2.rm.batch.effect.active()=="1"
      }else{
        actvalue$rm.batch.effect<-TRUE
      }

      #3. s3.plot.combatch
      if(!is.null(s3.plot.combatch.active())){
        actvalue$plot.combatch<-s3.plot.combatch.active()=="1"
      }else{
        actvalue$plot.combatch<-TRUE
      }


      #4. s4.cluster.method
      if(!is.null(s4.cluster.method.active())){
        actvalue$cluster.method<-ifelse(s4.cluster.method.active()=="1","RF","SVM")
      }else{
        actvalue$cluster.method<-"SVM"
      }

      #5. s5.kernel
      if(!is.null(s5.kernel.active())){
        actvalue$kernel<-ifelse(s5.kernel.active()=="1","radial",
                                ifelse(s5.kernel.active()=="2","linear",
                                       ifelse(s5.kernel.active()=="3",
                                              "polynomial","sigmoid")))
      }else{
        actvalue$kernel<-"radial"
      }

      #6. s6.plot.surv.curve
      if(!is.null(s6.plot.surv.curve.active())){
        actvalue$plot.surv.curve<-s6.plot.surv.curve.active()=="1"
      }else{
        actvalue$plot.surv.curve<-TRUE
      }

      #7. ntree
      if(!is.null(ntree.active())){
        actvalue$ntree<-ntree.active()
      }else{
        actvalue$ntree<-1000
      }


      #8. nodesize
      if(!is.null(nodesize.active())){
        actvalue$nodesize<-nodesize.active()
      }else{
        actvalue$nodesize<-3
      }

      #9. mtry
      if(!is.null(mtry.active())){
        actvalue$mtry<-mtry.active()
      }else{
        actvalue$mtry<-5
      }

      #10. cost
      if(!is.null(cost.active())){
        actvalue$cost<-cost.active()
      }else{
        actvalue$cost<-4
      }


      #11. survcurve.break.x.by

      if(!is.null(survcurve.break.x.by.active())){
        actvalue$survcurve.break.x.by<-survcurve.break.x.by.active()
      }else{
        actvalue$survcurve.break.x.by<-12
      }

      #12. file.PartSamples.Pred Whether to select parts of samples to perform prediction analysis
      if(!is.null(input$file.PartSamples.Pred)){
        actvalue$PartSamples.Pred.Path<-file.PartSamples.Pred.path1()
      }else{
          actvalue$PartSamples.Pred.Path<-NULL
      }

      #13. plot.samples 获取选择器输入框的值
      input.tumor.exp<-data.table::fread(file = actvalue$input.exp.path,
                                      data.table = F,showProgress = T)
      gene.column.guess<-grep(x=colnames(input.tumor.exp),pattern = "gene|feature|symbol",
           ignore.case = T,value = T)
      input.tumor.exp<-input.tumor.exp%>%tibble::column_to_rownames(var = gene.column.guess)

      #如果只分析部分samples，则需要进一步添加筛选条件：
      if(!is.null(input$file.PartSamples.Pred)){
        PartSamples.Pred<-data.table::fread(actvalue$PartSamples.Pred.Path,
                                            data.table = F,header = F)[,1]
        input.tumor.exp<-input.tumor.exp[PartSamples.Pred]
      }

      #如果只输入了3个样本，则不能选择绘制1:10个样本的生存曲线，所以需要使用head(1:num.max,10)
      num.max <- dim(input.tumor.exp)[2]

      if(!is.null(plot.samples.active())){
        p.sam.order.loc.end<-regexpr(plot.samples.active(),
                                     pattern = "(No\\.[0-9]+)")%>%
          attr("match.length")
        actvalue$plot.samples<-as.numeric(
          substr(plot.samples.active(),4,p.sam.order.loc.end))
      }else{
        actvalue$plot.samples<-head(1:num.max,10)
      }

      # 14. legend.size
      actvalue$leg.size<-leg.size.active()

      #15. sel.geneset #勾选基因集
      if(!is.null(sel.geneset.active())){
        actvalue$sel.geneset<-sel.geneset.active()
      }else{
        actvalue$sel.geneset<-"2"#默认为"2":Immunegene
      }

      #16. selgene.uplod #fileinput基因集
      if(!is.null(input$selgene.uplod)){
        actvalue$selgene.uplod<-selgene.uplod.active()
      }else{
        actvalue$selgene.uplod<-NULL
      }

      #17. text.input.gene.explore #输入text兴趣基因
      if(!is.null(text.input.gene.explore.active())&length(text.input.gene.explore.active())>0){
        actvalue$text.input.gene.explore<-text.input.gene.explore.active()
      }else{
        actvalue$text.input.gene.explore<-NULL
      }

      #18. select.hp.gene
      select.hp.gene.active
      if(!is.null(select.hp.gene.active())){
        actvalue$select.hp.gene<-select.hp.gene.active()
      }else{
        actvalue$select.hp.gene<-"1"
      }


      #19. survtime.unit #输入生存曲线时间单位
      if(!is.null(survtime.unit.active())){
        actvalue$survtime.unit<-survtime.unit.active()
      }else{
        actvalue$survtime.unit<-"months"
      }

      #20. hp.sample #热图是否选择sample进行绘制，是否显示sample名称
      if(!is.null(hp.sample.active())){
        actvalue$hp.sample<-hp.sample.active()
      }else{
        actvalue$hp.sample<-NULL
      }

      #21. limmaNorm.train Whether to perform limma::normalizeBetweenArrays()in the train cohort.
      if(!is.null(limmaNorm.train.active())){
        actvalue$train.exp.limma.normlize<-ifelse(limmaNorm.train.active()=="1",TRUE,FALSE)
      }else{
        actvalue$train.exp.limma.normlize<-FALSE
      }

      #22. limmaNorm.input Whether to perform limma::normalizeBetweenArrays()in the input cohort.
      if(!is.null(limmaNorm.input.active())){
        actvalue$input.exp.limma.normlize<-ifelse(limmaNorm.input.active()=="1",TRUE,FALSE)
      }else{
        actvalue$input.exp.limma.normlize<-FALSE
      }


      #23. RandomSamTrain Whether to use random samples to fit training models in the training cohort.
      if(!is.null(RandomSamTrain.active())){
        actvalue$RandomSamTrain<-ifelse(RandomSamTrain.active()=="1",TRUE,FALSE)
      }else{
        actvalue$RandomSamTrain<-FALSE
      }

      #24. TrainSeed If RandomSamTrain is TRUE, use a random seed to control the results of random sampling.
      if(!is.null(TrainSeed.active())){
        actvalue$TrainSeed<-TrainSeed.active()
      }else{
        actvalue$TrainSeed<-1234
      }

      #25. SamplingProb If RandomSamTrain is TRUE, provide a sampling probability. From 0.01 to 1
      SamplingProb.active<-eventReactive(input$run,{input$SamplingProb})
      if(!is.null(SamplingProb.active())){
        actvalue$SamplingProb<-SamplingProb.active()
      }else{
        actvalue$SamplingProb<-0.7
      }

      #perform the main program####
      actvalue$results<-oncoClassSurv(
        train.exp.path = actvalue$train.exp.path,
        input.exp.path = actvalue$input.exp.path,
        train.clin.path =actvalue$train.clin.path,
        input.clin.path = actvalue$input.clin.path,
        train_cluster.feature.path=actvalue$train_cluster.feature.path,
        train_survival.feature.path=actvalue$train_survival.feature.path,

        task=actvalue$task,
        rm.batch.effect=actvalue$rm.batch.effect,
        plot.combatch=actvalue$plot.combatch,
        cluster.method=actvalue$cluster.method,

        train.exp.limma.normlize =actvalue$train.exp.limma.normlize ,
        input.exp.limma.normlize =actvalue$input.exp.limma.normlize,

        random.sample.train=actvalue$RandomSamTrain,
        train.seeds=actvalue$TrainSeed,
        random.prob=actvalue$SamplingProb,
        PartSamples.Pred.Path=actvalue$PartSamples.Pred.Path,

        kernel = actvalue$kernel,
        cost = actvalue$cost,
        ntree = actvalue$ntree,
        nodesize = actvalue$nodesize,
        mtry = actvalue$mtry,

        plot.surv.curve=actvalue$plot.surv.curve,
        survcurve.break.x.by =actvalue$survcurve.break.x.by,
        plot.samples=actvalue$plot.samples,
        survtime.unit=actvalue$survtime.unit,

        print.combat.plots=FALSE,
        surv.t.custom=NULL,
        print.survplot = FALSE,
        show.message=FALSE)

      #显示警告信息（如果输入表达谱存在负值，需要进行警告）
       if(!is.null(actvalue$results$Neg.expression.warning)){
         showNotification(actvalue$results$Neg.expression.warning,
                          duration = 5,type = c("warning"))
       }

      if(actvalue$task%in%c(1,3)&actvalue$cluster.method=="SVM"){

        #for下文输出Class.pred表格
        Class.pred<-
          actvalue$results$svm.cluster$svm.cluster.pred%>%
          dplyr::mutate(Order=1:nrow(.))%>%
          dplyr::select(ncol(.),1:(ncol(.)-1))

        res_class.1<-actvalue$results$svm.cluster$svm.cluster.pred%>%
          dplyr::rename(Subtype=svm.cluster.pred)

        res_class.2<-actvalue$results$svm.cluster$svm.cluster.pred%>%
          dplyr::rename(Subtype=svm.cluster.pred)%>%
          dplyr::summarise(Count=as.numeric(table(Subtype)),
                           Allcount=sum(Count),
                           Proportion=scales::percent(Count/Allcount,0.01),
                           Subtype=levels(Subtype))

        fillcols<-mycol[1:length(unique(res_class.2$Subtype))]
        class.summary.p<-res_class.1 %>%
          ggplot(aes(factor(Subtype), fill = factor(Subtype))) +
          # 使用stat_count()函数来计算每个分组的计数，并且作为y坐标
          geom_bar(position = "stack", stat = "count", aes(y = ..count..)) +
          scale_fill_manual(values = fillcols) +
          guides(fill = guide_legend(title = "Subtypes")) +
          labs(title = paste0("Prediction for subtypes"),
               x = paste0("Subtypes", " (", actvalue$results$svm.cluster$cluster.method, ")"),
               y = "Patients count") +
          theme_bw() +
          theme(axis.text = element_text(color = "black"),
                axis.line = element_line(linewidth = 0),
                panel.border = element_rect(color="black",
                                            fill ="transparent")) +
          # 添加文字标签
          geom_text(data=res_class.2,aes(x=factor(Subtype),y=Count/2,
                                         label = paste(Count, Proportion, sep = "\n")), # 指定要显示的文本，用换行符分隔Count和Proportion
                    #position = position_stack(vjust = 0.5), # 指定文本的位置，居中于每个条形
                    size = 4, # 指定文本的大小
                    color = "white") # 指定文本的颜色

        ##3.2heatmap####
        #(3.2.1) genes.读取输入的基因集
        #list("Default" = 1, "Upload" = 2, "TextInput"=3)
        if(actvalue$select.hp.gene=="1"){
          if(actvalue$sel.geneset=="1"){
            cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
            input_gene<-cluster_markergenes
          }else{
            immunegenes<-readRDS(file = "./data/immunegenes.rds")
            immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
            input_gene<-immunegenes
          }
        }else{
          if(actvalue$select.hp.gene=="2"){
            if(!is.null(actvalue$selgene.uplod)){
              input_gene<-actvalue$selgene.uplod
            }else{
              if(actvalue$sel.geneset=="1"){
              cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
              input_gene<-cluster_markergenes
            }else{
              immunegenes<-readRDS(file = "./data/immunegenes.rds")
              immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
              input_gene<-immunegenes
            }}
          }else{
            if(actvalue$select.hp.gene=="3"){
              if(!is.null(actvalue$text.input.gene.explore)){
                input_gene<-actvalue$text.input.gene.explore
              }else{
                if(actvalue$sel.geneset=="1"){
                cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
                input_gene<-cluster_markergenes
              }else{
                immunegenes<-readRDS(file = "./data/immunegenes.rds")
                immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
                input_gene<-immunegenes
              }}
            }
          }
        }

        #显示通知信息1：
        if(actvalue$select.hp.gene=="2"&is.null(actvalue$selgene.uplod)){
          showNotification("No gene upload! Use default gene set.", duration = 5,type = c("warning"))
        }
        #显示通知信息2：
        if(actvalue$select.hp.gene=="2"&!is.null(actvalue$selgene.uplod)){
          check.hp_genes<-intersect(actvalue$selgene.uplod,rownames(input.tumor.exp))
          if(length(check.hp_genes)==0){
            showNotification("None upload gene matched! Use default gene set.", duration = 5,type = c("warning"))
          }
        }

        #显示通知信息3：
        if(actvalue$select.hp.gene=="3"&is.null(actvalue$text.input.gene.explore)){
          showNotification("No gene input! Use default gene set.", duration = 5,type = c("warning"))
        }
        #显示通知信息4：
        if(actvalue$select.hp.gene=="3"&!is.null(actvalue$text.input.gene.explore)){
          check.hp_genes<-intersect(actvalue$text.input.gene.explore,rownames(input.tumor.exp))
          if(length(check.hp_genes)==0){
            showNotification("None interesting gene matched! Use default gene set.", duration = 5,type = c("warning"))
          }
        }

        #(3.2.2) input.tumor.exp基因
        rownames(input.tumor.exp)<-gsub(x=rownames(input.tumor.exp),pattern = "-",replacement = "_")

        #(3.2.3) Common immune genes (输入基因集与input.tumor.exp基因取交集用于绘制热图)
        hp_genes<-intersect(input_gene,rownames(input.tumor.exp))

        #finial check hp_genes avoiding heat map with null gene. ###
        #If it happens (hp_genes is null), use the default gene set (immune or cluster-specific genes)
        if(length(hp_genes)==0){
          if(actvalue$sel.geneset=="1"){
            cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
            input_gene<-cluster_markergenes
          }else{
            immunegenes<-readRDS(file = "./data/immunegenes.rds")
            immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
            input_gene<-immunegenes
          }
          hp_genes<-intersect(input_gene,rownames(input.tumor.exp))
          if(length(hp_genes)==0){
            showNotification("Error: None gene matched, check your application!", duration = 5,type = c("error"))
          }
        }

        na.genes<-setdiff(input_gene,rownames(input.tumor.exp))
        if(length(na.genes)>0){
          na.warning<-paste0("Warning: The following genes were not found in the input expression profile file. Please check your input: ",
                 paste(na.genes,collapse = ","))
          showNotification(na.warning, duration = 5,type = c("warning"))
        }

        #(3.2.4) heatmap

        if(!is.null(actvalue$hp.sample)){
          if("1"%in%actvalue$hp.sample){
            sam.subtype.order<-res_class.1[actvalue$plot.samples,]%>%#与生存曲线的plot.sample保持一致
              dplyr::arrange(Subtype)%>%dplyr::select(-ID)#同时也是列注释
          }else{
            sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#此行代码意义为展示所有sample的热图。
          }
          if("2"%in%actvalue$hp.sample){
            show.hp.sample.name<-TRUE
          }else{
            show.hp.sample.name<-FALSE
          }
        }else{
          sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#此行代码意义为展示所有sample的热图。
          show.hp.sample.name<-FALSE
        }
        #sam.subtype.order同时也是列注释。
        #sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#同时也是列注释
        input.tumor.exp_sam.subtype.order<-input.tumor.exp[hp_genes,rownames(sam.subtype.order)]%>%
          as.data.frame()
        rownames(input.tumor.exp_sam.subtype.order)<-hp_genes
        colnames(input.tumor.exp_sam.subtype.order)<-rownames(sam.subtype.order)

        input.tumor.exp_sam.subtype.order.log2<-log2(input.tumor.exp_sam.subtype.order+1)
        #列为gene时进行scale()归一化，如果只输入了1个sample，则不进行scale()
        if(length(rownames(sam.subtype.order))>=2){
          input.tumor.exp_sam.subtype.order.log2.t<-t(scale(t(input.tumor.exp_sam.subtype.order.log2)))
        }else{
          input.tumor.exp_sam.subtype.order.log2.t<-input.tumor.exp_sam.subtype.order.log2

        }

        #ComplexHeatmap::Heatmap绘制热图###
        #注释颜色
        ann_fillcolors<-fillcols
        names(ann_fillcolors)<-unique(res_class.2$Subtype)

        #重设热图基因标签的size，不能过大过小
        hmp.g.size<-ifelse(150/length(hp_genes)>15,15,
                           ifelse(150/length(hp_genes)>3,150/length(hp_genes),3))
        #Warning: Error in hclust: 用群集时必需有n >= 2的对象
        cluster_rows<-ifelse(length(hp_genes)>=2, TRUE,FALSE)

        cluster.Top <- HeatmapAnnotation(Subtype= sam.subtype.order$Subtype,

                                         annotation_legend_param=list(labels_gp = gpar(fontsize=10),
                                                                      title_gp=gpar(fontsize =10,fontface="bold"),
                                                                      ncol=1),
                                         border = FALSE,
                                         col=list(Subtype= ann_fillcolors) ,
                                         show_annotation_name = TRUE ,
                                         annotation_name_side="right",
                                         annotation_name_gp = gpar(fontsize = 10))

        hmp<-Heatmap(input.tumor.exp_sam.subtype.order.log2.t, ## 对基因表达量进行限定，也可以进行转置t后scale进行过滤
                     name='Gene' ,#热图侧名
                     top_annotation = cluster.Top,#顶部注释
                     cluster_rows = cluster_rows,# 是否对行进行聚类
                     show_row_dend = TRUE,#是否显示聚类树
                     col=colorRamp2(c(-1,0,1) ,c('#008B8B' ,'#F5F5F5', '#DC143C')),#颜色steelblue
                     color_space ="RGB",
                     row_dend_width = unit(6,"mm"),
                     cluster_columns = FALSE,
                     border = FALSE,
                     row_order=NULL,
                     row_names_side = 'right' ,
                     column_order=NULL,
                     show_column_names = show.hp.sample.name,
                     show_row_names = TRUE ,
                     row_names_gp =gpar(fontsize = hmp.g.size),
                     gap=unit(4,"mm"),
                     column_title = NULL,
                     column_title_gp=gpar(fontsize = 10),
                     row_title_gp = gpar(fontsize=10) ,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param=list(
                       labels_gp = gpar(fontsize = 10) ,
                       title_gp =gpar(fontsize = 10,
                                      fontface = "bold"),
                       legend_direction="vertical",
                       title_position="topleft"))

        gg.hmp <- ggplotify::as.ggplot(hmp)
      }


      if(actvalue$task%in%c(1,3)&actvalue$cluster.method=="RF"){
        #for下文输出Class.pred表格
        Class.pred<-
          actvalue$results$rf.cluster$rf.cluster.pred%>%
          dplyr::mutate(Order=1:nrow(.))%>%
          dplyr::select(ncol(.),1:(ncol(.)-1))

        res_class.1<-actvalue$results$rf.cluster$rf.cluster.pred%>%
          dplyr::rename(Subtype=rf.cluster.pred)

        res_class.2<-actvalue$results$rf.cluster$rf.cluster.pred%>%
          dplyr::rename(Subtype=rf.cluster.pred)%>%
          dplyr::summarise(Count=as.numeric(table(Subtype)),
                           Allcount=sum(Count),
                           Proportion=scales::percent(Count/Allcount,0.01),
                           Subtype=levels(Subtype))

        fillcols<-mycol[1:length(unique(res_class.2$Subtype))]
        class.summary.p<-res_class.1 %>%
          ggplot(aes(factor(Subtype), fill = factor(Subtype))) +
          # 使用stat_count()函数来计算每个分组的计数，并且作为y坐标
          geom_bar(position = "stack", stat = "count", aes(y = ..count..)) +
          scale_fill_manual(values = fillcols) +
          guides(fill = guide_legend(title = "Subtypes")) +
          labs(title = paste0("Prediction for subtypes"),
               x = paste0("Subtypes", " (", actvalue$results$rf.cluster$cluster.method, ")"),
               y = "Patients count") +
          theme_bw() +
          theme(axis.text = element_text(color = "black"),
                axis.line = element_line(linewidth = 0),
                panel.border = element_rect(color="black",
                                            fill ="transparent")) +
          # 添加文字标签
          geom_text(data=res_class.2,aes(x=factor(Subtype),y=Count/2,
                                         label = paste(Count, Proportion, sep = "\n")), # 指定要显示的文本，用换行符分隔Count和Proportion
                    #position = position_stack(vjust = 0.5), # 指定文本的位置，居中于每个条形
                    size = 4, # 指定文本的大小
                    color = "white") # 指定文本的颜色

        ##3.2heatmap####
        #(3.2.1) genes.读取输入的基因集
        #list("Default" = 1, "Upload" = 2, "TextInput"=3)
        if(actvalue$select.hp.gene=="1"){
          if(actvalue$sel.geneset=="1"){
            cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
            input_gene<-cluster_markergenes
          }else{
            immunegenes<-readRDS(file = "./data/immunegenes.rds")
            immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
            input_gene<-immunegenes
          }
        }else{
          if(actvalue$select.hp.gene=="2"){
            if(!is.null(actvalue$selgene.uplod)){
              input_gene<-actvalue$selgene.uplod
            }else{
              if(actvalue$sel.geneset=="1"){
                cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
                input_gene<-cluster_markergenes
              }else{
                immunegenes<-readRDS(file = "./data/immunegenes.rds")
                immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
                input_gene<-immunegenes
              }}
          }else{
            if(actvalue$select.hp.gene=="3"){
              if(!is.null(actvalue$text.input.gene.explore)){
                input_gene<-actvalue$text.input.gene.explore
              }else{
                if(actvalue$sel.geneset=="1"){
                  cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
                  input_gene<-cluster_markergenes
                }else{
                  immunegenes<-readRDS(file = "./data/immunegenes.rds")
                  immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
                  input_gene<-immunegenes
                }}
            }
          }
        }

        #显示通知信息1：
        if(actvalue$select.hp.gene=="2"&is.null(actvalue$selgene.uplod)){
          showNotification("No gene upload! Use default gene set.", duration = 5,type = c("warning"))
        }
        #显示通知信息2：
        if(actvalue$select.hp.gene=="2"&!is.null(actvalue$selgene.uplod)){
          check.hp_genes<-intersect(actvalue$selgene.uplod,rownames(input.tumor.exp))
          if(length(check.hp_genes)==0){
            showNotification("None upload gene matched! Use default gene set.", duration = 5,type = c("warning"))
          }
        }

        #显示通知信息3：
        if(actvalue$select.hp.gene=="3"&is.null(actvalue$text.input.gene.explore)){
          showNotification("No gene input! Use default gene set.", duration = 5,type = c("warning"))
        }
        #显示通知信息4：
        if(actvalue$select.hp.gene=="3"&!is.null(actvalue$text.input.gene.explore)){
          check.hp_genes<-intersect(actvalue$text.input.gene.explore,rownames(input.tumor.exp))
          if(length(check.hp_genes)==0){
            showNotification("None interesting gene matched! Use default gene set.", duration = 5,type = c("warning"))
          }
        }

        #(3.2.2) input.tumor.exp基因
        rownames(input.tumor.exp)<-gsub(x=rownames(input.tumor.exp),pattern = "-",replacement = "_")

        #(3.2.3) Common immune genes (输入基因集与input.tumor.exp基因取交集用于绘制热图)
        hp_genes<-intersect(input_gene,rownames(input.tumor.exp))

        #finial check hp_genes avoiding heat map with null gene. ###
        #If it happens (hp_genes is null), use the default gene set (immune or cluster-specific genes)
        if(length(hp_genes)==0){
          if(actvalue$sel.geneset=="1"){
            cluster_markergenes <-readRDS(actvalue$train_cluster.feature.path)
            input_gene<-cluster_markergenes
          }else{
            immunegenes<-readRDS(file = "./data/immunegenes.rds")
            immunegenes<-gsub(x=immunegenes,pattern = "-",replacement = "_")
            input_gene<-immunegenes
          }
          hp_genes<-intersect(input_gene,rownames(input.tumor.exp))
          if(length(hp_genes)==0){
            showNotification("Error: None gene matched, check your application!", duration = 5,type = c("error"))
          }
        }

        na.genes<-setdiff(input_gene,rownames(input.tumor.exp))
        if(length(na.genes)>0){
          na.warning<-paste0("Warning: The following genes were not found in the input expression profile file. Please check your input: ",
                             paste(na.genes,collapse = ","))
          showNotification(na.warning, duration = 5)
        }

        #(3.2.4) heatmap
        if(!is.null(actvalue$hp.sample)){
          if("1"%in%actvalue$hp.sample){
            sam.subtype.order<-res_class.1[actvalue$plot.samples,]%>%#与生存曲线的plot.sample保持一致
              dplyr::arrange(Subtype)%>%dplyr::select(-ID)#同时也是列注释
          }else{
            sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#此行代码意义为展示所有sample的热图。
          }
          if("2"%in%actvalue$hp.sample){
            show.hp.sample.name<-TRUE
          }else{
            show.hp.sample.name<-FALSE
          }
        }else{
          sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#此行代码意义为展示所有sample的热图。
          show.hp.sample.name<-FALSE
        }
        #sam.subtype.order同时也是列注释。
        #sam.subtype.order<-res_class.1%>%dplyr::arrange(Subtype)%>%dplyr::select(-ID)#同时也是列注释
        input.tumor.exp_sam.subtype.order<-input.tumor.exp[hp_genes,rownames(sam.subtype.order)]%>%
          as.data.frame()
        rownames(input.tumor.exp_sam.subtype.order)<-hp_genes
        colnames(input.tumor.exp_sam.subtype.order)<-rownames(sam.subtype.order)

        input.tumor.exp_sam.subtype.order.log2<-log2(input.tumor.exp_sam.subtype.order+1)
        #列为gene时进行scale()归一化，如果只输入了1个sample，则不进行scale()
        if(length(rownames(sam.subtype.order))>=2){
          input.tumor.exp_sam.subtype.order.log2.t<-t(scale(t(input.tumor.exp_sam.subtype.order.log2)))
        }else{
          input.tumor.exp_sam.subtype.order.log2.t<-input.tumor.exp_sam.subtype.order.log2

        }

        #ComplexHeatmap::Heatmap绘制热图###
        #注释颜色
        ann_fillcolors<-fillcols
        names(ann_fillcolors)<-unique(res_class.2$Subtype)

        #重设热图基因标签的size，不能过大过小
        hmp.g.size<-ifelse(150/length(hp_genes)>15,15,
                           ifelse(150/length(hp_genes)>3,150/length(hp_genes),3))
        #Warning: Error in hclust: 用群集时必需有n >= 2的对象
        cluster_rows<-ifelse(length(hp_genes)>=2, TRUE,FALSE)

        cluster.Top <- HeatmapAnnotation(Subtype= sam.subtype.order$Subtype,

                                         annotation_legend_param=list(labels_gp = gpar(fontsize=10),
                                                                      title_gp=gpar(fontsize =10,fontface="bold"),
                                                                      ncol=1),
                                         border = FALSE,
                                         col=list(Subtype= ann_fillcolors) ,
                                         show_annotation_name = TRUE ,
                                         annotation_name_side="right",
                                         annotation_name_gp = gpar(fontsize = 10))

        hmp<-Heatmap(input.tumor.exp_sam.subtype.order.log2.t, ## 对基因表达量进行限定，也可以进行转置t后scale进行过滤
                     name='Gene' ,#热图侧名
                     top_annotation = cluster.Top,#顶部注释
                     cluster_rows = cluster_rows,# 是否对行进行聚类
                     show_row_dend = TRUE,#是否显示聚类树
                     col=colorRamp2(c(-1,0,1) ,c('#008B8B' ,'#F5F5F5', '#DC143C')),#颜色steelblue
                     color_space ="RGB",
                     row_dend_width = unit(6,"mm"),
                     cluster_columns = FALSE,
                     border = FALSE,
                     row_order=NULL,
                     row_names_side = 'right' ,
                     column_order=NULL,
                     show_column_names = show.hp.sample.name,
                     show_row_names = TRUE ,
                     row_names_gp =gpar(fontsize = hmp.g.size),
                     gap=unit(4,"mm"),
                     column_title = NULL,
                     column_title_gp=gpar(fontsize = 10),
                     row_title_gp = gpar(fontsize=10) ,
                     show_heatmap_legend = TRUE,
                     heatmap_legend_param=list(
                       labels_gp = gpar(fontsize = 10) ,
                       title_gp =gpar(fontsize = 10,
                                      fontface = "bold"),
                       legend_direction="vertical",
                       title_position="topleft"))

        gg.hmp <- ggplotify::as.ggplot(hmp)
      }


      if(actvalue$task%in%c(1,3)&actvalue$cluster.method=="SVM"){
        output$class.method <- renderText({
          paste0("The current classfication method: ",actvalue$results$svm.cluster$cluster.method)  # 根据你的输入生成文本
        })

      }

      if(actvalue$task%in%c(1,3)&actvalue$cluster.method=="RF"){
        output$class.method <- renderText({
          paste0("The current classfication method: ",actvalue$results$rf.cluster$cluster.method)  # 根据你的输入生成文本
        })
      }

      if(actvalue$task==2){
        output$class.method <- renderText({
          NULL
        })
      }

      #task包含分型
      if(actvalue$task%in%c(1,3)){
        #亚型统计
        output$Plot_class.summary <- renderPlot({
          class.summary.p
        })
        #不同分型的基因表达热图
        output$Plot_hp <- renderPlot({
          gg.hmp
        })
      }else{
        #亚型统计
        output$Plot_class.summary <- renderPlot({
          NULL
        })
        #不同分型的基因表达热图
        output$Plot_hp <- renderPlot({
          NULL
        })

      }

      #task包含预后
      if(actvalue$task%in%c(2,3)){
        #生存曲线
        surv.plots<-actvalue$results$ggsurv.curve$ggsurv.curve$plot+
          ggplot2::theme(legend.text = element_text(size = actvalue$leg.size))
        output$ggsurv.curve <- renderPlot({
          surv.plots
        })
      }else{
        output$ggsurv.curve <- renderPlot({
          NULL
        })
      }

      #task包含分型
      if(actvalue$task%in%c(1,3)){
        #分型表格
        output$class.table <- renderTable(
          Class.pred
        )
      }else{
        output$class.table <- renderTable(
          NULL
        )
      }

      #task包含预后
      if(actvalue$task%in%c(2,3)){
        #预后表格
        output$surv.table <- renderTable(
          actvalue$results$surv.probablity
        )
      }else{
        output$surv.table <- renderTable(
          NULL
        )
      }

      #绘制overview:combat PCA
      combat.pca<-actvalue$results$original_combat.plots
      output$Combat <- renderPlot({
        combat.pca
      })

    ## Our dataset
    xtime<-paste0(substr(Sys.time(),1,10)," ",
                  paste0(unlist(strsplit(substring(Sys.time(),12),split = ":")),
                         c("h","m","s"),
                         collapse = ""))

    #task包含分型
    if(actvalue$task%in%c(1,3)){

      output$download.class <- downloadHandler(
        filename = function() {
          paste("Classification_Prediction_Results_", xtime, ".csv", sep="")
        },
        content = function(file) {
          write.csv(Class.pred, file)
        }
      )}else{
      output$download.class <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )
    }

    #task包含预后
    if(actvalue$task%in%c(2,3)){
      surv.probablity<-actvalue$results$surv.probablity
      output$download.surv <- downloadHandler(
        filename = function() {
          paste("Prognosis_Prediction_Results_", xtime, ".csv", sep="")
        },
        content = function(file) {
          write.csv(surv.probablity, file)
        }
      )
    }else{
      output$download.surv <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )
    }

    #task包含分型
    if(actvalue$task%in%c(1,3)){
      #download.class.count
      output$download.class.count <- downloadHandler(
        filename = function() {
          paste("Classification_Prediction_Barplot_", xtime, ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          print(class.summary.p)#输出pdf图片必须加print()，否则输出图片报错。
          dev.off()
        }
      )

      #download.class.hp
      output$download.class.hp <- downloadHandler(
        filename = function() {
          paste("Classification_Prediction_Expression_Heat_Map_", xtime, ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          print(gg.hmp)
          dev.off()
        }
      )
    }else{
      #download.class.count
      output$download.class.count <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )

      #download.class.hp
      output$download.class.hp <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )
    }

    #task包含预后
    if(actvalue$task%in%c(2,3)){
      #download.survplot
      output$download.survplot <- downloadHandler(
        filename = function() {
          paste("Prognosis_Prediction_Survival_Curves_", xtime, ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          print(surv.plots)
          dev.off()
        }
      )
    }else{
      #download.survplot
      output$download.survplot <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )
    }

    #download.combat
    if(!is.null(combat.pca)){
      output$download.combat <- downloadHandler(
        filename = function() {
          paste("Remove_Batch_Effect_PCAplot_", xtime, ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          print(combat.pca)
          dev.off()
        }
      )
    }else{
      output$download.combat <- downloadHandler(
        filename = function() {
          "Blank"
        },
        content = function(file) {
          NULL
        }
      )
    }

    #初始全部计算之后，统计sample数量，sample名称，sample分型情况，并更新UI中选择器的sample选择范围
    input.sample.names<-colnames(input.tumor.exp)

    #如果task=1或3，则把患者的分型情况也更新到选择器中
    if(actvalue$task%in%c(1,3)){
      # 更新选择器输入框的选项，使其等于最小值和最大值之间的序列
      updatePickerInput(session, "plot.samples", choices = paste0("No.",1:num.max,": ",input.sample.names,": ",Class.pred[,3]),
                        selected = input$plot.samples)
      #selected = input$plot.samples保证选择器中选择的数值不被实时清空
    }else{
      # 更新选择器输入框的选项，使其等于最小值和最大值之间的序列
      updatePickerInput(session, "plot.samples", choices = paste0("No.",1:num.max,": ",input.sample.names),
                        selected = input$plot.samples)
      #selected = input$plot.samples保证选择器中选择的数值不被实时清空
    }

    hidePageSpinner()

  })

  # Reset parameters
  observeEvent(input$reset_button, {
      js$resetPage()
  })

})
