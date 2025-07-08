#' This function draw table
#'
#' @param df data, gvar:group,varlist:variate,skewvar:Skew variable,cat_var:cat_var,out_docx:filename.
#' @return table.
#' @export
#' 


#更新fisher
#修复Z统计量
data_statistic <- function (df, gvar= NULL, varlist = NULL, skewvar = NULL,
                            cat_var=NULL,out_docx,path,
                            pnormtest =0.05){
  start_time <- Sys.time()
  if(!require("stringr")) install.packages("stringr")
  if(!require("dplyr")) install.packages("dplyr")
  library("stringr")
  library("dplyr")

  tt1 <- 0; zz1 <- 0
  tt2 <- 0; zz2 <- 0
  ff2 <- 0; kw2 <- 0
  cc2 <- 0; FH2 <- 0
  
  str_var <- function(temp_variable){
    temp_variable <- gsub('_', ' ', temp_variable, fixed=TRUE)
    temp_variable <- gsub('.', ' ', temp_variable, fixed=TRUE)
    temp_variable <- stringr::str_squish(string = temp_variable)
    if (stringr::str_to_upper(substr(temp_variable, 1, 1))!=
        stringr::str_trunc(temp_variable,1,side = c("right"),ellipsis="")){
      stringr::str_to_title(temp_variable)
    }
    else {temp_variable}
  }

  ####去掉分组缺失####
  if (length(gvar>0)){
    dfs <- subset(df,!is.na(df[,gvar]))
    dfs[, gvar] <- as.factor(dfs[, gvar])
    if (length(table(dfs[, gvar])) > 2) {
      print(paste0("本次分析分组变量",str_to_title(gvar),"为多分类变量"),quote=F)
    }
    name_t <- paste0("Total"," (n = ",nrow(dfs),")")
    name_g <- paste0(names(table(dfs[, gvar]))," (n = ",table(dfs[, gvar]),")")
  }
  if (length(gvar)==0){
    dfs <- df
    name_t <- paste0("Total"," (n = ",nrow(dfs),")")
  }
  
  if (length(cat_var)>0){
    for (f in 1:length(cat_var)){
      dfs[,cat_var[f]] <- factor(dfs[,cat_var[f]])
    }
  }
  
  
  table <- list()
  for (i in 1:length(varlist)){
    ####开始循环####
    if (class(dfs[, varlist[i]]) == "factor" | class(dfs[, varlist[i]]) == "character"){
      if (length(levels(factor(dfs[, varlist[i]]))) >10){
        print(paste0("变量：",varlist[i],"分类水平超过了10，请确认数据！"),quote=F)
      }
      dfs[,varlist[i]] <- factor(dfs[,varlist[i]])
      dfa <- subset(dfs,!is.na(dfs[,varlist[i]]))
      
      tableTol <- table(dfa[, varlist[i]])
      per <- prop.table(tableTol)
      prop=paste(as.data.frame(tableTol)[,"Freq"], " (",round(as.data.frame(per)[,"Freq"] * 100, 2), ")", sep = "")%>% 
        as.data.frame()%>% as_tibble()
      if (length(gvar)>0){
        ####分类的资料####
        table.sub <- table(dfa[, varlist[i]], dfa[, gvar])
        options (warn=-1)
        expected_num <- chisq.test(dfa[, varlist[i]], dfa[, gvar])$expected
        ####满足期望频数的类型
        e5 <- ifelse(expected_num>=5,0,1)#检查频数是否<5
        e1 <- ifelse(expected_num>=1,0,1)#检查频数是否<1
        
        prob <- list()
        for (g in 1:length(table(dfa[,gvar]))){
          prob[g]=paste0(table.sub[,g]," (",sprintf('%0.2f',(table.sub[,g]/sum(table.sub[,g]))*100),")") %>% list
        }
        var_grp <- do.call(cbind,prob)%>% as.data.frame()%>% as_tibble()
        if (nrow(table.sub) == 1) {
          P = 1
          statistic = NULL
        }
        else {
          ####卡方检验优化 
          if (length(levels(factor(dfa[, varlist[i]])))==2 & length(levels(factor(dfa[, gvar])))==2){ 
            #### 2*2
            if (sum(expected_num)>=40 & sum(e5)==0){
              cc2 <- cc2 + 1
              print(paste0("正在分析:",str_var(varlist[i]),"进行卡方检验"),quote=F)
              #n ≥ 40，所有T ≥ 5，选择Pearson卡方检验 (非连续校正)
              Statistic <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = F)$statistic
              P <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = F)$p.value
            }
            else if (sum(expected_num)>=40 & sum(e1)==0 & sum(e5)>0){
              cc2 <- cc2 + 1
              print(paste0("正在分析:",str_var(varlist[i]),"进行校正卡方检验"),quote=F)
              #n ≥ 40，至少一个单元格1 ≤ T＜5，选择连续校正卡方检验；
              Statistic <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = T)$statistic
              P <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = T)$p.value
            }
            else if (sum(expected_num)<40 | sum(e1)>0){
              FH2 <- FH2 + 1
              print(paste0("正在分析:",str_var(varlist[i]),"进行Fisher确切概率法"),quote=F)
              #n＜40或T＜1，选择Fisher精确检验
              Statistic <- NULL
              P <- fisher.test(dfa[, varlist[i]], dfa[, gvar])$p.value
            }
          }
          if (length(levels(factor(dfa[, varlist[i]])))*length(levels(factor(dfa[, gvar])))>4){
            ####R*C
            if (sum(expected_num)>=40 & 
                sum(e5)/length(levels(factor(dfa[, varlist[i]])))*length(levels(factor(dfa[, gvar])))<=0.2){
              cc2 <- cc2 + 1
              print(paste0("正在分析:",str_var(varlist[i]),"进行卡方检验"),quote=F)
              #n ≥ 40，所有T ≥ 5，选择Pearson卡方检验 (非连续校正)
              Statistic <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = F)$statistic
              P <- chisq.test(dfa[, varlist[i]], dfa[, gvar],correct = F)$p.value
            }
            else{
              #n＜40或T＜1，选择Fisher精确检验
              FH2 <- FH2 + 1
              print(paste0("正在分析:",str_var(varlist[i]),"进行Fisher确切概率法"),quote=F)
              Statistic <- NULL
              P <- fisher.test(dfa[, varlist[i]], dfa[, gvar],simulate.p.value=TRUE)$p.value
            }
            options(warn=1)
          }
        }
        options(warn=1)
        newline <- c(paste(str_var(varlist[i]), ", n (%)", sep = ""),
                     rep("", length(table(dfa[, gvar]))+1),
                     ifelse(is.null(Statistic), "-", paste0("χ²=",sprintf('%0.3f',Statistic))),
                     ifelse(P <0.001, "<.001",sprintf('%0.3f',P)))
        table_temp <- cbind(Variable = paste0("  ",levels(factor(dfs[, varlist[i]]))), 
                            setNames(prop, name_t)%>% as.data.frame()%>% as_tibble(),
                            setNames(var_grp, name_g)%>% as.data.frame()%>% as_tibble(),
                            Statistic = "",
                            P = "")%>% as.data.frame()%>% as_tibble()
        rm(prop,var_grp)
        table[i] <- rbind(newline,table_temp) %>% list()
        rm(table_temp )
      }
      if (length(gvar)==0){
        newline <- c(paste(str_var(varlist[i]), ", n (%)", sep = ""),
                     rep("",1))
        table_temp <- cbind(Variable = paste0("  ",levels(factor(dfs[, varlist[i]]))), 
                            setNames(prop, name_t)%>% as.data.frame()%>% as_tibble())%>% 
          as.data.frame()%>% as_tibble()
        table[i] <- rbind(newline,table_temp) %>% list()
        rm(table_temp )
      }
    }
    else if (class(dfs[, varlist[i]]) == "integer" | class(dfs[, varlist[i]]) == "numeric"){
      library("nortest")
      dfa <- subset(dfs,!is.na(dfs[,varlist[i]]))
      dfa$subvar <-  (dfa[,which(colnames(dfa)==varlist[i])]) %>% as.numeric()
      if ((ad.test(dfa[, varlist[i]])$p.value >= pnormtest) | (!(varlist[i] %in% skewvar))) {
        tt1 <- tt1 + 1
        MS <- dfa %>%
          dplyr::summarise(mean=sprintf("%0.2f",mean(subvar)),
                           sd=sprintf("%0.2f",sd(subvar)),
                           ms=paste0(mean," ± ",sd)) %>% select(ms) %>% 
          t %>% as.data.frame()%>% as_tibble()
        if (length(gvar)>0){
          MS_g <- dfa %>% group_by(dfa[,gvar]) %>%
            dplyr::summarise(mean=sprintf("%0.2f",mean(subvar)),
                             sd=sprintf("%0.2f",sd(subvar)),
                             ms=paste0(mean," ± ",sd)) %>% select(ms) %>% 
            t %>% as.data.frame()%>% as_tibble()
          if (length(table(dfa[, gvar])) == 2){
            #### t ####
            print(paste0("正在分析:",str_var(varlist[i]),"进行t检验"),quote=F)
            tt2 <- tt2 + 1
            library("car")
            var_test <- leveneTest(dfa[, varlist[i]],dfa[, gvar])
            if (var_test$`Pr(>F)`[1]>0.05){
              Statistic <- t.test(dfa[, varlist[i]] ~ dfa[, gvar],var.equal  =T)$statistic
              p <- t.test(dfa[, varlist[i]] ~ dfa[, gvar],var.equal  =T)$p.value
            }
            else if (var_test$`Pr(>F)`[1]<=0.05){
              Statistic <- t.test(dfa[, varlist[i]] ~ dfa[, gvar],var.equal =F)$statistic
              p <- t.test(dfa[, varlist[i]] ~ dfa[, gvar],var.equal =F)$p.value
            }
            table[i] <- cbind(Variable =paste0(str_var(varlist[i]),", ","Mean ± SD"),
                              setNames(MS, name_t)%>% as.data.frame()%>% as_tibble(),  #合计
                              setNames(MS_g, name_g)%>% as.data.frame()%>% as_tibble(), #分组 
                              Statistic =paste0("t=",sprintf('%0.3f',Statistic)),
                              P = ifelse(p < 0.001, "<.001", sprintf('%0.3f',p))) %>% list()
          }
          else if (length(table(dfa[, gvar])) > 2){
            #### F ####
            ff2 <- ff2 + 1
            print(paste0("正在分析:",str_var(varlist[i]),"进行方差分析"),quote=F)
            ANOVA <- aov(dfa[, varlist[i]]~dfa[, gvar], dfa)
            Statistic <- summary(ANOVA)[[1]][1,4]
            p <- summary(ANOVA)[[1]][1,5]
            table[i] <- cbind(Variable =paste0(str_var(varlist[i]),", ","Mean ± SD"),
                              setNames(MS, name_t)%>% as.data.frame()%>% as_tibble(),  #合计
                              setNames(MS_g, name_g)%>% as.data.frame()%>% as_tibble(), #分组 
                              Statistic =paste0("F=",sprintf('%0.3f',Statistic)),
                              P = ifelse(p < 0.001,"<.001", sprintf('%0.3f',p))) %>% list()
            rm(MS,MS_g )
          }
        }
        if (length(gvar)==0){
          table[i] <- cbind(Variable =paste0(str_var(varlist[i]),", ","Mean ± SD"),
                            setNames(MS, name_t)%>% 
                              as.data.frame()%>% as_tibble()) %>% list()
        }
      }
      else if ((ad.test(dfa[, varlist[i]])$p.value < pnormtest) | (varlist[i] %in% skewvar)){
        zz1 <- zz1 + 1
        MQ <-dfa  %>% 
          dplyr::summarise(median=sprintf("%0.2f",quantile(subvar,0.5)),
                           Q1=sprintf("%0.2f",quantile(subvar,0.25)),
                           Q3=sprintf("%0.2f",quantile(subvar,0.75))) 
        MQ <- paste0(MQ[1]," (",MQ[2],", ",MQ[3],")")%>% as.data.frame()%>% as_tibble()
        MQ<- setNames(MQ, name_t)%>% as.data.frame()%>% as_tibble()
        if (length(gvar)>0){
          ####分组####
          qrange <- dfa %>% group_by(dfa[,gvar]) %>%  
            dplyr::summarise(median=sprintf("%0.2f",quantile(subvar,0.5)),
                             Q1=sprintf("%0.2f",quantile(subvar,0.25)),
                             Q3=sprintf("%0.2f",quantile(subvar,0.75))) %>% t
          grp_stat <- list()
          for (g in 1:length(table(dfa[, gvar]))){
            grp_stat[g] <- paste0(qrange[2,g]," (",qrange[3,g],", ",qrange[4,g],")")%>% as.data.frame()%>% as_tibble()
          }
          var_grp <- do.call(cbind,grp_stat)%>% as.data.frame()%>% as_tibble()
          var_grp<- setNames(var_grp, name_g)%>% as.data.frame()%>% as_tibble()
          if (length(table(dfa[, gvar])) == 2){
            #### Z ####
            zz2 <- zz2 + 1
            print(paste0("正在分析:",str_var(varlist[i]),"进行秩和检验"),quote=F)
            p <- wilcox.test(dfa[, varlist[i]] ~ dfa[, gvar])$p.value
            Z=qnorm(p/2)
            table[i] <-cbind(Variable =paste0(str_var(varlist[i]),", ","M (Q₁, Q₃)"),
                             MQ,var_grp,
                             Statistic =paste0("Z=",sprintf('%0.3f',Z)),
                             P = ifelse(p <0.001, "<.001", sprintf('%0.3f',p))) %>% list()
            rm(MQ,var_grp )
          }
          else if (length(table(dfa[, gvar])) >= 2){
            #### kw ####
            kw2 <- kw2 + 1
            print(paste0("正在分析:",str_var(varlist[i]),"进行多组的秩和检验"),quote=F)
            a <- kruskal.test(dfa[, varlist[i]]~dfa[, gvar],data=dfa)
            KW <- a$statistic
            p <- a$p.value
            table[i]<-cbind(Variable =paste0(str_var(varlist[i]),", ","M (Q₁, Q₃)"),
                            MQ,var_grp,
                            Statistic =paste0("χ²=",sprintf('%0.3f',KW),"#"),
                            P = ifelse(p < 0.001, "<.001", sprintf('%0.3f',p))) %>% list()
            rm(MQ,var_grp )
          }
        }
        if (length(gvar)==0){
          table[i] <-cbind(Variable =paste0(str_var(varlist[i]),", ","M (Q₁, Q₃)"),
                           MQ) %>% list()
        }
      }
    }
  }
  final_temp <- do.call(rbind,table) %>% as_tibble()
  if (length(gvar>0)){
    final_temp <- final_temp %>% 
      mutate(id=ifelse(substr(P, 1, 1)=="<",1,
                       ifelse(as.numeric(gsub('<', '0', P, fixed=TRUE))<0.05 & P !="",1,0))) %>% as_tibble()
  }
  
  final <<- final_temp
  if(!require("flextable")) install.packages("flextable")
  library("flextable")
  
  ####描述脚注####
  tt1_txt <- "" ; zz1_txt <- ""
  if (tt1>0){tt1_txt <- ", SD: standard deviation"}
  if (zz1>0){zz1_txt <- ", M: Median, Q₁: 1st Quartile, Q₃: 3st Quartil"}
  ####方法学脚注####
  tt2_txt <- "" ; zz2_txt <- ""
  ff2_txt <- "" ; kw2_txt <- ""
  cc2_txt <- "" ; FH2_txt <- ""
  if (tt2>0){tt2_txt <- ", t: t-test"}
  if (zz2>0){zz2_txt <- ", Z: Mann-Whitney test"}
  if (ff2>0){ff2_txt <- ", F: ANOVA"}
  if (kw2>0){kw2_txt <- ", #: Kruskal-waills test"}
  if (cc2>0){cc2_txt <- ", χ²: Chi-square test"}
  if (FH2>0){FH2_txt <- ", -: Fisher exact"}

  name <- names(final)
  name <- c(name[!name %in% "id"])
  
  c <- flextable(
    final,
    col_keys = name
  )
  t_head <- names(final)[2]
  ####描述的脚注####
  txt <- paste0(tt1_txt,zz1_txt)#拼接
  txt_length <- nchar(txt, type='bytes')#计算长度
  if (txt_length>0){
    txt <- substr(txt,2,txt_length)#把第一个截掉
    c = add_footer_lines(c, values = txt)#插入脚注
  }
  ####方法学脚注####
  if (length(gvar>0)){
    x <- rep("",length(table(df[, gvar]))-1)
    txt2 <- paste0(tt2_txt,zz2_txt,ff2_txt,kw2_txt,cc2_txt,FH2_txt)#拼接
    txt2_length <- nchar(txt2, type='bytes')#计算长度
    txt2 <- substr(txt2,2,txt2_length)#把第一个截掉
    ITA_COL <- length(table(df[, gvar]))+4
    c = flextable::add_footer_lines(c, values = txt2) %>% #插入脚注
      flextable::italic(j=ITA_COL,part="header") %>% 
      flextable::color(color = "red",i = ~ id ==1,j=~P) %>%  #标红
      flextable::bold(bold = TRUE,i = ~ id ==1,j=~P) %>%  #加粗
      flextable::autofit(add_w = 0.1, add_h = 0.1, part = c("body", "header")) %>% 
      add_header_row(
        top = TRUE, 
        values = c("Variable",
                   t_head, 
                   str_var(gvar),
                   x, 
                   "Statistic", 
                   "P")) %>% 
      align(align = "center", j = c(3:(length(table(df[, gvar]))+2)), part = "header") %>%
      merge_at(i = 1:2, j = 1, part = "header") %>% 
      merge_at(i = 1:2, j = 2, part = "header") %>% 
      merge_at(i = 1:2, j = length(table(df[, gvar]))+3, part = "header") %>% 
      merge_at(i = 1:2, j = length(table(df[, gvar]))+4, part = "header") %>% 
      merge_at(i = 1, j = 3:(length(table(df[, gvar]))+2), part = "header")
  }
  c <- flextable::fontsize(c, size = 10.5, part = "body") %>% 
    flextable::padding(j=1, i = which(grepl("  ", final[[1]])), padding.left = 20) %>% 
    flextable::font(fontname="Times New Roman",part="all") %>% 
    flextable::set_table_properties(layout = "autofit")
  print(c)
  if (out_docx==T){
    flextable::save_as_docx(c,path=path)
  }
}