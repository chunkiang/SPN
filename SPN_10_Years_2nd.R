################### 导入数据，加载数据####
setwd("~/R/SPN_10_Years")
library("foreign", lib.loc="C:/Program Files/R/R-3.3.2/library")
SPN_10_Years<- read.dbf(file="D:/WHONET5/Output/SPN_10_Years_R_2nd.dbf")
mic <- SPN_10_Years[30:45]
str(mic)
names(mic)

##去除所有的非数字比如<= and >
##方法2：不需要任何修改，whonet的原始数据就能用
##方法4：先改小于等于号，再改所有的符号，最后整理为数字格式
for(j in 1:ncol(mic)){
  b <- NULL; a<-NULL
  for(i in 1:nrow(mic)){
    #b <- c(b,stringi::stri_replace_all_fixed(mic[i,j],c("<=."),c("0."), vectorize_all=FALSE))
    a <- stringi::stri_replace_all_fixed(mic[i,j],c("."),c("0."), vectorize_all=FALSE)
    a <- stringi::stri_replace_all_fixed(a,c("<="),c(""), vectorize_all=FALSE)
    b <- c(b,stringi::stri_trim_both(a, "\\p{N}"))
  }
  mic[,j] <- as.numeric(b);rm(a,b,i,j)
}


#SPN通用折点设置  16个药 可以判断是否是空值
#SPN通用折点设置  16个药 可以判断是否是空值使用前需要先格式化数据，数据只能是数字
breakpoint <- function(drug, mic){
  if(is.na(mic)==F){
    if(drug == "PEN"|drug == "PEN_NM"){
      if(mic <= 0.064)  return("S")
      if(mic >0.064 & mic < 2)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "AMC"|drug == "AMC_NM"){
      if(mic <= 2)       return("S")
      if(mic == 4)       return("I")
      else{
        return("R")
      }
    }
    #if(drug == "FEP"|drug == "FEP_NM"){#非脑膜炎折点
    #  if(mic <= 1)       return("S")
    #  if(mic == 2)       return("I")
    #  else{
    #    return("R")
    #  }
    #}
    if(drug == "CRO"|drug == "CRO_NM"){
      if(mic <= 1)  return("S")
      if(mic == 2)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "CXM"|drug == "CEC"|drug == "CXM_NM"|drug == "CEC_NM"){
      if(mic <= 1)  return("S")
      if(mic == 2)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "VAN"|drug == "VAN_NM"){
      if(mic <= 1)  return("S")
      #if(mic >= 4)       return("R")
      else{
        return(NA)#spn对这个药不可能出现耐药的情况，如果mic特别高则自动返回缺失值，保证敏感率计算任然为100%
      }
    }
    if(drug == "ERY"|drug == "ERY_NM"){
      if(mic <= 0.25)   return("S")
      if(mic >0.25 & mic < 1)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "AZM"|drug == "AZM_NM"){
      if(mic <= 0.5)    return("S")
      if(mic >0.5 & mic < 2)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "CLA"|drug == "CLA_NM"){
      if(mic <= 0.25)  return("S")
      if(mic >0.25 & mic < 1)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "TCY"|drug == "TCY_NM"){
      if(mic <= 1)  return("S")
      if(mic == 2)       return("I")
      else{
        return("S")
      }
    }
    if(drug == "LVX"|drug == "LVX_NM"){
      if(mic <= 2)  return("S")
      if(mic == 4)  return("I")
      else{
        return("R")
      }
    }
    if(drug == "MFX"|drug == "MFX_NM"){
      if(mic <= 1)  return("S")
      if(mic == 2)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "SXT"|drug == "SXT_NM"){
      if(mic <= 0.5)  return("S")
      if(mic >0.5 & mic < 4)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "CHL"|drug == "CHL_NM"){
      if(mic <= 4)       return("S")
      else{
        return("R")
      }
    }
    if(drug == "CLI"|drug == "CLI_NM"){
      if(mic <= 0.25)  return("S")
      if(mic >0.25 & mic < 1)       return("I")
      else{
        return("R")
      }
    }
    if(drug == "LNZ"|drug == "LNZ_NM"){
      if(mic <= 2)  return("S")
      #if(mic >= 4)       return("R")
      else{
        return(NA)#spn对这个药不可能出现耐药的情况，如果mic特别高则自动返回缺失值，保证敏感率计算任然为100%
      }
    }
    else{
      return(NA)
    } ##如果不是spn的药物则返回缺失值
  }else{
    return(NA)
  }#如果没有mic数值则返回缺失值
}

##将mic转换为ris，产生一个以 药物名称_sir 的list，里面全是sir。
##并且产生drug_name，存放所有 药物名称_sir list的名称。使用前一定要用str()看一下数据，必须全部都是num才能用
##产生一个新的data名称为sir，里面放置的是每种药物的sir。
##同时产生value list，名称为 药物名称_sir。
##使用前一定要用str()看一下数据，必须全部都是num才能用
##如果使用了whonet数据清扫除去所有的“<=”,">",".0"则可以直接使用
drug.name.list <-NULL; sir <- data.frame(matrix(ncol = ncol(mic), nrow = nrow(mic) ))	#产生空的药物名称变量，名称为sir的data，行列为mic data的行列
for(j in 1:ncol(mic)){																	#从第一列到最后一列
  interpretation <- NULL;																#interpretation 就是 “S”、“I”、“R”
  for(i in 1:nrow(mic)){																#从给定列的第一行开始，到最后一行
    interpretation <- c(interpretation, breakpoint(names(mic[j]), mic[i,j]))			#从给定列的第一行开始，到最后一行运行结束后，形成一个SSSRRRIISS的list
    name <- paste(names(mic[j]),"_SIR",sep = "")										#以此列的变量名称生产新的变量名称，名称后缀有“_SIR”,新变量名称为 name
    assign(name,as.factor(interpretation))												#把生成的sir list 赋值给新生成的变量
    sir[i,j] <- breakpoint(names(mic[j]), mic[i,j])										#把breakpoint函数的返回值传递给sir表的指定位置
    #print(name); print(i)   															#这个是监测哪里出错时用的，可以返回处理哪个变量，在第几行出错
  }
  drug.name.list <- c(drug.name.list,name)												#生成所有药物的list
  #print(drug.name.list)																	#打印已处理的药物的列表，视觉效果很cool
  sir[,j]<- as.factor(sir[,j])															#把sir表中的字符型变量转变为factor后，就可以用table()统计了
  rm(i,j,interpretation,name)																	#删除中间变量
}
colnames(sir) <- drug.name.list															#给sir表的变量名称赋值												


##生成mic结果和sir结果肩并肩的显示
mic.sir <- data.frame(matrix(ncol = 2*ncol(mic), nrow = nrow(mic) ))
mic.sir.name <- NULL
mic.name <- names(mic)
sir.name <- names(sir)
for(j in 1:ncol(mic)){
  for(i in 1:nrow(mic)){
    mic.sir[i,(j*2-1)] <- (mic[i,j])
    mic.sir[i,(j*2)] <- as.character(sir[i,j])
  }
  mic.sir[,(j*2)] <- as.factor(mic.sir[,(j*2)])
  mic.sir.name <- c(mic.sir.name,mic.name[j],sir.name[j])
}
colnames(mic.sir) <- mic.sir.name



##这个是生成一个第一列是抗生素，第二列是sir，第三列是mic的表，这个表可以用table直接计算sir的比例
dim(sir)
drug.sir<-data.frame(matrix(ncol = 3, nrow = nrow(mic)*(ncol(mic)-0 ) )); colnames(drug.sir) <- c("drug","sir","mic") 
#如果第一列是菌株号，那么列数要减一

#把抗菌素名称列入第一例
drug.name.list <- NULL
for(j in 1:ncol(sir)){
  for(i in 1:nrow(sir)){
    drug.name.list <- c(drug.name.list,names(sir[j]))
  }
}
drug.sir[,1] <- as.factor(drug.name.list)


#把sir放入第二列
drug.sir.list <- NULL
for(j in 1:ncol(sir)){
  for(i in 1:nrow(sir)){
    drug.sir.list <- c(drug.sir.list,as.character(sir[i,j]))
  }
}
drug.sir.list <- as.factor(drug.sir.list)
levels(drug.sir.list)
drug.sir.list <- factor(drug.sir.list, levels(drug.sir.list)[c(3,1,2)])## reordering the "RIS"
drug.sir[,2] <- as.factor(drug.sir.list)


#把mic放入第三列
drug.mic.list <- NULL
for(j in 1:ncol(mic)){
  for(i in 1:nrow(mic)){
    drug.mic.list <- c(drug.mic.list,mic[i,j])
  }
}
drug.sir[,3] <- as.numeric(drug.mic.list)

##加入年份变量
SPN_10_Years$SPEC_DATE[1:100]#to inspect the format of year.
levels(as.factor(SPN_10_Years$SPEC_DATE))
year.by1 <- as.factor(substr(SPN_10_Years$SPEC_DATE,1,4))
drug.sir$year <- rep(year.by1, time =16)
##分析不同药物的敏感性
table(drug.sir$drug,drug.sir$sir)
addmargins(prop.table(table(drug.sir$drug,drug.sir$sir),1))*100

sircol <- c(
  "R"="#f46d43",
  "I"="#74add1",
  "S"="#4575b4"
)


ggplot(subset(drug.sir, !is.na(sir)), aes(x = drug, fill = sir))+ ##去除sir缺失值
  #facet_wrap(~drug)+
  geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="S_I_R" )


ggplot(subset(drug.sir, !is.na(sir)), aes(x = drug, fill=sir))+
  geom_bar(position="fill")



##不同药物在不同年份耐药率的变化
drug.vs.years <- ggplot(subset(drug.sir, !is.na(sir)&!is.na(year)), aes(x = year, fill=sir))+
  facet_wrap(~drug)+
  geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="Rsistance" )


##不同药物在不同年份耐药率的变化 五岁以下儿童
drug.vs.yearsless5 <- ggplot(subset(drug.sir[drug.sir$age<=5,], !is.na(sir)&!is.na(year)), aes(x = year, fill=sir))+
    facet_wrap(~drug)+
    geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="Rsistance" )
drug.vs.yearsless5
