################### �������ݣ���������####
setwd("~/R/SPN_10_Years")
library("foreign", lib.loc="C:/Program Files/R/R-3.3.2/library")
SPN_10_Years<- read.dbf(file="D:/WHONET5/Output/SPN_10_Years_R_2nd.dbf")
mic <- SPN_10_Years[30:45]
str(mic)
names(mic)

##ȥ�����еķ����ֱ���<= and >
##����2������Ҫ�κ��޸ģ�whonet��ԭʼ���ݾ�����
##����4���ȸ�С�ڵ��ںţ��ٸ����еķ��ţ��������Ϊ���ָ�ʽ
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


#SPNͨ���۵�����  16��ҩ �����ж��Ƿ��ǿ�ֵ
#SPNͨ���۵�����  16��ҩ �����ж��Ƿ��ǿ�ֵʹ��ǰ��Ҫ�ȸ�ʽ�����ݣ�����ֻ��������
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
    #if(drug == "FEP"|drug == "FEP_NM"){#����Ĥ���۵�
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
        return(NA)#spn�����ҩ�����ܳ�����ҩ����������mic�ر�����Զ�����ȱʧֵ����֤�����ʼ�����ȻΪ100%
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
        return(NA)#spn�����ҩ�����ܳ�����ҩ����������mic�ر�����Զ�����ȱʧֵ����֤�����ʼ�����ȻΪ100%
      }
    }
    else{
      return(NA)
    } ##�������spn��ҩ���򷵻�ȱʧֵ
  }else{
    return(NA)
  }#���û��mic��ֵ�򷵻�ȱʧֵ
}

##��micת��Ϊris������һ���� ҩ������_sir ��list������ȫ��sir��
##���Ҳ���drug_name��������� ҩ������_sir list�����ơ�ʹ��ǰһ��Ҫ��str()��һ�����ݣ�����ȫ������num������
##����һ���µ�data����Ϊsir��������õ���ÿ��ҩ���sir��
##ͬʱ����value list������Ϊ ҩ������_sir��
##ʹ��ǰһ��Ҫ��str()��һ�����ݣ�����ȫ������num������
##���ʹ����whonet������ɨ��ȥ���еġ�<=��,">",".0"�����ֱ��ʹ��
drug.name.list <-NULL; sir <- data.frame(matrix(ncol = ncol(mic), nrow = nrow(mic) ))	#�����յ�ҩ�����Ʊ���������Ϊsir��data������Ϊmic data������
for(j in 1:ncol(mic)){																	#�ӵ�һ�е����һ��
  interpretation <- NULL;																#interpretation ���� ��S������I������R��
  for(i in 1:nrow(mic)){																#�Ӹ����еĵ�һ�п�ʼ�������һ��
    interpretation <- c(interpretation, breakpoint(names(mic[j]), mic[i,j]))			#�Ӹ����еĵ�һ�п�ʼ�������һ�����н������γ�һ��SSSRRRIISS��list
    name <- paste(names(mic[j]),"_SIR",sep = "")										#�Դ��еı������������µı������ƣ����ƺ�׺�С�_SIR��,�±�������Ϊ name
    assign(name,as.factor(interpretation))												#�����ɵ�sir list ��ֵ�������ɵı���
    sir[i,j] <- breakpoint(names(mic[j]), mic[i,j])										#��breakpoint�����ķ���ֵ���ݸ�sir����ָ��λ��
    #print(name); print(i)   															#����Ǽ���������ʱ�õģ����Է��ش����ĸ��������ڵڼ��г���
  }
  drug.name.list <- c(drug.name.list,name)												#��������ҩ���list
  #print(drug.name.list)																	#��ӡ�Ѵ�����ҩ����б����Ӿ�Ч����cool
  sir[,j]<- as.factor(sir[,j])															#��sir���е��ַ��ͱ���ת��Ϊfactor�󣬾Ϳ�����table()ͳ����
  rm(i,j,interpretation,name)																	#ɾ���м����
}
colnames(sir) <- drug.name.list															#��sir���ı������Ƹ�ֵ												


##����mic�����sir����粢�����ʾ
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



##���������һ����һ���ǿ����أ��ڶ�����sir����������mic�ı��������������tableֱ�Ӽ���sir�ı���
dim(sir)
drug.sir<-data.frame(matrix(ncol = 3, nrow = nrow(mic)*(ncol(mic)-0 ) )); colnames(drug.sir) <- c("drug","sir","mic") 
#�����һ���Ǿ���ţ���ô����Ҫ��һ

#�ѿ��������������һ��
drug.name.list <- NULL
for(j in 1:ncol(sir)){
  for(i in 1:nrow(sir)){
    drug.name.list <- c(drug.name.list,names(sir[j]))
  }
}
drug.sir[,1] <- as.factor(drug.name.list)


#��sir����ڶ���
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


#��mic���������
drug.mic.list <- NULL
for(j in 1:ncol(mic)){
  for(i in 1:nrow(mic)){
    drug.mic.list <- c(drug.mic.list,mic[i,j])
  }
}
drug.sir[,3] <- as.numeric(drug.mic.list)

##������ݱ���
SPN_10_Years$SPEC_DATE[1:100]#to inspect the format of year.
levels(as.factor(SPN_10_Years$SPEC_DATE))
year.by1 <- as.factor(substr(SPN_10_Years$SPEC_DATE,1,4))
drug.sir$year <- rep(year.by1, time =16)
##������ͬҩ���������
table(drug.sir$drug,drug.sir$sir)
addmargins(prop.table(table(drug.sir$drug,drug.sir$sir),1))*100

sircol <- c(
  "R"="#f46d43",
  "I"="#74add1",
  "S"="#4575b4"
)


ggplot(subset(drug.sir, !is.na(sir)), aes(x = drug, fill = sir))+ ##ȥ��sirȱʧֵ
  #facet_wrap(~drug)+
  geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="S_I_R" )


ggplot(subset(drug.sir, !is.na(sir)), aes(x = drug, fill=sir))+
  geom_bar(position="fill")



##��ͬҩ���ڲ�ͬ�����ҩ�ʵı仯
drug.vs.years <- ggplot(subset(drug.sir, !is.na(sir)&!is.na(year)), aes(x = year, fill=sir))+
  facet_wrap(~drug)+
  geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="Rsistance" )


##��ͬҩ���ڲ�ͬ�����ҩ�ʵı仯 �������¶�ͯ
drug.vs.yearsless5 <- ggplot(subset(drug.sir[drug.sir$age<=5,], !is.na(sir)&!is.na(year)), aes(x = year, fill=sir))+
    facet_wrap(~drug)+
    geom_bar(position="fill")+
  scale_fill_manual(values = sircol,
                    name="Rsistance" )
drug.vs.yearsless5