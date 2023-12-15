library(foreign)
library(tidyverse)
library(kableExtra)
library(ltm)
library(foreign)
library(MatchIt)
## Internal Consistency
options(digits=2)
data <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_rev.sav",to.data.frame=TRUE)
data1 <- data.frame(id = (1:nrow(data)), Age = data$Age)
data1$Age <- factor(data1$Age)
rlne <- subset(data,select = grep("(rlne)\\d", names(data), value=TRUE))
snwe <- subset(data,select = grep("(snwe)\\d", names(data), value=TRUE))
sege <- subset(data,select = grep("(sege)\\d", names(data), value=TRUE))
lswe <- subset(data,select = grep("(lswe)\\d", names(data), value=TRUE))
rhme <- subset(data,select = grep("(rhme)\\d", names(data), value=TRUE))
lske <- subset(data,select = grep("(lske)\\d", names(data), value=TRUE))
df_list <-list(rlne,snwe,sege,lswe,rhme,lske)
alphas <- list()
for (i in 1:length(df_list)){
  get_alpha <- function(x) {  
      raw_alpha <-
          psych::alpha(df_list[[i]][data1[data1$Age == x, 1], ],check.keys = TRUE)$total[1]
          size <- nrow(df_list[[i]][data1[data1$Age == x, 1], ])
          SEM <- 15* sqrt(1-raw_alpha)
   data.frame(
          age = x,
          n = size,
          alpha = raw_alpha,
          SEM = SEM)
   }
  alphas[[i]] <- do.call(rbind,lapply(levels(data1$Age),get_alpha))
}
data_processor <- function (data,name){
  names(data) <- name
  df <- do.call("rbind",data)
  df[["scales"]] <- rep(names(data), sapply(data, nrow))
  rownames(df) <- NULL
  df_wide <- reshape(df,idvar=c("age","n"),v.names = c("raw_alpha","raw_alpha.1"),
                     timevar = "scales",direction="wide")
  names(df_wide) <- c("age","n","r","sem","r","sem","r","sem",
                    "r","sem","r","sem","r","sem")
  return(df_wide)
}
kbl(data_processor(alphas,c("rlne","snwe","sege","lswe","rhme","lske")),booktabs = T,
    caption = "Table 1") %>%
  kable_classic(latex_options="striped",full_width=F) %>%
  add_header_above(c(" " = 2 , "ERNL-E" = 2, "SPW-E" = 2, "ESEG-E" = 2,
                     "LSW-E" = 2, "RHY-E" = 2, "LSK-E" = 2),bold = T)

## Rasch IRT
data <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODS_Child.sav",to.data.frame=TRUE)
data1 <- data.frame(id = (1:nrow(data)), Age = data$Age)
data1$Age <-factor(data1$Age)
wrf <- subset(data,select = grep("WRF", names(data), value=TRUE))
qrf <- subset(data,select = grep("QRF", names(data), value=TRUE))
df_list <-list(wrf,qrf)
get_alphas <- list()
for (i in 1:length(df_list)){
  get_alpha <- function(x) {
    fit2 <- rasch(df_list[[i]][data1[data1$Age == x, 1], ])
    eap <- factor.scores(fit2, method="EAP")$score.dat
    e <- mean(eap$se.z1^2,na.rm=TRUE)
    s <- var(eap$z1,na.rm=TRUE)
    rasch_alpha <- 1-(e/(s+e))
    SEM <- 15 * sqrt(1-rasch_alpha)
    size <- nrow(df_list[[i]][data1[data1$Age  == x, 1], ])
    data.frame(
      age = as.character(x),
      n = size,
      rasch_alpha = rasch_alpha,
      sem = SEM
    )
  }
  get_alphas[[i]] <- do.call(rbind,lapply(levels(data1$Age),get_alpha))
}
data_processor <- function (data,name){
  names(data) <- name
  df <- do.call("rbind",data)
  df[["scales"]] <- rep(names(data), sapply(data, nrow))
  rownames(df) <- NULL
  df_wide <- reshape(df,idvar=c("age","n"),v.names = c("rasch_alpha","sem"),
                     timevar = "scales",direction="wide")
  names(df_wide) <- c("age","n","r","sem","r","sem")
  return(df_wide)
}
kbl(data_processor(get_alphas,c("wrf","qrf")),booktabs = T,caption = "Table 2") %>%
  kable_classic(latex_options="striped",full_width=F) %>%
  add_header_above(c(" " = 2 ,"WRF-S" = 2,"QRF-S" = 2),bold = T)

## Split-half Reliability 
data <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_rev.sav",to.data.frame=TRUE)
data_stand <- data[which(data$standsample == 1),]
data2 <- data.frame(id = (1:nrow(data_stand)), Age = data_stand$Age)
data2$Age <- factor(data2$Age)
rlne <- subset(data_stand,select = grep("(rlne)\\d", names(data_stand), value=TRUE))
snwe <- subset(data_stand,select = grep("(snwe)\\d", names(data_stand), value=TRUE))
snwe_new <- snwe[,c("snwe1","snwe2", "snwe3","snwe4","snwe5","snwe6","snwe11","snwe9",
                "snwe10","snwe7","snwe8","snwe13","snwe14","snwe12","snwe15","snwe20",
                "snwe18","snwe19","snwe24","snwe21","snwe16","snwe29","snwe17","snwe31",
                "snwe26","snwe30","snwe28","snwe35","snwe33","snwe34","snwe27","snwe32")] 
# sege raw data and reorder
sege <- subset(data_stand,select = grep("(sege)\\d", names(data_stand), value=TRUE))
sege_new <- sege[,c("sege1","sege2","sege3","sege4","sege5","sege6","sege7","sege9",
                    "sege8","sege10","sege13","sege11","sege12","sege16","sege15",
                    "sege14","sege17","sege18","sege19","sege20","sege21","sege23",
                    "sege22","sege25","sege24")]
# lswe raw data and reorder
lswe <- subset(data_stand,select = grep("(lswe)\\d", names(data_stand), value=TRUE))
lswe_new <- lswe[,c("lswe1","lswe2","lswe3","lswe4","lswe5","lswe6","lswe9","lswe7",
                    "lswe10","lswe8","lswe11","lswe14","lswe12","lswe13","lswe16",
                    "lswe18","lswe17","lswe19","lswe23","lswe21","lswe22","lswe25",
                    "lswe26","lswe31","lswe32","lswe30","lswe27","lswe33","lswe29",
                    "lswe34","lswe35","lswe28","lswe24","lswe38","lswe36","lswe37",
                    "lswe39","lswe40")]
# rhme raw data and reorder
rhme <- subset(data_stand,select = grep("(rhme)\\d", names(data_stand), value=TRUE))
rhme_new <- rhme[,c("rhme2","rhme5","rhme3","rhme4","rhme10","rhme1","rhme7","rhme8",
                    "rhme9","rhme6","rhme11","rhme16","rhme13","rhme18","rhme19",
                    "rhme14","rhme17","rhme15","rhme26","rhme22","rhme12","rhme20",
                    "rhme23","rhme25","rhme30","rhme21","rhme24","rhme28","rhme29",
                    "rhme27")]

# lske raw data and reorder
lske <- subset(data_stand,select = grep("(lske)\\d", names(data_stand), value=TRUE))
lske_new <- lske[,c("lske2","lske4", "lske1","lske5","lske3","lske6","lske8","lske7",
                    "lske9","lske10","lske12","lske13","lske14","lske15","lske16",
                    "lske19","lske18","lske20","lske17","lske23","lske22","lske21",
                    "lske24","lske26","lske27","lske30","lske32","lske31","lske28",
                    "lske29","lske33","lske34","lske35")] 
df_list <-list(rlne,snwe_new,sege_new,lswe_new,rhme_new,lske_new)
get_alphas <- list()
for (i in 1:length(df_list)){
  get_alpha <- function(x) {
    a <- df_list[[i]][data2[data2$Age == x, 1], ][,seq(1, ncol(df_list[[i]]), 2)]
    b <- df_list[[i]][data2[data2$Age == x, 1], ][,seq(2, ncol(df_list[[i]]), 2)]
    split <- list(a,b)
    split_sum <- lapply(split, function(x) data.frame(total = rowSums(x, na.rm=T)))
    split_alpha <- cor.test(split_sum[[1]]$total,split_sum[[2]]$total)$estimate
    full_split <- 2*split_alpha/(1+split_alpha)
    SEM <- 15 * sqrt(1-full_split)
    size <- nrow(df_list[[i]][data2[data2$Age == x, 1], ])
    data.frame(
      age = as.numeric(x),
      n = size,
      split_alpha = trunc(100 * full_split)/100,
      sem = SEM
    )
  }
  get_alphas[[i]] <- do.call(rbind,lapply(levels(data2$Age),get_alpha))
}
data_processor <- function (data,name){
  names(data) <- name
  df <- do.call("rbind",data)
  df[["scales"]] <- rep(names(data), sapply(data, nrow))
  rownames(df) <- NULL
  df_wide <- reshape(df,idvar=c("age","n"),v.names = c("split_alpha","sem"),
                     timevar = "scales",direction="wide")
  names(df_wide) <-  c("age","n","r","sem","r","sem","r","sem",
                    "r","sem","r","sem","r","sem")
  return(df_wide)
}
kbl(data_processor(get_alphas,c("rlne","snwe","sege","lswe","rhme","lske")),booktabs = T,
    caption = "Table 3") %>%
  kable_classic(latex_options="striped",full_width=F) %>%
  add_header_above(c(" " = 2 , "ERNL-E" = 2, "SPW-E" = 2, "ESEG-E" = 2,
                     "LSW-E" = 2, "RHY-E" = 2, "LSK-E" = 2),bold = T)
                        
## Linear Combination Reliability
library(foreign)
data <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_composites.sav",
                  to.data.frame=TRUE)
data$Grades <- ifelse((data$GradeSemester == 'K_Fall' | data$GradeSemester == 'K_Spring'), "K", 
                   ifelse((data$GradeSemester == '1_Fall' | data$GradeSemester == '1_Spring'), "1",
                   ifelse((data$GradeSemester == '2_Fall' | data$GradeSemester == '2_Spring'), "2",
                           0)))
data_stand <- data[which(data$standsample == 1),]
data1 <- data.frame(id = (1:nrow(data_stand)), Grade = data_stand$Grades)
data1$Grade <-factor(data1$Grade)
ELP_grade_ss <- subset(data_stand,select = grep("rhme_SS|rlne_SS|sege_SS", names(data_stand), value=TRUE))
df_list_ELP <- list(ELP_grade_ss)
alphas_ELP <- list()
ci <- list()
for (i in 1:length(df_list_ELP)){
get_alpha <- function(x) {
    var_cov <- cov(df_list_ELP[[i]][data1[data1$Grade == x, 1], ],use="everything")
    sum_c <- sum(var_cov,na.rm = T)
    diag_c <- diag(var_cov)
    if (any(x == "K")){
    raw_alpha <- 1 - (sum(diag_c,na.rm = T)-(sum(c(1,0.94,0.95)*diag_c,na.rm = T)))/sum_c
    }else if (any(x == 1)){
    raw_alpha <- 1 - (sum(diag_c,na.rm = T)-(sum(c(0.99,0.90,0.91)*diag_c,na.rm = T)))/sum_c  
    }else
    raw_alpha <- 1 - (sum(diag_c,na.rm = T)-(sum(c(0.99,0.83,0.91)*diag_c,na.rm = T)))/sum_c
    size <- nrow(df_list_ELP[[i]][data1[data1$Grade == x, 1], ])
    SEM <- 15* sqrt(1-raw_alpha)
  data.frame(
    grade = x,
    n = size,
    alpha = raw_alpha,
    SEM = SEM
  )
}
alphas_ELP[[i]] <- do.call(rbind,lapply(levels(data1$Grade)[c(3,1:2)],get_alpha))
}
                        
## Correlation Analysis   
data <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_final.sav",to.data.frame=TRUE)
final <- data[,c("rlne_SS","snwe_SS","sege_SS","lswe_SS","rhme_SS","lske_SS","pv_SS",
                 "lw_ss","qrf_ss","wrf_ss")]
names(final)<- c("ERNL-E","SPW-E","ESEG-E","LSW-E","RHY-E","LSK-E","PV-S","LWC-S","QRF-S","WRF-S")
corr_output <- cor(final,use = "pairwise.complete.obs")

## Effect Size
data_PT <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_PT.sav",to.data.frame=TRUE)
data1_PT <- data.frame(id = (1:nrow(data_PT)), Age = data_PT$Age)
data1_PT$Age <- factor(data1_PT$Age)
parents_new <- subset(data_PT,select = grep("(p)\\d|p.+tscore$", names(data_PT), value=TRUE))
teachers_new <- subset(data_PT,select = grep("(t)\\d|t.+tscore$", names(data_PT), value=TRUE))
parents_new[,1:(ncol(parents_new)-1)] <- ifelse(parents_new[,1:(ncol(parents_new)-1)] 
                                                == "Strongly Disagree", 1, 
                             ifelse(parents_new[,1:(ncol(parents_new)-1)] == "Disagree", 2,
                             ifelse(parents_new[,1:(ncol(parents_new)-1)] == "Agree", 3,
                             ifelse(parents_new[,1:(ncol(parents_new)-1)] == "Strongly Agree", 4, 0))))
teachers_new[,1:(ncol(teachers_new)-1)] <- ifelse(teachers_new[,1:(ncol(teachers_new)-1)] 
                                                  == "Strongly Disagree", 1, 
                               ifelse(teachers_new[,1:(ncol(teachers_new)-1)] == "Disagree", 2,
                               ifelse(teachers_new[,1:(ncol(teachers_new)-1)] == "Agree", 3,
                               ifelse(teachers_new[,1:(ncol(teachers_new)-1)] == "Strongly Agree", 4, 0))))
pt_list_new <- list(parents_new,teachers_new)
alphas_PT <- list()
for (i in 1:length(pt_list_new)){
get_alpha_PT <- function(x) {  
  raw_alpha_PT <-
    round(psych::alpha(pt_list_new[[i]][,1:ncol(pt_list_new[[i]])-1][data1_PT[data1_PT$Age == x, 1], ],
                       check.keys = TRUE)$total[1],digits = 2)
    size <- sum(!is.na(pt_list_new[[i]][,ncol(pt_list_new[[i]])][data1_PT[data1_PT$Age == x, 1]]))
    overlap <- as.data.frame(cbind(pt_list_new[[1]][,ncol(pt_list_new[[1]])][data1_PT[data1_PT$Age == x, 1]],
                               pt_list_new[[2]][,ncol(pt_list_new[[2]])][data1_PT[data1_PT$Age == x, 1]]))
    overlap_new <- overlap[!with(overlap,is.na(V1)|is.na(V2)),]
    overlap_size <- nrow(overlap_new)
    t_scores <- pt_list_new[[i]][,ncol(pt_list_new[[i]])][data1_PT[data1_PT$Age == x, 1]]
    t_avg <- round(mean(t_scores,na.rm = TRUE),digits = 2)
    t_sd <- round(sd(t_scores,na.rm = TRUE),digits = 1) 
  data.frame(
    age = as.numeric(x),
    n = size,
    alpha = raw_alpha_PT,
    mean = t_avg,
    sd = t_sd,
    n_overlap = overlap_size)
  }
alphas_PT[[i]] <- do.call(rbind,lapply(levels(data1_PT$Age),get_alpha_PT))
}
d <- round(abs((alphas_PT[[1]]["mean"]-alphas_PT[[2]]["mean"])/sqrt(((alphas_PT[[1]]["sd"])^2 
            + (alphas_PT[[2]]["sd"])^2)/2)), digits = 2)
alphas_es <- Map(cbind, alphas_PT , effect_size = d)
data_processor <- function (data,name){
  names(data) <- name
  df <- do.call("rbind",data)
  df[["scales"]] <- rep(names(data), sapply(data, nrow))
  rownames(df) <- NULL
  df_wide <- reshape(df,idvar=c("age","n_overlap"),v.names = c("n", "raw_alpha","mean",
                                "sd","effect_size"),timevar = "scales",direction="wide")[,-7]
  names(df_wide) <- c("age","n","alpha","mean","sd","n","alpha","mean",
                                  "sd","n_overlap","effect size")
  return(df_wide)
}
kbl(data_processor(alphas_es,c("parents","teachers")),booktabs = T,caption = "Table 4") %>%
  kable_classic(latex_options="striped",full_width=F) %>%
  add_header_above(c(" " = 1 , "parents" = 4, "teachers" = 4," " = 1," " = 1),bold = T)

## Matched Sample
data_final <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_final.sav",
                        to.data.frame=TRUE)
m.out1 <- matchit(Face_Masks ~ Age + Gender + Ethnicity +  SES, data = data_final,
                  method = "nearest", distance = "glm")
matched.TODE <- match.data(m.out1)
tblFun <- function(x){
  tbl <- table(x)
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('n','percent(%)')
  res
}
data_split <- split(matched.TODE, matched.TODE$Face_Masks)
result <- list()
final <- list()
for (i in 1:length(data_split)){
  result[[i]] <- lapply(data_split[[i]][c("Gender","SES" ,"Ethnicity","Age")], function(x) tblFun(x))
  final[[i]] <- do.call("rbind", result[[i]])
  mnm <- do.call(cbind,final)
}
kbl(mnm, caption = "Table 5") %>%
  #kable_paper("striped", full_width = F) %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  add_header_above(c(" " = 1, "non_mask" = 2, "mask" = 2))%>%
  pack_rows("Gender", 1, 3) %>%
  pack_rows("SES", 4, 7) %>%
  pack_rows("Ethnicity", 8, 14) %>%
  pack_rows("Age", 15, 18)
                        
## Frequency Tables
data_final <- read.spss("https://raw.github.com/wpscode/TOD/master/TOD%20Data/TODE_final.sav",
                        to.data.frame=TRUE) 
data_stand <- subset(data,standsample == 1)
tblFun <- function(x){
  tbl <- table(x)
  res <- cbind(tbl,round(prop.table(tbl)*100,2))
  colnames(res) <- c('n','percent(%)')
  res
}
ft <- do.call("rbind",lapply(data_stand[c("Gender","SES" ,"Ethnicity","region",
                                        "Age","GradeSemester")],function(x) tblFun(x)))
ft_new <- ft[!(apply(ft, 1, function(y) any(y == 0))),]
cd <- data.frame(census = c(format(round(51.0, 1), nsmall = 1),format(round(49.0, 1), nsmall = 1),
                            11.5,26.1,30.3,32.2,4.8,13.5,50.2,0.7,0.2,5.2,
                            25.4,15.8,21.2,38.6,24.4,"","","","","","","","","",""))
names(cd) <- c("census(%)")
rownames(cd) <- rownames(ft_new)
sample_nation <- cbind(ft_new,cd)
kbl(sample_nation, caption = "Frequency table") %>%
  kable_classic(full_width = F, html_font = "Cambria")%>%
  pack_rows("Gender", 1, 2) %>%
  pack_rows("SES", 3, 6) %>%
  pack_rows("Ethnicity", 7, 13) %>%
  pack_rows("Region", 14, 17) %>%
  pack_rows("Age", 18, 21) %>%
  pack_rows("Grade", 22, 27)
