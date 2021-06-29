ReadData <- function(path){
  data <- read.csv(path)
  data <- data[data$CycID != "#N/A",]
  return(data)
}
GetDataFeatures <- function(expr){
  expr_mean <- log(expr$meta2d_Base)
  expr_amp <- log(expr$meta2d_AMP)
  expr_ramp <- log(expr$meta2d_rAMP)
  ind_bhq005 <- expr$meta2d_BH.Q<0.05
  ind_cir = 'ind_bhq005'
  aa0 <- cbind(expr_mean,expr_ramp,ind_bhq005)
  aa <-aa0[order(aa0[,'expr_ramp'],decreasing = T),];thresh_ind = 
    which(aa[,ind_cir]==1)[round(sum(aa[,ind_cir])*0.95)];aa = aa[1:thresh_ind,]
  
  table1 = c(nrow(aa0),nrow(aa),round(nrow(aa)/nrow(aa0),2),
             sum(aa0[,ind_cir]),sum(aa[,ind_cir]),round(sum(aa[,ind_cir])/sum(aa0[,ind_cir]),2))
  m = sum(aa0[,ind_cir])
  n = nrow(aa0) - sum(aa0[,ind_cir])
  k = sum(aa0[,ind_cir])
  aa = aa[order(aa[,'expr_mean']),]
  dd = exp(aa[,'expr_mean'])
  ee = replicate(1000,{
    # print(a <<- a+1)
    ind_simu = sample(1:nrow(aa),size = sum(aa[,ind_cir]),prob = dd)
    sum(aa[ind_simu,ind_cir])
  })
  dhydata = dhyper(0:m,m,n,k)
  v=m^2/(m+n)
  return_list <- list(ee, table1, dhydata,v)
  return(return_list)
}




GetDataFeatures <- function(expr){
  expr_mean <- log(expr$meta2d_Base)
  expr_amp <- log(expr$meta2d_AMP)
  expr_ramp <- log(expr$meta2d_rAMP)
  ind_bhq005 <- expr$meta2d_BH.Q<0.05
  ind_cir = 'ind_bhq005'
  aa0 <- cbind(expr_mean,expr_ramp,ind_bhq005)
  aa <-aa0[order(aa0[,'expr_ramp'],decreasing = T),];thresh_ind = 
    which(aa[,ind_cir]==1)[round(sum(aa[,ind_cir])*0.95)];aa = aa[1:thresh_ind,]
  
  table1 = c(nrow(aa0),nrow(aa),round(nrow(aa)/nrow(aa0),2),
             sum(aa0[,ind_cir]),sum(aa[,ind_cir]),round(sum(aa[,ind_cir])/sum(aa0[,ind_cir]),2))
  m = sum(aa0[,ind_cir])
  n = nrow(aa0) - sum(aa0[,ind_cir])
  k = sum(aa0[,ind_cir])
  aa = aa[order(aa[,'expr_mean']),]
  dd = exp(aa[,'expr_mean'])
  ee = replicate(1000,{
    # print(a <<- a+1)
    ind_simu = sample(1:nrow(aa),size = sum(aa[,ind_cir]),prob = dd)
    sum(aa[ind_simu,ind_cir])
  })
  dhydata = dhyper(0:m,m,n,k)
  v=m^2/(m+n)
  return_list <- list(ee, table1, dhydata,v)
  return(return_list)
}




GetDataFeaturesRange <- function(expr){
  expr_max <- apply(expr[,24:35], 1, max)
  expr_min <- apply(expr[,24:35], 1, min)
  expr_range <- expr_max - expr_min
  expr_mean <- log(expr$meta2d_Base)
  expr_amp <- log(expr$meta2d_AMP)
  expr_ramp <- log(expr$meta2d_rAMP)
  ind_bhq005 <- expr$meta2d_BH.Q<0.05
  ind_cir = 'ind_bhq005'
  aa0 <- cbind(expr_mean,expr_range,ind_bhq005)
  aa <-aa0[order(aa0[,'expr_range'],decreasing = T),];thresh_ind = 
    which(aa[,ind_cir]==1)[round(sum(aa[,ind_cir])*0.95)];aa = aa[1:thresh_ind,]
  
  table1 = c(nrow(aa0),nrow(aa),round(nrow(aa)/nrow(aa0),2),
             sum(aa0[,ind_cir]),sum(aa[,ind_cir]),round(sum(aa[,ind_cir])/sum(aa0[,ind_cir]),2))
  m = sum(aa0[,ind_cir])
  n = nrow(aa0) - sum(aa0[,ind_cir])
  k = sum(aa0[,ind_cir])
  aa = aa[order(aa[,'expr_range']),]
  dd = exp(aa[,'expr_range'])
  ee = replicate(1000,{
    # print(a <<- a+1)
    ind_simu = sample(1:nrow(aa),size = sum(aa[,ind_cir]),prob = dd)
    sum(aa[ind_simu,ind_cir])
  })
  dhydata = dhyper(0:m,m,n,k)
  v=m^2/(m+n)
  return_list <- list(ee, table1, dhydata,v)
  return(return_list)
}



f_completion = function(data,method='median'){
  id_na = is.na(data)
  if(method=='median'){data[id_na] = median(data,na.rm = T)}
  return(data)
}





GetRnaInfo <- function(species, codingseq, transcript_length, parameters,data){
  tissuename= c("Adr","Aorta","BFAT","BS","Cere","Heart",
                "Hypo","Kidney","Liver","Lung","Mus","WFAT")
  tissuenamefull= c("Adrenal Gland","Aorta","Brown fat","Brainstem","Cerebellum","Heart",
                    "Hypothalamus","Kidney","Liver","Lung","Skeletal muscle","White fat")
  RNAtype=c("5utr","3utr","fulllength","coding","intron","protein")
  if(species == "yest"){
    v_yeast = 69
    cmain_yeast = 18.785*10^9
    cgrow_yeast = 2468.2*10^9
    ntotr_yeast = 60*10^3
    ntotp_yeast = 82.2*10^6
    Ne_yeast = 10^7
    ntotr = ntotr_yeast
    ntotp = ntotp_yeast
    Ne = Ne_yeast
  }else{
    v_mouse = 2040
    cgrow_mouse = -10^9
    cmain_mouse = -1^9
    ntotr_mouse = 173000
    ntotp_mouse = 3.133*10^9
    Ne_mouse = c(25000,120000)
    ntotr = ntotr_mouse
    ntotp = ntotp_mouse
    Ne = Ne_mouse
  }
  
  # 基本符号CC=sumCM，LPRE，D=ln2/HL，R=expr，F ---------------------
  # RNA --------------------------------------------------------------
  # mRNA的序列与长度，sequence_list_DNA，sequence_list_name_DNA，LMAT_RNA
  con = file(codingseq, "r")
  sequence_list_DNA = c()
  sequence_list_name_DNA = c()
  for(i in 1:13384){
    #print(i)
    line= readLines(con,1)
    sequence_list_name_DNA = c(sequence_list_name_DNA,line)
    line= readLines(con,1)
    sequence_list_DNA = c(sequence_list_DNA, line)
  }
  close(con)
  sequence_list_name_DNA = substr(sequence_list_name_DNA,2,nchar(sequence_list_name_DNA))
  LMAT_RNA = nchar(sequence_list_DNA)
  names(LMAT_RNA) = sequence_list_name_DNA
  # CC_RNA，向量 -------------------------------------------------
  CM_RNA = c(49.7,46.7,50.7,45.7)
  names(CM_RNA) = c('A','C','G','T')
  CC_RNA = c()
  for (i in 1:length(sequence_list_DNA)) {
    tmp = table(strsplit(sequence_list_DNA[i],''))
    CC_RNA = c(CC_RNA,sum(tmp*CM_RNA[names(tmp)]))
  }
  names(CC_RNA) = sequence_list_name_DNA
  # LPRE_RNA，向量 ----------------------------------------------------
  tmp = read.table(transcript_length,
                   header = T,sep = '\t')
  LPRE_RNA = tmp$Transcript.end..bp.-tmp$Transcript.start..bp.+1
  names(LPRE_RNA) = tmp$Gene.stable.ID
  # D_RNA，向量 ---------------------------------------------------------
  tmp = read.table(parameters,header = T)
  HL_RNA = tmp$mRNA_halflive
  D_RNA = log(2)/HL_RNA
  names(D_RNA) = tmp$ID
  # R_RNA，矩阵 ---------------------------------------------------
  tmp = read.csv(data)
  R_RNA = tmp[,24:43]
  rownames(R_RNA) = tmp$CycID
  # 去掉含有0的行
  R_RNA = R_RNA[rowSums(R_RNA == 0) == 0,]
  # 统一id顺序，genename ----------------------------------------------
  # union
  genename = Reduce(union,list(v1 = names(CC_RNA),v2 = names(D_RNA),v3 = names(LPRE_RNA),v4 = rownames(R_RNA)))
  # 排序
  genename = as.character(sort(genename))
  # 其他符号都按这个排
  CC_RNA = CC_RNA[genename]
  D_RNA = D_RNA[genename]
  LPRE_RNA = LPRE_RNA[genename]
  R_RNA = R_RNA[genename,]
  # 赋name
  names(CC_RNA) = genename
  names(D_RNA) = genename
  names(LPRE_RNA) = genename
  rownames(R_RNA) = genename
  # 补全NA
  CC_RNA = f_completion(CC_RNA) 
  D_RNA = f_completion(D_RNA) 
  LPRE_RNA = f_completion(LPRE_RNA) 
  R_RNA_rowname = rownames(R_RNA)
  R_RNA_colname = colnames(R_RNA)
  for (i in 1:ncol(R_RNA)) {
    R_RNA[,i] = f_completion(R_RNA[,i])
  }
  # F_RNA，矩阵 ------------------------------------------------------------------
  # F_S_RNA
  F_S_RNA = R_RNA
  for (i in 1:ncol(R_RNA)) {
    F_S_RNA[i] = R_RNA[i]*D_RNA
  }
  
  result_list <- list(F_S_RNA, CC_RNA, LPRE_RNA, R_RNA, genename)
  
  return(result_list)
}




GetProteinInfo <- function(species, protein_file, parameters, mRNA_expression, protein_abundance, R_RNA, genename){
  # PRO --------------------------------------------------------------
  # PRO的序列与长度，sequence_list_name_PRO，sequence_list_PRO，LMAT_PRO
  if(species == "yest"){
    v_yeast = 69
    cmain_yeast = 18.785*10^9
    cgrow_yeast = 2468.2*10^9
    ntotr_yeast = 60*10^3
    ntotp_yeast = 82.2*10^6
    Ne_yeast = 10^7
    ntotr = ntotr_yeast
    ntotp = ntotp_yeast
    Ne = Ne_yeast
  }else{
    v_mouse = 2040
    cgrow_mouse = -10^9
    cmain_mouse = -1^9
    ntotr_mouse = 173000
    ntotp_mouse = 3.133*10^9
    Ne_mouse = c(25000,120000)
    ntotr = ntotr_mouse
    ntotp = ntotp_mouse
    Ne = Ne_mouse
  }
  con = file(protein_file, "r")
  sequence_list_name_PRO = c()
  sequence_list_PRO = c()
  for(i in 1:13384){
    #print(i)
    line= readLines(con,1)
    sequence_list_name_PRO = c(sequence_list_name_PRO,line)
    line= readLines(con,1)
    sequence_list_PRO = c(sequence_list_PRO, line)
  }
  close(con)
  sequence_list_name_PRO = substr(sequence_list_name_PRO,2,nchar(sequence_list_name_PRO))
  # 去掉末尾*号
  sequence_list_PRO = substr(sequence_list_PRO,1,nchar(sequence_list_PRO)-1)
  LMAT_PRO = nchar(sequence_list_PRO)
  names(LMAT_PRO) = sequence_list_name_PRO
  # CC_PRO，向量 ---------------------------------------------------------
  CM_PRO = c(23.55,10,23,13,11,23,11,12,10,33,29,22,27,31,46,16,10,17,63,44,20)
  names(CM_PRO) = c('X','A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')
  CC_PRO = c()
  for (i in 1:length(sequence_list_PRO)) {
    tmp = table(strsplit(sequence_list_PRO[i],''))
    CC_PRO = c(CC_PRO,sum(tmp*CM_PRO[names(tmp)]))
  }
  names(CC_PRO) = sequence_list_name_PRO
  # LPRE_PRO，向量 ----------------------------------------------
  LPRE_PRO = LMAT_PRO
  # D_PRO，向量 ----------------------------------------------------
  tmp = read.table(parameters,header = T)
  HL_PRO = tmp$Protein_halflive
  D_PRO = log(2)/HL_PRO
  names(D_PRO) = tmp$ID
  # R_PRO，矩阵 -------------------------------------------------------------
  # 以下计算方式二选一：若所依据实验数据中RNA与PRO相对关系可信则用方式三，否则用方式四
  # 方式三，通过一次实验中的数据得到R_RNA与R_PRO的比例关系，并用到全体时刻中
  tmp = read.table(mRNA_expression,header = T)
  tmp1 = tmp$mRNA_Transcription_level
  names(tmp1) = tmp$Name
  tmp = read.table(protein_abundance,header = T)
  tmp2 = tmp$Protein_abundance
  names(tmp2) = toupper(tmp$ORF)
  tmp_id = names(tmp1)
  tmp = tmp2[tmp_id]/tmp1[tmp_id]
  tmp = na.omit(tmp) # 得到比例tmp，然后与R_RNA相乘
  R_PRO = R_RNA
  for (i in 1:ncol(R_RNA)) {
    # 因为R_RNA是按照genename排的序，所以直接使用genename作为标号即可
    R_PRO[i] = R_RNA[i]*tmp[genename]
  }
  for (i in 1:ncol(R_PRO)) {
    R_PRO[,i] = f_completion(R_PRO[,i])
  }
  # 方式四，在方式三的基础上再限制sumR_RNA与sumR_PRO的比例
  # yeast中暂用该方式
  # sumR_RNA/sumR_PRO = ntotr/ntotp
  # 由于R_RNA与R_PRO存在NA和20等数据缺失问题，它们不再代表全体，所以该限制应在tmp1与tmp2阶段做出
  # 即sumtmp1/sumtmp2 = ntotr/ntotp
  # 此处两种方法，一种在tmp1与tmp2对齐之前保持比例，一种在对其并去NA之后保持比例。暂采用后一种
  tmp = read.table(mRNA_expression,header = T)
  tmp1 = tmp$mRNA_Transcription_level
  names(tmp1) = tmp$Name
  tmp = read.table(protein_abundance,header = T)
  tmp2 = tmp$Protein_abundance
  names(tmp2) = toupper(tmp$ORF)
  tmp_id = names(tmp1)
  tmp = tmp2[tmp_id]/tmp1[tmp_id]
  tmp = na.omit(tmp)
  tmp = tmp*(ntotp/ntotr)*sum(tmp1[names(tmp)])/sum(tmp2[names(tmp)]) # 得到比例tmp，然后与R_RNA相乘
  tmp = tmp[genename]
  names(tmp)= genename
  tmp = f_completion(tmp)
  R_PRO = R_RNA
  for (i in 1:ncol(R_RNA)) {
    # 因为R_RNA是按照genename排的序，所以直接使用genename作为标号即可
    R_PRO[,i] = R_RNA[,i]*tmp[genename]
  }
  for (i in 1:ncol(R_PRO)) {
    R_PRO[,i] = f_completion(R_PRO[,i])
  }
  # 统一id顺序----------------------------------------------------
  # 直接按照之前的RNA部分得到的genename排,注意会出现很多NA
  CC_PRO = CC_PRO[genename]
  LPRE_PRO = LPRE_PRO[genename]
  D_PRO = D_PRO[genename]
  R_PRO = R_PRO[genename,]
  # 赋name
  names(CC_PRO) = genename
  names(LPRE_PRO) = genename
  names(D_PRO) = genename
  rownames(R_PRO) = genename
  # 补全
  CC_PRO = f_completion(CC_PRO)
  LPRE_PRO = f_completion(LPRE_PRO)
  D_PRO = f_completion(D_PRO)
  for (i in 1:ncol(R_PRO)) {
    R_PRO[,i] = f_completion(R_PRO[,i])
  }
  # F_PRO，矩阵 ------------------------------------------------------------------
  # F_S_PRO
  F_S_PRO = R_PRO
  for (i in 1:ncol(R_PRO)) {
    F_S_PRO[i] = R_PRO[i]*D_PRO
  }
  result_list <- list(F_S_PRO, CC_PRO, LPRE_PRO)
}



GetRnaProteinInfo <- function(j){
  tissuename= c("Adr","Aorta","BFAT","BS","Cere","Heart",
                "Hypo","Kidney","Liver","Lung","Mus","WFAT")
  tissuenamefull= c("Adrenal Gland","Aorta","Brown fat","Brainstem","Cerebellum","Heart",
                    "Hypothalamus","Kidney","Liver","Lung","Skeletal muscle","White fat")
  RNAtype=c("5utr","3utr","fulllength","coding","intron","protein")
  
  v_yeast = 69
  cmain_yeast = 18.785*10^9
  cgrow_yeast = 2468.2*10^9
  ntotr_yeast = 60*10^3
  ntotp_yeast = 82.2*10^6
  Ne_yeast = 10^7
  
  v_mouse = 2040
  cgrow_mouse = -10^9
  cmain_mouse = -1^9
  ntotr_mouse = 173000
  ntotp_mouse = 3.133*10^9
  Ne_mouse = c(25000,120000)
  rowrange <- function(data){
    unlist(apply(data,1,max)) - unlist(apply(data, 1, min))
  }
  
  # 小鼠 --------------------------------
  ntotr = ntotr_mouse
  ntotp = ntotp_mouse
  Ne = Ne_mouse
  Ne = mean(Ne)
  # 基本符号CC=sumCM，LPRE，D=ln2/HL，R=expr，F ---------------------
  # RNA --------------------------------------------------------------
  # mRNA的序列与长度，sequence_list_DNA，sequence_list_name_DNA，LMAT_RNA
  load("ContainedData/EnergyCost/sequenceDNA.RData")
  sequence_list_name_DNA = sequenceNAME
  sequence_list_name_DNA = unlist(strsplit(sequence_list_name_DNA,'\\|'))[2*(1:length(sequence_list_name_DNA))]
  sequence_list_DNA = sequenceDNA
  names(sequence_list_DNA) = sequence_list_name_DNA
  LMAT_RNA = nchar(sequence_list_DNA)
  names(LMAT_RNA) = sequence_list_name_DNA
  rm(sequenceDNA,sequenceNAME)
  # CC_RNA，向量 -------------------------------------------------
  # 已保存CM_RNA与CC_RNA，可直接读取
  load("ContainedData/EnergyCost/CC_RNA.RData")
  # LPRE_RNA，向量 ----------------------------------------------------
  LPRE_RNA = LMAT_RNA
  # D_RNA，向量 ---------------------------------------------------------
  tmp = read.table("ContainedData/EnergyCost/Mouse_data2.txt",header = T)
  HL_RNA = tmp$mRNA_half.life_h
  D_RNA = log(2)/HL_RNA
  names(D_RNA) = tmp$ID
  D_RNA = na.omit(D_RNA)
  # R_RNA，矩阵 ---------------------------------------------------
  
  tmp = read.csv(paste0("ContainedData/EnergyCost/expression_amplitude/",
                        tissuename[j],'_JTK.csv'))
  tmp_id = tmp$X.1 != "#N/A"
  tmp = tmp[tmp_id,]
  R_RNA = tmp[,8:31]
  rownames(R_RNA) = tmp$X
  # 去掉含有0的行
  R_RNA = R_RNA[rowSums(R_RNA == 0) == 0,]
  # 统一id顺序，genename ----------------------------------------------
  # id对应表
  tmp = read.csv("ContainedData/EnergyCost/ampID.csv")
  # 转换成同一种id，即ensembl转录本id
  # 将R_RNA行名换成ensembl转录本id
  tmp_id = match(rownames(R_RNA),tmp$Probe)
  R_RNA = R_RNA[!is.na(tmp_id),]
  tmp_id = na.omit(tmp_id)
  rownames(R_RNA) = tmp[tmp_id,1]
  # union
  genename = tmp$Ensembl_TransID
  # 排序
  genename = as.character(sort(genename))
  # 其他符号都按这个排
  CC_RNA = CC_RNA[genename]
  D_RNA = D_RNA[genename]
  LPRE_RNA = LPRE_RNA[genename]
  R_RNA = R_RNA[genename,]
  names(CC_RNA) = genename
  names(D_RNA) = genename
  names(LPRE_RNA) = genename
  rownames(R_RNA) = genename
  library(dplyr)
  id_na = CC_RNA %>% is.na()
  CC_RNA[id_na] = median(CC_RNA,na.rm = T)
  id_na = D_RNA %>% is.na()
  D_RNA[id_na] = median(D_RNA,na.rm = T)
  id_na = LPRE_RNA %>% is.na()
  LPRE_RNA[id_na] = median(LPRE_RNA,na.rm = T)
  id_na = R_RNA %>% is.na()
  R_RNA[id_na] = median(unlist(R_RNA),na.rm = T)
  # F_RNA，矩阵 ------------------------------------------------------------------
  # F_S_RNA
  # R*D
  F_S_RNA = R_RNA
  for (i in 1:ncol(R_RNA)) {
    F_S_RNA[i] = R_RNA[i]*D_RNA
  }
  # F_US1_RNA
  # F_US中暂用此方式
  # (R*D+(R2-R1)/(T2-T1))
  # 两端补齐，暂用
  # 两端缺失，在两端补齐结果中去掉两端即可
  tmp = R_RNA[3:24] - R_RNA[1:22]
  tmp = cbind(R_RNA[2]-R_RNA[24],tmp,R_RNA[1]-R_RNA[23])
  tmp = tmp/4
  colnames(tmp) = colnames(R_RNA)
  F_US1_RNA = F_S_RNA + tmp
  # 负数的处理
  # 若一行中负数个数<=2，则直接把负数变为0.1
  tmp_id = F_US1_RNA<0
  tmp_id1 = rowSums(tmp_id)<=2
  tmp_id2 = rowSums(tmp_id) >2
  tmp_id[tmp_id2,] = F
  F_US1_RNA[tmp_id] = 0.1
  # 若一行中负数个数>2，则该行整体加常数至最低值达到0.1
  tmp = apply(F_US1_RNA, 1, min)
  tmp = (-tmp)+0.1
  tmp[tmp_id1] = 0
  for (i in 1:24) {
    F_US1_RNA[i] = F_US1_RNA[i] + tmp
  }
  # # F_US2_RNA
  # # 暂不用此方式
  # # (R2+R1)/2*D+(R2-R1)/(T2-T1)
  # # 其对应的时刻为T1与T2中间，此处用T1表示
  # tmp = (R_RNA[c(2:24,1)] - R_RNA)/2
  # F_US2_RNA = (F_S_RNA[c(2:24,1)] + F_S_RNA)/2 +tmp
  # colnames(F_US2_RNA) = colnames(R_RNA)
  # # 负数的处理
  # # 若一行中负数个数<=2，则直接把负数变为0.1
  # tmp_id = F_US2_RNA<0
  # tmp_id1 = rowSums(tmp_id)<=2
  # tmp_id2 = rowSums(tmp_id) >2
  # tmp_id[tmp_id2,] = F
  # F_US2_RNA[tmp_id] = 0.1
  # # 若一行中负数个数>2，则该行整体加常数至最低值达到0.1
  # tmp = apply(F_US2_RNA, 1, min)
  # tmp = (-tmp)+0.1
  # tmp[tmp_id1] = 0
  # for (i in 1:24) {
  #   F_US2_RNA[i] = F_US2_RNA[i] + tmp
  # }
  
  # PRO --------------------------------------------------------------
  # PRO的序列与长度，sequence_list_name_PRO，sequence_list_PRO，LMAT_PRO
  load(file = 'ContainedData/EnergyCost/sequence_list_PRO_mouse.RData')
  sequence_list_name_PRO = substr(sequence_list_name_PRO,21,nchar(sequence_list_name_PRO))
  # 去掉末尾*号
  sequence_list_PRO = substr(sequence_list_PRO,1,nchar(sequence_list_PRO)-1)
  names(sequence_list_PRO) = sequence_list_name_PRO
  LMAT_PRO = nchar(sequence_list_PRO)
  names(LMAT_PRO) = sequence_list_name_PRO
  # 去掉长度<5的序列
  tmp_id = LMAT_PRO<5
  sequence_list_name_PRO = sequence_list_name_PRO[!tmp_id]
  sequence_list_PRO = sequence_list_PRO[!tmp_id]
  LMAT_PRO = LMAT_PRO[!tmp_id]
  # CC_PRO，向量 ---------------------------------------------------------
  # 已保存CM_PRO与CC_PRO，可直接读取
  load("ContainedData/EnergyCost/CC_PRO.RData")
  # LPRE_PRO，向量 ----------------------------------------------
  LPRE_PRO = LMAT_PRO
  # D_PRO，向量 ----------------------------------------------------
  tmp = read.table("ContainedData/EnergyCost/Mouse_data2.txt",header = T)
  HL_PRO = tmp$Protein_half.life_h
  D_PRO = log(2)/HL_PRO
  names(D_PRO) = tmp$ID
  D_PRO = na.omit(D_PRO)
  # R_PRO，矩阵 -------------------------------------------------------------
  # # 以下计算方式二选一：若所依据实验数据中RNA与PRO相对关系可信则用方式三，否则用方式四
  # # 方式三，通过一次实验中的数据得到R_RNA与R_PRO的比例关系，并用到全体时刻中
  # # mouse中暂用该方式
  # # 若计算F_PRO时使用RIB方式，则不需要计算R_PRO
  # tmp = read.table("E:/RESEARCH_DATA/circadian/cost_model/Mouse_data2.txt",header = T)
  # tmp1 = tmp$mRNA_copy_number
  # names(tmp1) = tmp$ID
  # tmp = read.table("E:/RESEARCH_DATA/circadian/cost_model/Mouse_data2.txt",header = T)
  # tmp2 = tmp$Protein_copy_number
  # names(tmp2) = tmp$ID
  # tmp = tmp2/tmp1
  # tmp = na.omit(tmp) # 得到比例tmp，然后与R_RNA相乘
  # R_PRO = R_RNA
  # for (i in 1:ncol(R_RNA)) {
  #   # 因为R_RNA是按照genename排的序，所以直接使用genename作为标号即可
  #   R_PRO[i] = R_RNA[i]*tmp[genename]
  # }
  # # 方式四，在方式三的基础上再限制sumR_RNA与sumR_PRO的比例
  # # sumR_RNA/sumR_PRO = ntotr/ntotp
  # # 由于R_RNA与R_PRO存在NA和20等数据缺失问题，它们不再代表全体，所以该限制应在tmp1与tmp2阶段做出
  # # 即sumtmp1/sumtmp2 = ntotr/ntotp
  # # 此处两种方法，一种在tmp1与tmp2对齐之前保持比例，一种在对其并去NA之后保持比例。暂用后一种
  # tmp = read.table("E:/RESEARCH_DATA/expression_level/yeast_mRNA_expression.txt",header = T)
  # tmp1 = tmp$mRNA_Transcription_level
  # names(tmp1) = tmp$Name
  # tmp = read.table("E:/RESEARCH_DATA/expression_level/yeast_protein_abundance.txt",header = T)
  # tmp2 = tmp$Protein_abundance
  # names(tmp2) = toupper(tmp$ORF)
  # tmp_id = names(tmp1)
  # tmp = tmp2[tmp_id]/tmp1[tmp_id]
  # tmp = na.omit(tmp)
  # tmp = tmp*(ntotp/ntotr)*sum(tmp1[names(tmp)])/sum(tmp2[names(tmp)]) # 得到比例tmp，然后与R_RNA相乘
  # R_PRO = R_RNA
  # for (i in 1:ncol(R_RNA)) {
  #   # 因为R_RNA是按照genename排的序，所以直接使用genename作为标号即可
  #   R_PRO[i] = R_RNA[i]*tmp[genename]
  # }
  # 统一id顺序----------------------------------------------------
  # 直接按照之前的RNA部分得到的genename排,注意会出现很多NA
  CC_PRO = CC_PRO[genename]
  LPRE_PRO = LPRE_PRO[genename]
  D_PRO = D_PRO[genename]
  # R_PRO = R_PRO[genename,]
  names(CC_PRO) = genename
  names(LPRE_PRO) = genename
  names(D_PRO) = genename
  # rownames(R_PRO) = genename
  library(dplyr)
  id_na = CC_PRO %>% is.na()
  CC_PRO[id_na] = median(CC_PRO,na.rm = T)
  id_na = LPRE_PRO %>% is.na()
  LPRE_PRO[id_na] = median(LPRE_PRO,na.rm = T)
  id_na = D_PRO %>% is.na()
  D_PRO[id_na] = median(D_PRO,na.rm = T)
  #id_na = R_PRO %>% is.na()
  # R_PRO[id_na] = median(unlist(R_PRO),na.rm = T)
  # F_PRO，矩阵 ------------------------------------------------------------------
  
  # # F_S_PRO
  # # R*D
  # F_S_PRO = R_PRO
  # for (i in 1:ncol(R_PRO)) {
  #   F_S_PRO[i] = R_PRO[i]*D_PRO
  # }
  
  # # F_US1_PRO
  # # F_US中暂用此方式
  # # (R*D+(R2-R1)/(T2-T1))
  # # 两端补齐，暂用
  # # 两端缺失，在两端补齐结果中去掉两端即可
  # tmp = R_PRO[3:24] - R_PRO[1:22]
  # tmp = cbind(R_PRO[2]-R_PRO[24],tmp,R_PRO[1]-R_PRO[23])
  # tmp = tmp/4
  # colnames(tmp) = colnames(R_PRO)
  # F_US1_PRO = F_S_PRO + tmp
  # # 负数的处理
  # # 若一行中负数个数<=2，则直接把负数变为0.1
  # tmp_id = F_US1_PRO<0
  # tmp_id1 = rowSums(tmp_id)<=2
  # tmp_id2 = rowSums(tmp_id) >2
  # tmp_id[tmp_id2,] = F
  # F_US1_PRO[tmp_id] = 0.1
  # # 若一行中负数个数>2，则该行整体加常数至最低值达到0.1
  # tmp = apply(F_US1_PRO, 1, min)
  # tmp = (-tmp)+0.1
  # tmp[tmp_id1] = 0
  # for (i in 1:24) {
  #   F_US1_PRO[i] = F_US1_PRO[i] + tmp
  # }
  
  # # F_US2_PRO
  # # 暂不用此方式
  # # (R2+R1)/2*D+(R2-R1)/(T2-T1)
  # # 其对应的时刻为T1与T2中间，此处用T1表示
  # tmp = (R_PRO[c(2:24,1)] - R_PRO)/2
  # F_US2_PRO = (F_S_PRO[c(2:24,1)] + F_S_PRO)/2 +tmp
  # colnames(F_US2_PRO) = colnames(R_PRO)
  # # 负数的处理
  # # 若一行中负数个数<=2，则直接把负数变为0.1
  # tmp_id = F_US2_PRO<0
  # tmp_id1 = rowSums(tmp_id)<=2
  # tmp_id2 = rowSums(tmp_id) >2
  # tmp_id[tmp_id2,] = F
  # F_US2_PRO[tmp_id] = 0.1
  # # 若一行中负数个数>2，则该行整体加常数至最低值达到0.1
  # tmp = apply(F_US2_PRO, 1, min)
  # tmp = (-tmp)+0.1
  # tmp[tmp_id1] = 0
  # for (i in 1:24) {
  #   F_US2_PRO[i] = F_US2_PRO[i] + tmp
  # }
  
  # F_RP_PRO
  # (V_RIB/1000)*(R_RPF/R_RNA)*R_RNA
  # V_RIB
  V_RIB = 39
  V_RIB = V_RIB*3600
  # R_RPF
  # RP(ribosome profiling)与TR(total RNA)的rpkm值
  rpzt = read.csv(
    "ContainedData/EnergyCost/RPZTrpkm.csv",
    row.names = 1)
  trzt = read.csv(
    "ContainedData/EnergyCost/TRZTrpkm.csv",
    row.names = 1)
  trzt = trzt[rownames(rpzt),]
  tmp = rpzt/trzt
  tmp = na.omit(tmp)
  tmp = rowMeans(tmp)
  # 去掉Inf
  tmp_id = tmp == Inf
  tmp = tmp[!tmp_id]
  # 小鼠ensembl转录本id与ensembl基因名id的对应表
  load("ContainedData/EnergyCost/sequenceDNA.RData")
  genename_ENSMUSG = unlist(strsplit(sequenceNAME,'\\|'))[2*(1:length(sequenceNAME))-1]
  genename_ENSMUST = unlist(strsplit(sequenceNAME,'\\|'))[2*(1:length(sequenceNAME))]
  rm(sequenceDNA,sequenceNAME)
  # id转换
  tmp_id1 = match(genename,genename_ENSMUST)
  tmp_id2 = genename_ENSMUSG[tmp_id1]
  tmp = tmp[tmp_id2]
  rpzt = rpzt[tmp_id2,]
  trzt = trzt[tmp_id2,]
  # 此处可以先对rp和tr进行id转换再计算tmp
  names(tmp) = genename
  # tmp填充na
  id_na = tmp %>% is.na()
  tmp[id_na] = median(tmp,na.rm = T)
  # 计算
  F_RP_PRO = (V_RIB/1000)*tmp*R_RNA
  
  
  resultList <- list(F_US1_RNA, CC_RNA, LPRE_RNA, F_RP_PRO, CC_PRO, LPRE_PRO, genename)
  return(resultList)
  
  
}




GetOverlapInTwoPeaks <- function(j, SG, genename){
  
  tissuename= c("Adr","Aorta","BFAT","BS","Cere","Heart",
                "Hypo","Kidney","Liver","Lung","Mus","WFAT")
  
  tmp = read.csv(paste0("ContainedData/EnergyCost/expression_amplitude/",
                        tissuename[j],'_JTK.csv'))
  tmp_id = tmp$X.1 != "#N/A"
  tmp = tmp[tmp_id,]
  R_RNA_withinfo = tmp
  rownames(R_RNA_withinfo) = tmp$X
  # 去掉含有0的行
  R_RNA_withinfo = R_RNA_withinfo[rowSums(R_RNA_withinfo[,8:31] == 0) == 0,]
  # id对应表
  tmp = read.csv("ContainedData/EnergyCost/ampID.csv")
  # 转换成同一种id，即ensembl转录本id
  # 将R_RNA_withinfo行名换成ensembl转录本id
  tmp_id = match(rownames(R_RNA_withinfo),tmp$Probe)
  R_RNA_withinfo = R_RNA_withinfo[!is.na(tmp_id),]
  tmp_id = na.omit(tmp_id)
  rownames(R_RNA_withinfo) = tmp[tmp_id,1]
  # 其他符号都按这个排
  R_RNA_withinfo = R_RNA_withinfo[genename,]
  
  # cir与SG>S0重合关系 
  tmp_id_cir = R_RNA_withinfo$BH.Q<0.05
  tmp_id_cir_mor = R_RNA_withinfo$BH.Q<0.05 & R_RNA_withinfo$LAG<=12
  tmp_id_cir_aft = R_RNA_withinfo$BH.Q<0.05 & R_RNA_withinfo$LAG> 12
  SG_max_wholeday = apply(SG, 1, max)
  SG_min_wholeday = apply(SG, 1, min)
  SG_max_min_wholeday = SG_max_wholeday - SG_min_wholeday
  
  Ne = c(25000,120000)
  Ne = mean(Ne)
  S0 = 4/Ne
  tmp_id_select_wholeday = SG_max_min_wholeday > S0
  # 所有可计算数量
  a = sum(!is.na(tmp_id_cir) & !is.na(tmp_id_select_wholeday))
  # 节律数量
  b = sum(tmp_id_cir & !is.na(tmp_id_select_wholeday))
  b1 = sum(tmp_id_cir_mor & !is.na(tmp_id_select_wholeday))
  b2 = sum(tmp_id_cir_aft & !is.na(tmp_id_select_wholeday))
  # 选择数量
  c = sum(tmp_id_select_wholeday,na.rm = T)
  # 节律与选择重合数量
  d = sum(tmp_id_cir & tmp_id_select_wholeday,na.rm = T)
  d1 = sum(tmp_id_cir_mor & tmp_id_select_wholeday,na.rm = T)
  d2 = sum(tmp_id_cir_aft & tmp_id_select_wholeday,na.rm = T)
  
  result_list <- list(a,b,b1,b2,c,d,d1,d2)
  return(result_list)
}




















library(dplyr)
library(tidyr)

p_fun = function(mean,sd2){mean/sd2}
r_fun = function(mean,sd2){mean^2/(sd2-mean)}
point_pcc = function(x,y,xlab='',ylab='',main='',textx=0,text1y=0,text2y=0,cex=1){
  plot(x,y,
       xlab = xlab, ylab = ylab, main = main,
       pch = 20)
  text_pvalue = cor.test(x,y)$p.value %>% 
    round(digits = 3)
  text_cor = cor.test(x,y)$estimate %>% 
    round(digits = 3)
  text(textx,text1y,paste0('pcc p = ',text_pvalue),
       pos = 4, cex = cex,offset = 0)
  text(textx,text2y,paste0('pcc r = ',text_cor),
       pos = 4,cex = cex,offset = 0)
}
n_level_num = function(data,id,n=5){
  data_num = length(data)
  id_ordered = id[order(data)]
  id_group = round((data_num/n)*1:n)
  id_group_1 = c(0,id_group[-n])+1
  id_group_2 = id_group
  inlevel_num = 1:n
  for (i in 1:n) {
    inlevel_num[i] = sum(id_ordered[id_group_1[i]:id_group_2[i]])
  }
  return(inlevel_num)
}
f_goterm_genenumber_dist = function(goclass='all',level_go='max number',funcnumpergene = F){
  edge_go_dags = read.table("ContainedData/Simulation/edge_go_dags.txt",
                            stringsAsFactors = F)
  vertex_go_dags = read.table("ContainedData/Simulation/vertex_go_dags.txt",
                              stringsAsFactors = F)
  {
    mus_gomap = read.table("ContainedData/Simulation/mus_gomap.txt",
                           stringsAsFactors = F)
    # mgi2ens读入与处理
    # 包org.Hs.eg.db内容不全更新不足，不使用; 使用从MGI网站下载的id对应数据
    mgi2ens = read.table("ContainedData/Simulation/MGI_to_ENSEMBL.txt",
                         header = T,sep = "\t",quote = "",stringsAsFactors = F)
    mgi2ens = mgi2ens[-2]
    names(mgi2ens)[c(1,2)] = c("ENSMUST_ID","MGI_ID")
    mgi2ens = mgi2ens[mgi2ens$MGI_ID!="No associated gene",]
    mgi2ens = mgi2ens[!duplicated(mgi2ens$MGI_ID),]
  }
  mgi_id_set = mus_gomap[,1] %>% unique()
  {
    termseed = c("GO:0008150","GO:0003674","GO:0005575")
    termi = termseed
    termrank = matrix(0,nrow = nrow(vertex_go_dags),ncol = 20)
    rownames(termrank) = vertex_go_dags[,1]
    k=1
    while(length(termi)!=0){
      termrank[termi,k] = k
      indbool2 = edge_go_dags[,2] %in% termi
      termi = edge_go_dags[indbool2,1]
      k=k+1
    }
    termrank_l = apply(termrank,1,max)
    print(table(termrank_l))
    termrank_s = apply(termrank,1,function(x){min(x[x>0])})
    print(table(termrank_s))
  }
  
  
  
  
  
  
  if (goclass =='all'){ 
    id =T}  else if(goclass=='BP'){ 
      id = vertex_go_dags[,2] =='P'}  else if(goclass=='CC'){ 
        id = vertex_go_dags[,2] =='C'}  else if(goclass=='MF'){
          id = vertex_go_dags[,2] =='F'}
  termrank_l_sub = termrank_l[id]
  table(termrank_l_sub) %>% print()
  if (level_go == 'max number') {level_go = table(termrank_l_sub) %>% which.max()}
  ind_tmp = termrank_l_sub==level_go
  termrank_l_sub_level = termrank_l_sub[ind_tmp] %>% names()
  gene_num_per_func = c()
  tmp_genelist = c()
  for (i in 1:length(termrank_l_sub_level)) {
    #print(i)
    {
      termseed = termrank_l_sub_level[i]
      termi = termseed
      indbool1 = F
      while(length(termi)!=0){
        indbool2 = edge_go_dags[,2] %in% termi
        indbool1 = indbool1 | indbool2
        termi = edge_go_dags[indbool2,1]
      }
      termall = union(edge_go_dags[indbool1,1],termseed)
      mgi_id_set_termseed = unique(mus_gomap[mus_gomap[,3] %in% termall, 1])
    }
    gene_num_per_func = c(gene_num_per_func,length(mgi_id_set_termseed))
    tmp_genelist = c(tmp_genelist,mgi_id_set_termseed)
  }
  tmp_hash = cbind(rep(termrank_l_sub_level,gene_num_per_func),tmp_genelist)
  table(termrank_l_sub) %>% print()
  gene_num_per_func = gene_num_per_func[gene_num_per_func>0]
  func_num_per_gene = tmp_hash[,2] %>% table()
  paste0('GO level: ',level_go,'\n') %>% cat()
  paste0('term number: ',length(gene_num_per_func),'\n') %>% cat()
  paste0('gene number: ',length(func_num_per_gene),'\n') %>% cat()
  # 一个term对应多少gene，的分布
  # 是个什么形状
  plot(table(gene_num_per_func),
       xlab = 'gene number',ylab = 'function number',main = paste0('GO-',goclass,' level ',level_go),
       log = 'x')
  # 双对数坐标
  tmp_table = log(gene_num_per_func) %>% table() %>% log()
  tmp_table = tmp_table[tmp_table>log(1)]
  tmp_table_name = tmp_table %>% names() %>% as.numeric()
  plot(tmp_table,
       xlab = 'log(gene number)',ylab = 'log(function number)',main = paste0('GO-',goclass,' level ',level_go),
       xaxt = 'n')
  axis(side = 1,at = tmp_table_name,labels = round(tmp_table_name,digits = 2))
  # 根据拟合的直线求拟合曲线
  plot(tmp_table_name,tmp_table,
       xlab = 'log(gene number)',ylab = 'log(function number)',main = paste0('GO-',goclass,' level ',level_go),
       pch = 20,
       yaxt = 'n')
  axis(side = 2,at = 0:10,labels = 0:10)
  # 拟合直线
  tmp_mean = aggregate(tmp_table_name,tmp_table %>% as.vector() %>% list(),mean)
  plot(tmp_mean[,2],tmp_mean[,1],
       xlab = 'log(gene number)',ylab = 'log(function number)',main = paste0('GO-',goclass,' level ',level_go),
       pch = 20)
  tmp_lm = lm(tmp_mean[,1] ~ tmp_mean[,2])
  tmp_lm
  abline(tmp_lm)
  # 拟合曲线
  plot(table(gene_num_per_func),
       xlab = 'gene number',ylab = 'function number',main = paste0('GO-',goclass,' level ',level_go),
       log = 'x')
  alpha = tmp_lm$coefficients[2]
  beta = exp(tmp_lm$coefficients[1])
  lines(1:3000,(1:3000)^alpha*beta)
  paste0('coefficients: ',alpha %>% round(digits = 3),', ',beta %>% round(digits = 3),'\n') %>% cat()
  if (funcnumpergene==T) {
    # 一个gene对应多少term，的分布
    # 双对数
    tmp_table = log(func_num_per_gene) %>% table() %>% log()
    tmp_table = tmp_table[tmp_table>log(1)]
    tmp_table_name = tmp_table %>% names() %>% as.numeric()
    plot(tmp_table,
         xlab = 'log(function number)',ylab = 'log(gene number)',main = paste0('GO-',goclass,' level ',level_go),
         xaxt = 'n')
    axis(side = 1,at = tmp_table_name,labels = round(tmp_table_name,digits = 2))
    # 拟合直线
    tmp_mean = aggregate(tmp_table_name,tmp_table %>% as.vector() %>% list(),mean)
    plot(tmp_mean[,2],tmp_mean[,1],
         xlab = 'log(function number)',ylab = 'log(gene number)',main = paste0('GO-',goclass,' level ',level_go),
         pch = 20)
    tmp_lm = lm(tmp_mean[,1] ~ tmp_mean[,2])
    tmp_lm
    abline(tmp_lm)
    # 拟合曲线
    plot(table(func_num_per_gene),
         xlab = 'function number',ylab = 'gene number',main = paste0('GO-',goclass,' level ',level_go),
         log = '')
    alpha = tmp_lm$coefficients[2]
    beta = exp(tmp_lm$coefficients[1])
    lines(1:3000,(1:3000)^alpha*beta)
    paste0('coefficients: ',alpha %>% round(digits = 3),', ',beta %>% round(digits = 3),'') %>% cat()
    # 用负二项拟合试一下
    paste0('\n','length(mgi_id_set):',length(mgi_id_set),
           '   0 number:',length(mgi_id_set)-length(func_num_per_gene)) %>% cat()
    func_num_per_gene_0 = c(func_num_per_gene,rep(0,length(mgi_id_set)-length(func_num_per_gene)))
    tmpr = r_fun(mean = mean(func_num_per_gene_0),sd2 = sd(func_num_per_gene_0)^2)
    tmpp = p_fun(mean = mean(func_num_per_gene_0),sd2 = sd(func_num_per_gene_0)^2)
    tmp_dnb = dnbinom(0:1000,tmpr,tmpp)
    # 把x和y取log
    plot(tmp_mean[,2],tmp_mean[,1],
         xlab = 'log(function number)',ylab = 'log(gene number)',main = paste0('GO-',goclass,' level ',level_go),
         pch = 20)
    lines(1:1000 %>% log(),(table(func_num_per_gene_0)[1]/tmp_dnb[1]*tmp_dnb[-1]) %>% log(),
          col='red')
    # 只把y取log
    plot(func_num_per_gene_0 %>% table() %>% names() %>% as.numeric(),
         table(func_num_per_gene_0) %>% as.vector(),
         xlab = 'function number',ylab = 'gene number',main = paste0('GO-',goclass,' level ',level_go),
         log = 'y',
         pch = 20)
    lines(0:1000,table(func_num_per_gene_0)[1]/tmp_dnb[1]*tmp_dnb,
          col='red')
    # y值相同的取均值
    tmp_table_ylog = func_num_per_gene_0 %>% table()
    tmp_table_name_ylog = tmp_table_ylog %>% names() %>% as.numeric()
    tmp_mean_ylog = aggregate(tmp_table_name_ylog,tmp_table_ylog %>% as.vector() %>% list(),mean)
    plot(tmp_mean_ylog[,2],tmp_mean_ylog[,1],
         xlab = 'function number',ylab = 'gene number',main = paste0('GO-',goclass,' level ',level_go),
         log = 'y',
         pch = 20)
    lines(0:1000,table(func_num_per_gene_0)[1]/tmp_dnb[1]*tmp_dnb,
          col='red')
    # 原柱状图并添加线
    plot(table(func_num_per_gene_0),
         xlab = 'function number',ylab = 'gene number',main = paste0('GO-',goclass,' level ',level_go),
         log = '')
    lines(0:1000,table(func_num_per_gene_0)[1]/tmp_dnb[1]*tmp_dnb,
          col='red')
    # 打印参数
    paste0('\n',
           'coefficients r:',tmpr %>% round(digits = 3),
           ', p:',tmpp %>% round(digits = 3),
           ', mean:',mean(func_num_per_gene_0) %>% round(digits = 3),
           ', sd2:',(sd(func_num_per_gene_0)^2) %>% round(digits = 3)) %>% 
      cat()
  }
}
f_gene_goterm_matrix = function(goclass='all',level_go='max number'){
  
  edge_go_dags = read.table("ContainedData/Simulation/edge_go_dags.txt",
                            stringsAsFactors = F)
  vertex_go_dags = read.table("ContainedData/Simulation/vertex_go_dags.txt",
                              stringsAsFactors = F)
  {
    mus_gomap = read.table("ContainedData/Simulation/mus_gomap.txt",
                           stringsAsFactors = F)
    # mgi2ens读入与处理
    # 包org.Hs.eg.db内容不全更新不足，不使用; 使用从MGI网站下载的id对应数据
    mgi2ens = read.table("ContainedData/Simulation/MGI_to_ENSEMBL.txt",
                         header = T,sep = "\t",quote = "",stringsAsFactors = F)
    mgi2ens = mgi2ens[-2]
    names(mgi2ens)[c(1,2)] = c("ENSMUST_ID","MGI_ID")
    mgi2ens = mgi2ens[mgi2ens$MGI_ID!="No associated gene",]
    mgi2ens = mgi2ens[!duplicated(mgi2ens$MGI_ID),]
  }
  mgi_id_set = mus_gomap[,1] %>% unique()
  {
    termseed = c("GO:0008150","GO:0003674","GO:0005575")
    termi = termseed
    termrank = matrix(0,nrow = nrow(vertex_go_dags),ncol = 20)
    rownames(termrank) = vertex_go_dags[,1]
    k=1
    while(length(termi)!=0){
      termrank[termi,k] = k
      indbool2 = edge_go_dags[,2] %in% termi
      termi = edge_go_dags[indbool2,1]
      k=k+1
    }
    termrank_l = apply(termrank,1,max)
    print(table(termrank_l))
    termrank_s = apply(termrank,1,function(x){min(x[x>0])})
    print(table(termrank_s))
  }
  
  
  if (goclass =='all'){ 
    id =T}  else if(goclass=='BP'){ 
      id = vertex_go_dags[,2] =='P'}  else if(goclass=='CC'){ 
        id = vertex_go_dags[,2] =='C'}  else if(goclass=='MF'){
          id = vertex_go_dags[,2] =='F'}
  termrank_l_sub = termrank_l[id]
  table(termrank_l_sub) %>% print()
  if (level_go == 'max number') {level_go = table(termrank_l_sub) %>% which.max()}
  ind_tmp = termrank_l_sub==level_go
  termrank_l_sub_level = termrank_l_sub[ind_tmp] %>% names()
  gene_num_per_func = c()
  tmp_genelist = c()
  gene_term_list = c()
  for (i in 1:length(termrank_l_sub_level)) {
    #print(i)
    {
      termseed = termrank_l_sub_level[i]
      termi = termseed
      indbool1 = F
      while(length(termi)!=0){
        indbool2 = edge_go_dags[,2] %in% termi
        indbool1 = indbool1 | indbool2
        termi = edge_go_dags[indbool2,1]
      }
      termall = union(edge_go_dags[indbool1,1],termseed)
      mgi_id_set_termseed = unique(mus_gomap[mus_gomap[,3] %in% termall, 1])
      if (length(mgi_id_set_termseed)>0) {
        gene_term_list = rbind(gene_term_list, cbind(mgi_id_set_termseed,termseed))
      }
    }
    gene_num_per_func = c(gene_num_per_func,length(mgi_id_set_termseed))
    tmp_genelist = c(tmp_genelist,mgi_id_set_termseed)
  }
  tmp_hash = cbind(rep(termrank_l_sub_level,gene_num_per_func),tmp_genelist)
  table(termrank_l_sub) %>% print()
  gene_num_per_func = gene_num_per_func[gene_num_per_func>0]
  func_num_per_gene = tmp_hash[,2] %>% table()
  paste0('GO level: ',level_go,'\n') %>% cat()
  paste0('term number: ',length(gene_num_per_func),'\n') %>% cat()
  paste0('gene number: ',length(func_num_per_gene),'\n') %>% cat()
  gene_term_matrix = table(gene_term_list[,1], gene_term_list[,2]) %>% as.matrix()
  paste0('term number: ',ncol(gene_term_matrix),'\n') %>% cat()
  paste0('gene number: ',nrow(gene_term_matrix),'\n') %>% cat()
  return(gene_term_matrix)
}