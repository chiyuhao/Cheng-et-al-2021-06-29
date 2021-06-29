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
rowrange = function(data){
  unlist(apply(data,1,max)) - unlist(apply(data, 1, min))
}
MouseCalculateCost <- function(){
  {
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
  i=9
  tmp = read.csv(paste0("ContainedData/EnergyCost/expression_amplitude/",
                        tissuename[i],'_JTK.csv'))
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
  genename = as.character(tmp$Ensembl_TransID)
  # 排序
  genename = sort(genename)
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
  # id_na = R_PRO %>% is.na()
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
  # 单个基因的能量成本CG，矩阵----------------------------------------
  # 公式一 ---------------------------------------------------------
  # 单体从头合成
  {
    # F_RNA = F_US1_RNA
    # F_PRO = F_RP_PRO
    # CG_RNA = F_RNA
    # for (i in 1:ncol(F_RNA)) {
    #   CG_RNA[i] = F_RNA[i]*CC_RNA
    # }
    # CG_PRO = F_PRO
    # for (i in 1:ncol(F_PRO)) {
    #   CG_PRO[i] = F_PRO[i]*CC_PRO
    # }
    # CG = CG_RNA+CG_PRO
  }
  # 公式二 ----------------------------------------------------
  # R*CC+2.17*F*T*LPRE
  # R*CC+5*F*L*T-R*L
  {
    # F_RNA = F_S_RNA
    # F_PRO = F_S_PRO
    # T_0 = 1
    # CG_RNA = F_RNA
    # for (i in 1:ncol(F_RNA)) {
    #   T_i = T0+2*i-2
    #   CG_RNA[i] = F_RNA[i]*CC_RNA + 2.17*F_RNA[i]*LPRE_RNA*T_i
    # }
    # CG_PRO = F_PRO
    # for (i in 1:ncol(F_PRO)) {
    #   T_i = T0+2*i-2
    #   CG_PRO[i] = F_PRO[i]*CC_PRO + 5*F_PRO[i]*LPRE_PRO*T_i - R_PRO[i]*LPRE_PRO
    # }
    # CG = CG_RNA+CG_PRO
  }
  # 公式三 ----------------------------------------------------------
  # F*(2.17*LPRE+sigma*CC)
  # F*(5*L+sigma*CC)
  {
    F_RNA = F_US1_RNA
    F_PRO = F_RP_PRO
    sigma_RNA = 0.5
    sigma_PRO = 0.5
    CG_RNA = F_RNA
    for (i in 1:ncol(F_RNA)) {
      CG_RNA[i] = F_RNA[i]*(2.17*LPRE_RNA + sigma_RNA*CC_RNA)
    }
    CG_PRO = F_PRO
    for (i in 1:ncol(F_PRO)) {
      CG_PRO[i] = F_PRO[i]*(5*LPRE_PRO + sigma_PRO*CC_PRO)
    }
    CG = CG_RNA+CG_PRO
  }
  # 细胞总成本CT，向量 --------------------------------------------------
  CT = 20*colSums(na.omit(CG_RNA))+colSums(na.omit(CG_PRO))
  # SG，矩阵 --------------------------------------------------------
  SG = CG
  for (i in 1:ncol(CG)) {
    SG[i] = CG[i]/CT[i]
  }
  return(SG)
}

















YeastCalculateCost <- function(){
  {
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
  }
  f_completion = function(data,method='median'){
    id_na = is.na(data)
    if(method=='median'){data[id_na] = median(data,na.rm = T)}
    return(data)
  }
  f_NamesortCompletion = function(dataVector,nameVector,method='median'){
    dataVector = dataVector[nameVector]
    names(dataVector) = nameVector
    dataVector = f_completion(dataVector)
    return(dataVector)
  }
  rowrange = function(data){
    unlist(apply(data,1,max)) - unlist(apply(data, 1, min))
  }
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
  
  # 酵母 --------------------------------
  ntotr = ntotr_yeast
  ntotp = ntotp_yeast
  Ne = Ne_yeast
  tmp = read.csv(
    'ContainedData/EnergyCost/meta2d_sample2.csv')
  genename = as.character(tmp$CycID %>% sort())
  # 基本符号CC=sumCM，LPRE，D=ln2/HL，R=expr，F ---------------------
  # RNA --------------------------------------------------------------
  # mRNA的序列与长度，sequence_list_DNA，sequence_list_name_DNA，LMAT_RNA
  con = file("ContainedData/EnergyCost/yeast_codingseq_clean_All2.txt", "r")
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
  tmp = read.table("ContainedData/EnergyCost/mart_export_transcript_length.txt",
                   header = T,sep = '\t')
  LPRE_RNA = tmp$Transcript.end..bp.-tmp$Transcript.start..bp.+1
  names(LPRE_RNA) = tmp$Gene.stable.ID
  # D_RNA，向量 ---------------------------------------------------------
  tmp = read.table("ContainedData/EnergyCost/yeast_parameters2.txt",header = T)
  HL_RNA = tmp$mRNA_halflive
  D_RNA = log(2)/HL_RNA
  names(D_RNA) = tmp$ID
  # R_RNA，矩阵 ---------------------------------------------------
  tmp = read.csv(
    "ContainedData/EnergyCost/meta2d_sample2.csv")
  R_RNA = tmp[,24:43]
  rownames(R_RNA) = tmp$CycID
  # 去掉含有0的行
  # R_RNA = R_RNA[rowSums(R_RNA == 0) == 0,]
  # 统一id顺序，genename ----------------------------------------------
  # union
  # genename = Reduce(union,list(v1 = names(CC_RNA),v2 = names(D_RNA),v3 = names(LPRE_RNA),v4 = rownames(R_RNA)))
  # 排序
  # genename = sort(genename)
  # 其他符号都按这个排
  CC_RNA = f_NamesortCompletion(CC_RNA,genename) 
  D_RNA = f_NamesortCompletion(D_RNA,genename) 
  LPRE_RNA = f_NamesortCompletion(LPRE_RNA,genename) 
  R_RNA = R_RNA[genename,]
  
  # F_RNA，矩阵 ------------------------------------------------------------------
  # F_S_RNA
  F_S_RNA = R_RNA
  for (i in 1:ncol(R_RNA)) {
    F_S_RNA[i] = R_RNA[i]*D_RNA
  }
  
  tmp = R_RNA[3:20] - R_RNA[1:18]
  tmp = cbind(R_RNA[2]-R_RNA[20],tmp,R_RNA[1]-R_RNA[19])
  tmp = tmp/(26/60)
  colnames(tmp) = colnames(R_RNA)
  F_US1_RNA = F_S_RNA + tmp
  
  tmp_id = F_US1_RNA<0
  tmp_id1 = rowSums(tmp_id)<=2
  tmp_id2 = rowSums(tmp_id) >2
  tmp_id[tmp_id2,] = F
  F_US1_RNA[tmp_id] = 0.1
  # 若一行中负数个数>2，则该行整体加常数至最低值达到0.1
  tmp = apply(F_US1_RNA, 1, min)
  tmp = (-tmp)+0.1
  tmp[tmp_id1] = 0
  for (i in 1:20) {
    F_US1_RNA[i] = F_US1_RNA[i] + tmp
  }
  
  
  # F_US1_RNA
  # F_US2_RNA
  # PRO --------------------------------------------------------------
  # PRO的序列与长度，sequence_list_name_PRO，sequence_list_PRO，LMAT_PRO
  con = file("ContainedData/EnergyCost/yeast_protein_clean_All2.txt", "r")
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
  tmp = read.table("ContainedData/EnergyCost/yeast_parameters2.txt",header = T)
  HL_PRO = tmp$Protein_halflive
  D_PRO = log(2)/HL_PRO
  names(D_PRO) = tmp$ID
  # R_PRO，矩阵 -------------------------------------------------------------
  # 以下计算方式二选一：若所依据实验数据中RNA与PRO相对关系可信则用方式三，否则用方式四
  # # 方式三，通过一次实验中的数据得到R_RNA与R_PRO的比例关系，并用到全体时刻中
  # tmp = read.table("E:/RESEARCH_DATA/expression_level/yeast_mRNA_expression.txt",header = T)
  # tmp1 = tmp$mRNA_Transcription_level
  # names(tmp1) = tmp$Name
  # tmp = read.table("E:/RESEARCH_DATA/expression_level/yeast_protein_abundance.txt",header = T)
  # tmp2 = tmp$Protein_abundance
  # names(tmp2) = toupper(tmp$ORF)
  # tmp_id = names(tmp1)
  # tmp = tmp2[tmp_id]/tmp1[tmp_id]
  # tmp = na.omit(tmp) # 得到比例tmp，然后与R_RNA相乘
  # R_PRO = R_RNA
  # for (i in 1:ncol(R_RNA)) {
  #   # 因为R_RNA是按照genename排的序，所以直接使用genename作为标号即可
  #   R_PRO[i] = R_RNA[i]*tmp[genename]
  # }
  # for (i in 1:ncol(R_PRO)) {
  #   R_PRO[,i] = f_completion(R_PRO[,i])
  # }
  
  # 方式四，在方式三的基础上再限制sumR_RNA与sumR_PRO的比例
  # yeast中暂用该方式
  # sumR_RNA/sumR_PRO = ntotr/ntotp
  # 由于R_RNA与R_PRO存在NA和20等数据缺失问题，它们不再代表全体，所以该限制应在tmp1与tmp2阶段做出
  # 即sumtmp1/sumtmp2 = ntotr/ntotp
  # 此处两种方法，一种在tmp1与tmp2对齐之前保持比例，一种在对其并去NA之后保持比例。暂采用后一种
  tmp = read.table("ContainedData/EnergyCost/yeast_mRNA_expression.txt",header = T)
  tmp1 = tmp$mRNA_Transcription_level
  names(tmp1) = tmp$Name
  tmp = read.table("ContainedData/EnergyCost/yeast_protein_abundance.txt",header = T)
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
  CC_PRO = f_NamesortCompletion(CC_PRO,genename)
  LPRE_PRO = f_NamesortCompletion(LPRE_PRO,genename)
  D_PRO = f_NamesortCompletion(D_PRO,genename)
  R_PRO = R_PRO[genename,]
  # F_PRO，矩阵 ------------------------------------------------------------------
  # F_S_PRO
  F_S_PRO = R_PRO
  for (i in 1:ncol(R_PRO)) {
    F_S_PRO[i] = R_PRO[i]*D_PRO
  }
  # F_US1_PRO
  # F_US2_PRO
  # F_RP_PRO
  
  # tmp = read.table("E:/RESEARCH_DATA/circadian/cost_model/yeast_ribo_profiling/GSE56622_ZidProcessedDataAll.txt",
  #                  sep = '\t',header = T,quote = "",fill = T)
  # tmp = tmp[,!(names(tmp) %in% c('X','X.1','X.2','X.3'))]
  # tmp = na.omit(tmp)
  #   # The ribosome occupancy along the mRNA was calculated by dividing the ribosome read counts 
  #   # at each base pair along the gene by the average mRNA reads per base pair for each gene.
  # mrnatopro = (unlist(tmp['Ribo.Glu_RPKM.BY4741...Replicate.1']))/
  #   (unlist(tmp['RNA.Glu_RPKM.BY4741.TotalRNA...Replicate.1']))
  # names(mrnatopro) = tmp[,1]
  # sssssssssss
  
  # 单个基因的能量成本CG，矩阵----------------------------------------
  # 公式一 ---------------------------------------------------------
  # # 单体从头合成
  # CG_RNA = F_S_RNA
  # for (i in 1:ncol(F_S_RNA)) {
  #   CG_RNA[i] = F_S_RNA[i]*CC_RNA
  # }
  # CG_PRO = F_S_PRO
  # for (i in 1:ncol(F_S_PRO)) {
  #   CG_PRO[i] = F_S_PRO[i]*CC_PRO
  # }
  # CG = CG_RNA+CG_PRO
  # 公式二 ----------------------------------------------------
  
  # 公式三 ----------------------------------------------------------
  # F*(2.17*LPRE+sigma*CC)
  # F*(5*L+sigma*CC)
  {
    F_RNA = F_US1_RNA
    F_PRO = F_S_PRO
    sigma_RNA = 0.5
    sigma_PRO = 0.5
    CG_RNA = F_RNA
    for (i in 1:ncol(F_RNA)) {
      CG_RNA[i] = F_RNA[i]*(2.17*LPRE_RNA + sigma_RNA*CC_RNA)
    }
    CG_PRO = F_PRO
    for (i in 1:ncol(F_PRO)) {
      CG_PRO[i] = F_PRO[i]*(5*LPRE_PRO + sigma_PRO*CC_PRO)
    }
    CG = CG_RNA+CG_PRO
  }
  # 细胞总成本CT，向量 --------------------------------------------------
  CT = 20*colSums(na.omit(CG_RNA))+colSums(na.omit(CG_PRO))
  # SG，矩阵 --------------------------------------------------------
  SG = CG
  for (i in 1:ncol(CG)) {
    SG[i] = CG[i]/CT[i]
  }
  result <- list()
  result[[1]] <- SG
  result[[2]] <- genename
  result[[3]] <- R_RNA
  return(result)
}