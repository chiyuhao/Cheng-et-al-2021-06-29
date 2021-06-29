library(pheatmap)
library(VennDiagram)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(dplyr)
library(forcats)
library(pheatmap)
library(VennDiagram)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(dplyr)
library(forcats)

get_or_gene <- function(){
species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','chlamydomonas','cyanobacteria')
species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','chlamydomonas','cyanobacteria')
gene_num_list <- data.frame(0,0)
colnames(gene_num_list) <- c("species","gene_num")
venn_data <- data.frame(0,0,0,0,0,0,0,0,0,0)
colnames(venn_data) <- c("species1","species2","num1","num2","overlap","num1_or","num2_or","overlap_or","background","or")
for(i in 1:(length(species) - 1))
{
  for(j in (i+1):length(species))
  {
    #Get orthologue gene in current 2 species
    current_file_name <- paste(paste(paste(paste("ContainedData/Plot_required_file/protein2gene_20210122/",species[i],sep=""),'__v__',sep=""), species[j],sep=""),".txt",sep="")
    #Determine whether the file name exists
    if (file.exists(current_file_name))
    {
      gene_list <- unique(read.table(current_file_name))
      #Get cycling gene in curent 2 species
      current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
      
      current_background_gene1_name <- paste0(paste0("ContainedData/Plot_required_file/background/", species[i]), "_background.txt")
      current_background_gene2_name <- paste0(paste0("ContainedData/Plot_required_file/background/", species[j]), "_background.txt")
      
      current_background_gene1 <-read.table(current_background_gene1_name)$V1
      current_background_gene2 <-read.table(current_background_gene2_name)$V1
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      gene1_num <- 0
      gene2_num <- 0
      overlap <- 0
      gene_num1_new <- 0
      gene_num2_new <- 0
      overlap_new <- 0
      background <- 0
      num <- 0
      gene1_total_num <- c()
      gene2_total_num <- 0
      or <- 0
      for(k in 1:ncol(gene1))
      {
        for(l in 1:ncol(gene2))
        {
          num <- num + 1
          if(which(species_sort == species[i]) < which(species_sort == species[j]))
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
            current_gene1_background <- gene_list[gene_list$V1 %in% current_background_gene1,]$V1
            current_gene2_background <- gene_list[gene_list$V2 %in% current_background_gene2,]$V1
            intersect_gene <- intersect(current_gene2_background, current_gene1_background)
            current_gene1_new <- unique(current_gene1[current_gene1 %in% intersect_gene])
            current_gene2_new <- unique(current_gene2[current_gene2 %in% intersect_gene])
          }else
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
            current_gene1_background <- gene_list[gene_list$V1 %in% current_background_gene1,]$V2
            current_gene2_background <- gene_list[gene_list$V2 %in% current_background_gene2,]$V2
            intersect_gene <- intersect(current_gene2_background, current_gene1_background)
            current_gene1_new <- unique(current_gene1[current_gene1 %in% intersect_gene])
            current_gene2_new <- unique(current_gene2[current_gene2 %in% intersect_gene])
            
            
             
            
            
            
          }
          
          gene_num1_new <- gene_num1_new + length(current_gene1_new)
          gene_num2_new <- gene_num2_new + length(current_gene2_new)
          overlap_new <- overlap_new + length(intersect(current_gene1_new, current_gene2_new))
          background <- background + length(intersect_gene)
          overlap <- overlap + length(intersect(current_gene1, current_gene2))
          id_cir1 <-  intersect_gene %in% current_gene1_new
          id_cir2 <-  intersect_gene %in% current_gene2_new
          
          res_fishertest = fisher.test(id_cir1,id_cir2)
          or <- or + res_fishertest$estimate
        }
      }
      gene1_num <- gene1_num / num
      gene2_num <- gene2_num / num
      overlap <- overlap / num
      or <- or / num
      gene_num1_new <- gene_num1_new / num
      gene_num2_new <- gene_num2_new /num
      overlap_new <- overlap_new /num
      background <- background / num
      temp_row <- c(species[i], species[j], gene1_num, gene2_num, overlap, gene_num1_new, gene_num2_new, overlap_new, background, or)
      venn_data <- rbind(venn_data, temp_row)
      gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene1_total_num))))
      gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene2_total_num))))
    }else
    {
      current_file_name <- paste(paste(paste(paste("ContainedData/Plot_required_file/protein2gene_20210122/",species[j],sep=""),'__v__',sep=""), species[i],sep=""),".txt",sep="")
      gene_list <- unique(read.table(current_file_name))
      current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
      current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
      gene1 <- read.csv(current_file_gene1_name,header = T)
      gene2 <- read.csv(current_file_gene2_name,header = T)
      
      
      current_background_gene1_name <- paste0(paste0("ContainedData/Plot_required_file/background/", species[j]), "_background.txt")
      current_background_gene2_name <- paste0(paste0("ContainedData/Plot_required_file/background/", species[i]), "_background.txt")
      current_background_gene1 <-read.table(current_background_gene1_name)$V1
      current_background_gene2 <-read.table(current_background_gene2_name)$V1
      
      gene1_num <- 0
      gene2_num <- 0
      overlap <- 0
      or <- 0
      num <- 0
      gene_num1_new <- 0
      gene_num2_new <- 0
      overlap_new <- 0
      background <- 0
      for(k in 1:ncol(gene1))
      {
        for(l in 1:ncol(gene2))
        {
          num <- num + 1
          if(which(species_sort == species[i]) > which(species_sort == species[j]))
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
            
            current_gene1_background <- gene_list[gene_list$V1 %in% current_background_gene1,]$V1
            current_gene2_background <- gene_list[gene_list$V2 %in% current_background_gene2,]$V1
            intersect_gene <- intersect(current_gene2_background, current_gene1_background)
            current_gene1_new <- unique(current_gene1[current_gene1 %in% intersect_gene])
            current_gene2_new <- unique(current_gene2[current_gene2 %in% intersect_gene])
            
          }else
          {
            current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
            current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
            gene1_num <- gene1_num + length(unique(current_gene1))
            gene2_num <- gene2_num + length(unique(current_gene2))
            
            
            current_gene1_background <- gene_list[gene_list$V1 %in% current_background_gene1,]$V2
            current_gene2_background <- gene_list[gene_list$V2 %in% current_background_gene2,]$V2
            intersect_gene <- intersect(current_gene2_background, current_gene1_background)
            current_gene1_new <- unique(current_gene1[current_gene1 %in% intersect_gene])
            current_gene2_new <- unique(current_gene2[current_gene2 %in% intersect_gene])
            

            
          }
          gene_num1_new <- gene_num1_new + length(current_gene1_new)
          gene_num2_new <- gene_num2_new + length(current_gene2_new)
          overlap_new <- overlap_new + length(intersect(current_gene1_new, current_gene2_new))
          background <- background + length(intersect_gene)
          
          
          
          overlap <- overlap + length(intersect(current_gene1, current_gene2))
          
          id_cir1 <-  intersect_gene %in% current_gene1_new
          id_cir2 <-  intersect_gene %in% current_gene2_new
          
          res_fishertest = fisher.test(id_cir1,id_cir2)
          or <- or + res_fishertest$estimate
          print(res_fishertest$p.value)
          print(res_fishertest$estimate)
        }
      }
      gene1_num <- gene1_num / num
      gene2_num <- gene2_num / num
      overlap <- overlap / num
      or <- or / num
      gene_num1_new <- gene_num1_new / num
      gene_num2_new <- gene_num2_new /num
      overlap_new <- overlap_new /num
      background <- background / num
      temp_row <- c(species[j], species[i], gene1_num, gene2_num, overlap, gene_num1_new, gene_num2_new, overlap_new, background, or)
      venn_data <- rbind(venn_data, temp_row)
      
    }
  }
}
return(venn_data)
}


get_overlap_gene <- function(){
  species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','chlamydomonas','cyanobacteria')
  species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','chlamydomonas','cyanobacteria')
  gene_num_list <- data.frame(0,0)
  colnames(gene_num_list) <- c("species","gene_num")
  venn_data <- data.frame(0,0,0,0,0)
  colnames(venn_data) <- c("species1","species2","num1","num2","overlap")
  for(i in 1:(length(species) - 1))
  {
    for(j in (i+1):length(species))
    {
      #Get orthologue gene in current 2 species
      current_file_name <- paste(paste(paste(paste("ContainedData/Plot_required_file/protein2gene_20210122/",species[i],sep=""),'__v__',sep=""), species[j],sep=""),".txt",sep="")
      #Determine whether the file name exists
      if (file.exists(current_file_name))
      {
        gene_list <- unique(read.table(current_file_name))
        #Get cycling gene in curent 2 species
        current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
        current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
        gene1 <- read.csv(current_file_gene1_name,header = T)
        gene2 <- read.csv(current_file_gene2_name,header = T)
        gene1_num <- 0
        gene2_num <- 0
        overlap <- 0
        num <- 0
        gene1_total_num <- c()
        gene2_total_num <- 0
        for(k in 1:ncol(gene1))
        {
          for(l in 1:ncol(gene2))
          {
            num <- num + 1
            if(which(species_sort == species[i]) < which(species_sort == species[j]))
            {
              current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
              current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
              gene1_num <- gene1_num + length(unique(current_gene1))
              gene2_num <- gene2_num + length(unique(current_gene2))
              
            }else
            {
              current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
              current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
              gene1_num <- gene1_num + length(unique(current_gene1))
              gene2_num <- gene2_num + length(unique(current_gene2))
              
            }
            overlap <- overlap + length(intersect(current_gene1, current_gene2))
            
          }
        }
        gene1_num <- gene1_num / num
        gene2_num <- gene2_num / num
        overlap <- overlap / num
        temp_row <- c(species[i], species[j], gene1_num, gene2_num, overlap)
        venn_data <- rbind(venn_data, temp_row)
        gene_num_list <- rbind(gene_num_list,c(species[i],length(unique(gene1_total_num))))
        gene_num_list <- rbind(gene_num_list,c(species[j],length(unique(gene2_total_num))))
      }else
      {
        current_file_name <- paste(paste(paste(paste("ContainedData/Plot_required_file/protein2gene_20210122/",species[j],sep=""),'__v__',sep=""), species[i],sep=""),".txt",sep="")
        gene_list <- unique(read.table(current_file_name))
        current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
        current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
        gene1 <- read.csv(current_file_gene1_name,header = T)
        gene2 <- read.csv(current_file_gene2_name,header = T)
        gene1_num <- 0
        gene2_num <- 0
        overlap <- 0
        num <- 0
        for(k in 1:ncol(gene1))
        {
          for(l in 1:ncol(gene2))
          {
            num <- num + 1
            if(which(species_sort == species[i]) > which(species_sort == species[j]))
            {
              current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V1
              current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V1
              gene1_num <- gene1_num + length(unique(current_gene1))
              gene2_num <- gene2_num + length(unique(current_gene2))
            }else
            {
              current_gene1 <- gene_list[gene_list$V1 %in% gene1[,k],]$V2
              current_gene2 <- gene_list[gene_list$V2 %in% gene2[,l],]$V2
              gene1_num <- gene1_num + length(unique(current_gene1))
              gene2_num <- gene2_num + length(unique(current_gene2))
            }
            overlap <- overlap + length(intersect(current_gene1, current_gene2))
            
          }
        }
        gene1_num <- gene1_num / num
        gene2_num <- gene2_num / num
        overlap <- overlap / num
        temp_row <- c(species[j], species[i], gene1_num, gene2_num, overlap)
        venn_data <- rbind(venn_data, temp_row)
        
      }
    }
  }
  return(venn_data)
}




















get_heatmap_function <- function(){
  plot_heatmap_list <- list()
  
  draw_heatmap_function <- function(exp_data_name, range1, range2, photo_name){
    #draw cycling gene expression heatmap
    #exp_data_name: name of expression file
    #range1: start column number of expression data
    #range2: end column number of expression data
    #photo_name: output file name
    exp_data<-read.csv(exp_data_name)
    exp_data <- exp_data[exp_data$meta2d_BH.Q < 0.05,]
    exp_data <- exp_data[order(exp_data$meta2d_phase),]
    data <- exp_data[,range1:range2]
    exp_data_scale<-apply(data,1,scale)
    exp_data_scale<-t(exp_data_scale)
    p <- as.ggplot(pheatmap(exp_data_scale,cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = F,show_rownames = F,legend=FALSE,border=NA,color = colorRampPalette(c("black","purple","purple","yellow","yellow","yellow"))(20)))
    plot_heatmap_list <<- c(plot_heatmap_list, list(p))
  }
  #draw heatmap of rhythmic expression in different species
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Human_blood_for_heatmap.csv',  19,28,"Human blood")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Human_skin_for_heatmap.csv',  19,22,"Human skin")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Mouse_LIV.csv',  24,35, "Mouse liver")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Mouse_KID.csv',  24,35, "Mouse kidney")
  
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Fly_old.csv',  24,35, "Fly old")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Fly_young.csv',  24,35, "Fly young")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Yeast_high.csv',  24,43,"Yeast high")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Yeast_low.csv',  24,47,"Yeast low")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Neurospora.csv',  24,35,"Neurospora")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Arabidopsis_SD LEAF.csv',  24,35,"Arabidopsis SD leaf")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Arabidopsis_SD M.csv',  24,35,"Arabidopsis SD M")
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Arabidopsis_SD VA.csv',  24,35,"Arabidopsis SD VA")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Chlamydomounas_for_heatmap.csv',  19,26,"Chlamydomounas")
  
  draw_heatmap_function('ContainedData/Meta2dResult/AllSpecies/Cyanobacteria_for_heatmap.csv',  19,28,"Cyanobacteria")
  return(plot_heatmap_list)
}





















get_cycling_gene_orthologue_distribution <- function(gene_name_table, gene_table){
  species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','cyanobacteria','chlamydomonas')
  species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','cyanobacteria','chlamydomonas')
  gene_num_list <- data.frame(0,0)
  colnames(gene_num_list) <- c("species","gene_num")
  
  venn_data <- data.frame(0,0,0,0,0)
  colnames(venn_data) <- c("species1","species2","num1","num2","overlap")
  notdup_gene <- list(c(),c(),c(),c(),c(),c(),c(),c())

  
  for(i in 1:(length(species)))
  {
    for(j in 1:length(species))
    {
      if(i != j)
      {
        current_file_name <- paste(paste(paste(paste("ContainedData/Plot_required_file/protein2gene_20210122/",species[i],sep=""),'__v__',sep=""), species[j],sep=""),".txt",sep="")
        if (file.exists(current_file_name))
        {
          gene_list <- unique(read.table(current_file_name))
          current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[i],sep=""),"_gene.csv",sep="")
          current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/circidian_gene/",species[j],sep=""),"_gene.csv",sep="")
          gene1 <- read.csv(current_file_gene1_name,header = T)
          gene2 <- read.csv(current_file_gene2_name,header = T)
          for(k in 1:length(gene_name_table[[i]])){
            if(gene_name_table[[i]][k] %in% gene_list[,1]){
              temp_gene_name <- as.character(gene_list[gene_list$V1 %in% gene_name_table[[i]][k],2])
              if(length(temp_gene_name) >= 1){
                flag <- 0
                for(m in 1:length(temp_gene_name)){
                  if(temp_gene_name[m] %in% gene_name_table[[j]]){
                    flag <- 1
                  }
                }
                if(flag == 1){
                  gene_table[[i]][k] <- gene_table[[i]][k] + 1
                }
                
              }
            }
          }
          for(k in 1:length(gene_name_table[[j]])){
            if(gene_name_table[[j]][k] %in% gene_list[,2]){
              temp_gene_name <- c()
              temp_gene_name <- gene_list[gene_list$V2 %in% gene_name_table[[j]][k],1]
              # for(l in 1:nrow(gene_list)){
              #   if(gene_list[l,1] == gene_name_table[[i]][k]){
              #     temp_gene_name[length(temp_gene_name) + 1] <- gene_list[l,1]
              #   }
              #   
              # }
              if(length(temp_gene_name) >= 1){
                flag <- 0
                for(m in 1:length(temp_gene_name)){
                  if(temp_gene_name[m] %in% gene_name_table[[i]]){
                    flag <- 1
                  }
                }
                if(flag == 1){
                  gene_table[[j]][k] <- gene_table[[j]][k] + 1
                }
                
              }
            }
          }
          
          
        }
      }
    }
  }      
  
  return(gene_table)
}



get_go_term_orthologue_function <- function(gene_name_table, gene_table){
  species <- c("Human","Mouse","Fly","Yeast","Neurospora","Arabidopsis","Green alga","Cyanobacteria")
  
  species_sort <- c('human','mouse','fly','arabidopsis','yeast','neurospora','cyanobacteria','chlamydomonas')
  #species <- c('human','mouse','fly','yeast','neurospora','arabidopsis','cyanobacteria','chlamydomonas')
  gene_num_list <- data.frame(0,0)
  colnames(gene_num_list) <- c("species","gene_num") 
  
  venn_data <- data.frame(0,0,0,0,0)
  colnames(venn_data) <- c("species1","species2","num1","num2","overlap")
  notdup_gene <- list(c(),c(),c(),c(),c(),c(),c(),c())
  
  for(i in 1:(length(species)))
  {
    for(j in 1:length(species))
    {
      if(i != j)
      {
        
        current_file_gene1_name <- paste(paste("ContainedData/Plot_required_file/tissue14/enrich/",species[i],sep=""),".csv",sep="")
        current_file_gene2_name <- paste(paste("ContainedData/Plot_required_file/tissue14/enrich/",species[j],sep=""),".csv",sep="")
        gene1 <- read.csv(current_file_gene1_name,header = T)
        gene2 <- read.csv(current_file_gene2_name,header = T)
        for(k in 1:length(gene_name_table[[i]])){
          #temp_gene_name <- c()
          temp_gene_name <- as.character(gene_name_table[[i]][k])
          # for(l in 1:nrow(gene_list)){
          #   if(gene_list[l,1] == gene_name_table[[i]][k]){
          #     temp_gene_name[length(temp_gene_name) + 1] <- gene_list[l,2]
          #   }
          #   
          # }
          if(length(temp_gene_name) >= 1){
            flag <- 0
            for(m in 1:length(temp_gene_name)){
              if(temp_gene_name[m] %in% gene_name_table[[j]]){
                flag <- 1
              }
            }
            if(flag == 1){
              gene_table[[i]][k] <- gene_table[[i]][k] + 1
            }
          }
        }

        # }
      }
    }
  }  
  return(gene_table)
}