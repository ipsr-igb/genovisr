library(dplyr)
library(meshr)
library(ggplot2)
library(maptools)

#ワーキングディレクトリを設定
setwd("/Users/okapi/Desktop/物置/")

#2021年の遺伝子型データを読み込み
data_72 <- read.csv("/Users/okapi/Desktop/修論/修論用genotype/NAなし/TP72_mcptaggr_gbscleanr_dosage_20230820.csv", row.names = 1) # Your file
rownames(data_72) <- paste0("21TP72_", sprintf("%03d",as.numeric(sub("sample", "",rownames(data_72)))))
data_72 <- data_72[order(rownames(data_72)),]
data_72 <- data_72[-c(93, 105, 184, 191, 192),]


data_73 <- read.csv("/Users/okapi/Desktop/修論/修論用genotype/NAなし/TP73_mcptaggr_gbscleanr_dosage_20230820.csv", row.names = 1) # Your file
rownames(data_73) <- paste0("21TP73_", sprintf("%03d",as.numeric(sub("sample", "",rownames(data_73)))))
data_73 <- data_73[order(rownames(data_73)),]
data_73 <- data_73[-c(179, 181, 184, 191, 192),]

#2022年の遺伝子型データを読み込み（22TPDPのTP27は除外）
data_22 <- read.csv("/Users/okapi/Desktop/修論/修論用genotype/NAなし/tp_mcptaggr_gbscleanr_dosage_filled_20230613.csv", row.names = 1) # Your file
data_22 <- data_22[-c(27, 40, 85),]

#両年のマーカー情報を抽出
marker_72 <- colnames(data_72)
markerdf_72 <- data.frame(marker = marker_72,
                          chr = sub("_.*", "", marker_72),
                          pos = as.numeric(sub(".*_", "", marker_72)))
marker_73 <- colnames(data_73)
markerdf_73 <- data.frame(marker = marker_73,
                          chr = sub("_.*", "", marker_73),
                          pos = as.numeric(sub(".*_", "", marker_73)))
marker_22 <- colnames(data_22)
markerdf_22 <- data.frame(marker = marker_22,
                          chr = sub("_.*", "", marker_22),
                          pos = as.numeric(sub(".*_", "", marker_22)))

################################################################################
#ここから2022年のマーカーの近傍に存在する2021年のマーカーを推定する
closemarker <- NULL
#01染色体から12染色体でループ
for (i in 1:12) {
  #実行する染色体を選択
  chr_i <- c(paste0("chr", sprintf("%02d", 1:12)))[i]
  tp72_i <- markerdf_72[markerdf_72$chr == chr_i,]  
  tp22_i <- markerdf_22[markerdf_22$chr == chr_i,]
  
  #当該染色体における、マーカーの物理位置を抽出
  tp72_pos_i <- tp72_i[,"pos"]
  closemarker_j <- NULL
  
  #2022年の当該染色体におけるマーカーの個数分ループさせる
  for (j in 1:length(rownames(tp72_i))) {
    #１つのマーカーと当該染色体に座乗する全てのマーカーとの物理地図を計算
    dis_j <- tp22_i[,"pos"] - tp72_pos_i[j] 
    #計算結果をもとに、当該マーカーと最も近傍に存在するマーカーを選択する
    closemarker_j <- c(closemarker_j, tp22_i[which.min(abs(dis_j))[1],"marker"])
  }
  #2022年のマーカーの近傍に存在する、2021のマーカーをベクトル化
  closemarker <- c(closemarker, closemarker_j)
}
#2022年のマーカー情報に2021年の近傍マーカーの情報を結合
markerdf_72 <- cbind(markerdf_72, closemarker)


#近傍マーカー間の距離を算出（後に最近傍マーカーが一致するマーカーを削除するため）
markerdf_72$dis <- as.numeric(sub(".*_", "", markerdf_72$closemarker))
markerdf_72$diff <- c(1, diff(markerdf_72$dis))

df_72 <- cbind(markerdf_72, t(data_72))
#2022年の遺伝子型データに推定した最近傍マーカーの情報を結合
newmarker <- NULL
for (i in 1:length(rownames(df_72))) {
  marker_i <- data_22[,as.character(closemarker[i])]
  newmarker <- cbind(newmarker, marker_i)
}
newmarker <- as.data.frame(t(newmarker))
colnames(newmarker) <- rownames(data_22)
rownames(newmarker) <- rownames(df_72)

#新しいマーカーの情報でデータフレームの情報を作り直し、
#2021年と2022年を結合した遺伝子型データファイル作成
df_fig <- cbind(df_72, newmarker)
#マーカーのbiningを行い、近傍マーカーが同一の領域を削除
df_fig <- df_fig[!df_fig$diff == 0,-(1:6)]


################################################################################
#同様の事柄を再度行い、TP73でも共通のマーカーを抽出する
data_22 <- as.data.frame(t(df_fig))
marker_22 <- colnames(data_22)
markerdf_22 <- data.frame(marker = marker_22,
                          chr = sub("_.*", "", marker_22),
                          pos = as.numeric(sub(".*_", "", marker_22)))

closemarker <- NULL
for (i in 1:12) {
  #実行する染色体を選択
  chr_i <- c(paste0("chr", sprintf("%02d", 1:12)))[i]
  tp73_i <- markerdf_73[markerdf_73$chr == chr_i,]  
  tp22_i <- markerdf_22[markerdf_22$chr == chr_i,]
  
  #当該染色体における、マーカーの物理位置を抽出
  tp73_pos_i <- tp73_i[,"pos"]
  closemarker_j <- NULL
  
  #2022年の当該染色体におけるマーカーの個数分ループさせる
  for (j in 1:length(rownames(tp73_i))) {
    #１つのマーカーと当該染色体に座乗する全てのマーカーとの物理地図を計算
    dis_j <- tp22_i[,"pos"] - tp73_pos_i[j] 
    #計算結果をもとに、当該マーカーと最も近傍に存在するマーカーを選択する
    closemarker_j <- c(closemarker_j, tp22_i[which.min(abs(dis_j))[1],"marker"])
  }
  #2022年のマーカーの近傍に存在する、2021のマーカーをベクトル化
  closemarker <- c(closemarker, closemarker_j)
}
#2022年のマーカー情報に2021年の近傍マーカーの情報を結合
markerdf_73 <- cbind(markerdf_73, closemarker)




#近傍マーカー間の距離を算出（後に最近傍マーカーが一致するマーカーを削除するため）
markerdf_73$dis <- as.numeric(sub(".*_", "", markerdf_73$closemarker))
markerdf_73$diff <- c(1, diff(markerdf_73$dis))

df_73 <- cbind(markerdf_73, t(data_73))
#2022年の遺伝子型データに推定した最近傍マーカーの情報を結合
newmarker <- NULL
for (i in 1:length(rownames(df_73))) {
  marker_i <- data_22[,as.character(closemarker[i])]
  newmarker <- cbind(newmarker, marker_i)
}
newmarker <- as.data.frame(t(newmarker))
colnames(newmarker) <- rownames(data_22)
rownames(newmarker) <- rownames(df_73)

#新しいマーカーの情報でデータフレームの情報を作り直し、
#2021年と2022年を結合した遺伝子型データファイル作成
df_fig <- cbind(df_73, newmarker)
#マーカーのbiningを行い、近傍マーカーが同一の領域を削除
df_fig <- df_fig[!df_fig$diff == 0,-(1:6)]
#データフレームの準備完了
################################################################################


################################################################################
#グラフ遺伝子型のdataframe作成

#marker情報を抽出
id <- rownames(df_fig)
chr <- sub("_.*", "", id)
pos <- as.numeric(sub(".*_", "", id))

#個体名情報から集団情報を抽出
pop <- sub("_.*", "", colnames(df_fig))
dp　<- !(1:length(pop) %in% c(which(pop == "21TP72", TRUE), which(pop == "21TP73", TRUE), which(pop == "TP102", TRUE)))
#集団毎に分割
df_72 <- df_fig[,which(pop == "21TP72", TRUE)]
df_73 <- df_fig[,which(pop == "21TP73", TRUE)]
df_102 <- df_fig[,which(pop == "TP102", TRUE)]
df_dp <- df_fig[,which(dp, TRUE)]

#集団単位でループさせてplot用のデータフレーム作出
df_bind <- NULL
for (input in c("21TP72", "21TP73", "22TP102", "22TPDP")) {
  
  #データフレームを選択
  if(input == "21TP72"){
    df_i <- df_72
    p_name <- "TP72"
  } else if(input == "21TP73"){
    df_i <- df_73
    p_name <- "TP73"
  } else if(input == "22TP102"){
    df_i <- df_102
    p_name <- "TP102"
  }else if(input == "22TPDP"){
    df_i <- df_dp
    p_name <- "TPDP"
  }else {
    stop("Invalid choice for input!")
  }

  #marker(row)毎に各plexの個体数を算出
  df_table <- apply(df_i, 1, function(x){
    return(table(factor(x, 0:4)))
    })
  #個体数の情報を割合に変換
  df_table <- apply(df_table, 2,function(x) {
    return (x / length(colnames(df_i)))
    })                #Check!!!!

df_table <- data.frame(ID = id, Chr = chr, Pos = pos, t(df_table))
colnames(df_table) <- c("ID",  "Chr",  "Pos", "quadruplex", "triplex", "duplex", "simplex", "nulliplex")

plex_order <- c("nulliplex","simplex","duplex","triplex","quadruplex")


df <- NULL
df <- as.data.frame(df)

for(i in c(paste0("chr", sprintf("%02d", 1:12)))){
  # data
  table_i <- filter(df_table,Chr == i)
  
  m_n <- as.numeric(nrow(table_i))
  m_n5 <- 5*m_n
  
  i_chr <- rep(i, times = m_n5)
  i_pos <- table_i %>% select(Pos)
  i_pos <- as.vector(t(i_pos))
  i_pos <- rep(i_pos, times = 5)
  i_pos <- (i_pos)/ as.numeric(rep("1000000",times = m_n5))
  i_dos0 <- table_i %>% select("nulliplex")
  i_dos1 <- table_i %>% select("simplex")
  i_dos2 <- table_i %>% select("duplex")
  i_dos3 <- table_i %>% select("triplex")
  i_dos4 <- table_i %>% select("quadruplex")
  i_dos <- c(t(i_dos4),t(i_dos3),t(i_dos2),t(i_dos1),t(i_dos0))
  i_dsrep <- c((rep("quadruplex", times = m_n)),(rep("triplex", times = m_n)),(rep("duplex", times = m_n)),(rep("simplex", times = m_n)),(rep("nulliplex", times = m_n)))
  table_i <- data.frame(Chr = i_chr,Pos = i_pos, Ratio = i_dos, Dosage = i_dsrep)
  
  last_pos <- table_i[as.numeric(nrow(table_i)),1]
  last_pos <- as.numeric(last_pos)
  
  df <- rbind(df, table_i)
}

 df$Population <- rep(p_name, length(rownames(df)))
 df_bind <- rbind(df_bind, df)
}
#dataframe 作成完了
################################################################################


################################################################################
#ここからplot

##################################################
#集団毎にplot
for (i in c("TP72", "TP73", "TP102", "TPDP")) {
  dfplot_i <- filter(df_bind, Population == i)
  p <- ggplot(dfplot_i, aes(x = Pos,y = Ratio,fill=factor(Dosage,plex_order))) + 
    geom_area(position="fill") +
    scale_fill_viridis_d() +
    scale_colour_viridis_d() +
    guides(fill="none") +
    scale_y_continuous(expand = c(0,0),
                       limits = c(0,1),
                       breaks = c(0.5, 1.0),
                       labels = c("0.5", "1.0")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(0,10,20,30,40),
                       labels = c("0","10M","20M","30M","40M")) +
    facet_grid(Population  ~ Chr, scales = "free", space = "free") +
    theme(panel.spacing.y = unit(0.12, "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text = element_text(size = 20),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 17))
  pdf(paste0(i, "_ChromosomeBind.pdf"), width = 25, height = 2)
  print(p)
  dev.off()
}
##################################################


##################################################
#4集団縦積みplot
p <- ggplot(df_bind, aes(x = Pos,y = Ratio,fill=factor(Dosage,plex_order))) + 
  geom_area(position="fill") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  guides(fill="none") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1),
                     breaks = c(0.5, 1.0),
                     labels = c("0.5", "1.0")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(0,10,20,30,40),
                     labels = c("0","10M","20M","30M","40M")) +
  facet_grid(Population  ~ Chr, scales = "free", space = "free") +
  theme(panel.spacing.y = unit(0.12, "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 17))
pdf("AllPopulation_ChromosomeBind.pdf", width = 25, height = 5.2)
print(p)
dev.off()
##################################################


##################################################
#4集団縦積み2段重ねplot作成用　
df_bind$Population <- factor(df_bind$Population,
                             levels = c("TP72", "TP73", "TP102", "TPDP"))
Chr01_06 <- c("chr01", "chr02", "chr03", "chr04", "chr05", "chr06")
Chr07_12 <- c("chr07", "chr08", "chr09", "chr10", "chr11", "chr12")
df_01_06 <- filter(df_bind, Chr == Chr01_06)
df_07_12 <- filter(df_bind, Chr == Chr07_12)

p1 <- ggplot(df_01_06, aes(x = Pos,y = Ratio,fill=factor(Dosage,plex_order))) +
  geom_area(position="fill") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  guides(fill="none") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1),
                     breaks = c(0.5, 1.0),
                     labels = c("0.5", "1.0")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(0,10,20,30,40),
                     labels = c("0","10M","20M","30M","40M")) +
  facet_grid(Population  ~ Chr, scales = "free", space = "free") +
  theme(panel.spacing.y = unit(0.12, "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 17))
pdf("Chromosome_bind_01_06.pdf", width = 20, height = 5.2)
print(p1)
dev.off()

p2 <- ggplot(df_07_12, aes(x = Pos,y = Ratio,fill=factor(Dosage,plex_order))) +
  geom_area(position="fill") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  guides(fill="none") +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1),
                     breaks = c(0.5, 1.0),
                     labels = c("0.5", "1.0")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(0,10,20,30,40),
                     labels = c("0","10M","20M","30M","40M")) +
  facet_grid(Population  ~ Chr, scales = "free", space = "free") +
  theme(panel.spacing.y = unit(0.12, "in"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 17))
pdf("Chromosome_bind07_12.pdf", width = 15.85, height = 5.2)
print(p2)
dev.off()
##################################################

