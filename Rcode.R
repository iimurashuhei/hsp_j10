
#データハンドリング ----

#前処理に必要なパッケージ読み込み
library(tidyverse)

#ローデータ読み込み
hsps_t1 <- read_csv("hsps_t1.csv", na = c(".", "")) #"hsps_t1.csv"のローデータ読み込み
bisbas_t1 <- read_csv("bisbas_t1.csv", na = c(".", "")) #"bisbas_t1.csv"のローデータ読み込み
tipi_t1 <- read_csv("tipi_t1.csv", na = c(".", "")) #"tipi_t1.csv"のローデータ読み込み
stai_t1 <- read_csv("stai_t1.csv", na = c(".", "")) #"stai_t1.csv"のローデータ読み込み
panas_t1 <- read_csv("panas_t1.csv", na = c(".", "")) #"panas_t1.csv"のローデータ読み込み
cesd_t1 <- read_csv("cesd_t1.csv", na = c(".", "")) #"cesd_t1.csv"のローデータ読み込み
hsps_t2 <- read_csv("hsps_t2.csv", na = c(".", "")) #"hsps_t2.csv"のローデータ読み込み

head(hsps_t1) #先頭6行確認
head(bisbas_t1) #先頭6行確認
head(tipi_t1) #先頭6行確認
head(stai_t1) #先頭6行確認
head(panas_t1) #先頭6行確認
head(cesd_t1) #先頭6行確認
head(hsps_t2) #先頭6行確認


#各ローデータの結合
data1 <- full_join(hsps_t1, bisbas_t1, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)
data2 <- full_join(data1, tipi_t1, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)
data3 <- full_join(data2, stai_t1, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)
data4 <- full_join(data3, panas_t1, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)
data5 <- full_join(data4, cesd_t1, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)
data6 <- full_join(data5, hsps_t2, by = "ID", na = c(".", "")) #full_join(データ1, データ2, by = 共通のキー変数名)

#csvで書き出し
write.csv(data6, file = "lowdata.csv", na = ".") #csvで書き出し


#結合済みローデータ読み込み
lowdata <- read_csv("lowdata.csv", na = c(".", "")) #ローデータ読み込み


# Study1: Exploratory Factor Analysis ----

#データ読み込み
library(tidyverse)
data_s1 <- read_csv("study1_rawdata.csv", na = c(".", "")) #"study1_rawdata.csv"のローデータ読み込み

#平行分析
library(psych)
library(MASS)
library(GPArotation)
names(data_s1)
fa.parallel(data_s1[,c(10:28)], fm = "ml", fa = "fa") #fm =最尤法推定, fa=主因子法抽出
  ##Parallel analysis suggests that the number of factors =  7  and the number of components =  NA 

#MAPテスト
VSS(data_s1[,c(10:28)], n = 10) #10因子まで推定
  ##The Velicer MAP achieves a minimum of 0.02  with  1  factors 

#固有値の算出
cor01 <- cor(data_s1[,c(10:28)])#相関係数算出
eigen01 <- eigen(cor01)$values　#固有値算出
eigen01
plot(eigen01, type="b", main="Scree Plot",xlab="Number", ylab="Eigenvalue")
  ## [1] 6.8014830 1.8076695 1.3126223 1.0392504 0.9693626 

#HSPS-J19の内的整合性
alpha(data_s1[, c(10:28)]) #alpha=.89
omega(data_s1[, c(10:28)],3,fm="ml") #omega hierarchical=.76

#EOEの内的整合性
alpha(data_s1[, c(12,15,17,20,21,24,25,26)]) #alpha=.83

#LSTの内的整合性
alpha(data_s1[, c(10,13,14,16,19,23,28)]) #alpha=.84

#AESの内的整合性
alpha(data_s1[, c(11,18,22,27)]) #alpha=.57


#3因子解のEFA（1回目）
factor3 <- fa(data_s1[, c(10:28)], nfactors = 3, rotate = "promax", fm = "ml", scores = TRUE)
print(factor3, sort = TRUE)

#3因子解のEFA（2回目）
factor3 <- fa(data_s1[, c(14,16:17,19,21:22,24:27)], nfactors = 3, rotate = "promax", fm = "ml", scores = TRUE)
print(factor3, sort = TRUE)

#1因子解のEFA（1回目）
factor1 <- fa(data_s1[, c(10:28)], nfactors = 1, rotate = "no", fm = "ml", scores = TRUE)
print(factor1, sort = TRUE)

#HSPS-J10の内的整合性
omega(data_s1[, c(14,16,17,19,21,22,24,25,26,27)],3,fm="ml") #alpha = 0.82

#EOE-J10の内的整合性
omega(data_s1[, c(24,25,26,17,21)],1,fm="ml") #alpha = 0.83

#LST-J10の内的整合性
omega(data_s1[, c(14,16,19)],1,fm="ml") #alpha = 0.79

#AES-J10の内的整合性
omega(data_s1[, c(27,22)],1,fm="ml") #alpha = 0.72


# Study2: Comfirmatory Factor Analysis and Convergent Validity ----

#データ読み込み
#ID0094を削除したデータを使用
library(tidyverse)
data_s2 <- read_csv("study2_rawdata.csv", na = c(".", "")) #"study2_rawdata.csv"のローデータ読み込み
names(data_s2)


#S2_1. データハンドリング -----

#HSPS-J10の下位尺度得点作成
data_s2 <- data_s2 %>% 
  dplyr::mutate(eoe_T1 = (hsp8_T1 + hsp11_T1 + hsp9_T1 + hsp15_T1 + hsp12_T1)/5, na.rm = TRUE) %>% #EOE平均得点
  dplyr::mutate(lst_T1 = (hsp5_T1 + hsp2_T1 + hsp1_T1)/3, na.rm = TRUE) %>% #LST平均得点
  dplyr::mutate(aes_T1 = (hsp18_T1 + hsp16_T1)/2, na.rm = TRUE) %>% #AES平均得点
  dplyr::mutate(hsp_T1 = (eoe_T1 + lst_T1 + aes_T1)/3, na.rm = TRUE) %>% #HSPS-J10平均得点
  dplyr::mutate(eoe_T2 = (hsp8_T2 + hsp11_T2 + hsp9_T2 + hsp15_T2 + hsp12_T2)/5, na.rm = TRUE) %>% #EOE平均得点
  dplyr::mutate(lst_T2 = (hsp5_T2 + hsp2_T2 + hsp1_T2)/3, na.rm = TRUE) %>% #LST平均得点
  dplyr::mutate(aes_T2 = (hsp18_T2 + hsp16_T2)/2, na.rm = TRUE) %>% #AES平均得点
  dplyr::mutate(hsp_T2 = (eoe_T2 + lst_T2 + aes_T2)/3, na.rm = TRUE) %>% #HSPS-J10平均得点
  dplyr::mutate(bis_T1 = BIS_total_T1/7, na.rm = TRUE) %>% #BIS平均得点
  dplyr::mutate(bas_T1 = BAS_total_T1/13, na.rm = TRUE) %>% #BAS平均得点
  dplyr::mutate(E_T1 = E_total_T1/2, na.rm = TRUE) %>% #E平均得点
  dplyr::mutate(A_T1 = A_total_T1/2, na.rm = TRUE) %>% #A平均得点
  dplyr::mutate(C_T1 = C_total_T1/2, na.rm = TRUE) %>% #C平均得点
  dplyr::mutate(N_T1 = N_total_T1/2, na.rm = TRUE) %>% #N平均得点
  dplyr::mutate(O_T1 = O_total_T1/2, na.rm = TRUE) %>% #O平均得点
  dplyr::mutate(STAI_T1 = STAI_total_T1/20, na.rm = TRUE) %>% #STAI平均得点
  dplyr::mutate(ces_T1 = CES_D_total_T1/20, na.rm = TRUE) %>% #CES-D平均得点
  dplyr::mutate(NA_T1 = NA_total_T1/8, na.rm = TRUE) %>% #NA平均得点
  dplyr::mutate(PA_T1 = PA_total_T1/8, na.rm = TRUE) %>% #NA平均得点
  dplyr::select(-na.rm) %>% #na.rmという謎の変数が勝手に作成されるので除外
  dplyr::select(-hsp_total_T1) %>% #HSPS-J19_T1の変数を除外
  dplyr::select(-hsp_total_T2) #HSPS-J19_T2の変数を除外
names(data_s2) #変数名確認

#S2_2. 度数分布とヒストグラム ----

#genderの度数分布とヒストグラム
gender <- dplyr::count(data_s2, Gender_T1)
knitr::kable(gender) #テーブル化
ggplot(data = data_s2, mapping = aes(x = Gender_T1, fill = factor(Gender_T1))) + geom_bar() #視覚化

#ageの度数分布とヒストグラム
age <- dplyr::count(data_s2, age_T1)
knitr::kable(age) #テーブル化
ggplot(data = data_s2, mapping = aes(x = age_T1, fill = factor(age_T1))) + geom_bar()  + guides(fill = "none") #視覚化

#HSP-J10_T1平均得点の度数分布とヒストグラム
hsp <- dplyr::count(data_s2, hsp_T1)
knitr::kable(hsp) #テーブル化
ggplot(data = data_s2, mapping = aes(x = hsp_T1, fill = factor(hsp_T1))) + 
  geom_histogram(binwidth = 0.3) + guides(fill = "none") #視覚化

#EOE_T1平均得点の度数分布とヒストグラム
eoe <- dplyr::count(data_s2, eoe_T1)
knitr::kable(eoe) #テーブル化
ggplot(data = data_s2, mapping = aes(x = eoe_T1, fill = factor(eoe_T1))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#LST_T1平均得点の度数分布とヒストグラム
lst <- dplyr::count(data_s2, lst_T1)
knitr::kable(lst) #テーブル化
ggplot(data = data_s2, mapping = aes(x = lst_T1, fill = factor(lst_T1))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#AES_T1平均得点の度数分布とヒストグラム
aes <- dplyr::count(data_s2, aes_T1)
knitr::kable(aes) #テーブル化
ggplot(data = data_s2, mapping = aes(x = aes_T1, fill = factor(aes_T1))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#BIS_total_T1合計得点の度数分布とヒストグラム
BIS_total_T1 <- dplyr::count(data_s2, BIS_total_T1)
knitr::kable(BIS_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = BIS_total_T1, fill = factor(BIS_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#BAS_total_T1合計得点の度数分布とヒストグラム
BAS_total_T1 <- dplyr::count(data_s2, BAS_total_T1)
knitr::kable(BAS_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = BAS_total_T1, fill = factor(BAS_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#E_total_T1合計得点の度数分布とヒストグラム
E_total_T1 <- dplyr::count(data_s2, E_total_T1)
knitr::kable(E_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = E_total_T1, fill = factor(E_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#A_total_T1合計得点の度数分布とヒストグラム
A_total_T1 <- dplyr::count(data_s2, A_total_T1)
knitr::kable(A_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = A_total_T1, fill = factor(A_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#C_total_T1合計得点の度数分布とヒストグラム
C_total_T1 <- dplyr::count(data_s2, C_total_T1)
knitr::kable(C_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = C_total_T1, fill = factor(C_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#N_total_T1合計得点の度数分布とヒストグラム
N_total_T1 <- dplyr::count(data_s2, N_total_T1)
knitr::kable(N_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = N_total_T1, fill = factor(N_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#O_total_T1合計得点の度数分布とヒストグラム
O_total_T1 <- dplyr::count(data_s2, O_total_T1)
knitr::kable(O_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = O_total_T1, fill = factor(O_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#STAI_total_T1合計得点の度数分布とヒストグラム
STAI_total_T1 <- dplyr::count(data_s2, STAI_total_T1)
knitr::kable(STAI_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = STAI_total_T1, fill = factor(STAI_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#NA_total_T1合計得点の度数分布とヒストグラム
NA_total_T1 <- dplyr::count(data_s2, NA_total_T1)
knitr::kable(NA_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = NA_total_T1, fill = factor(NA_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#PA_total_T1合計得点の度数分布とヒストグラム
PA_total_T1 <- dplyr::count(data_s2, PA_total_T1)
knitr::kable(PA_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = PA_total_T1, fill = factor(PA_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#CES_D_total_T1合計得点の度数分布とヒストグラム
CES_D_total_T1 <- dplyr::count(data_s2, CES_D_total_T1)
knitr::kable(CES_D_total_T1) #テーブル化
ggplot(data = data_s2, mapping = aes(x = CES_D_total_T1, fill = factor(CES_D_total_T1))) + 
  geom_histogram(binwidth = 1.0) + guides(fill = "none") #視覚化

#S2_3. 記述統計量 ----

#HSPS_T1の記述統計量
hsp_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(hsp.mean.T1 = mean (hsp_T1, na.rm = TRUE), #HSPS-J10_T1の平均値
                   hsp.sd.T1 = sd (hsp_T1, na.rm = TRUE), #HSPS-J10_T1のSD
                   eoe.mean.T1 = mean (eoe_T1, na.rm = TRUE), #EOEの平均値
                   eoe.sd.T1 = sd (eoe_T1, na.rm = TRUE), #EOEのSD
                   lst.mean.T1 = mean (lst_T1, na.rm = TRUE), #LSTの平均値
                   lst.sd.T1 = sd (lst_T1, na.rm = TRUE), #LSTのSD
                   aes.mean.T1 = mean (aes_T1, na.rm = TRUE), #AESの平均値
                   aes.sd.T1 = sd (aes_T1, na.rm = TRUE)) #EOEのSD
knitr::kable(hsp_t1_stat, digits = 2) #出力

#BISBASの記述統計量
bisbas_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(bis.mean.T1 = mean (bis_T1, na.rm = TRUE), #bis_T1の平均値
                   bis.sd.T1 = sd (bis_T1, na.rm = TRUE), #bis_T1のSD
                   bas.mean.T1 = mean (bas_T1, na.rm = TRUE), #bas_T1の平均値
                   bas.sd.T1 = sd (bas_T1, na.rm = TRUE)) #bas_T1のSD
knitr::kable(bisbas_t1_stat, digits = 2) #出力

#TIPIの記述統計量
tipi_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(E.mean.T1 = mean (E_T1, na.rm = TRUE), #E_T1の平均値
                   E.sd.T1 = sd (E_T1, na.rm = TRUE), #E_T1のSD
                   A.mean.T1 = mean (A_T1, na.rm = TRUE), #A_T1の平均値
                   A.sd.T1 = sd (A_T1, na.rm = TRUE), #A_T1のSD
                   C.mean.T1 = mean (C_T1, na.rm = TRUE), #C_T1の平均値
                   C.sd.T1 = sd (C_T1, na.rm = TRUE), #C_T1のSD
                   N.mean.T1 = mean (N_T1, na.rm = TRUE), #N_T1の平均値
                   N.sd.T1 = sd (N_T1, na.rm = TRUE), #N_T1のSD
                   O.mean.T1 = mean (O_T1, na.rm = TRUE), #O_T1の平均値
                   O.sd.T1 = sd (O_T1, na.rm = TRUE)) #O_T1のSD
knitr::kable(tipi_t1_stat, digits = 2) #出力

#STAIの記述統計量
stai_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(STAI.mean.T1 = mean (STAI_T1, na.rm = TRUE), #STAI_T1の平均値
                   STAI.sd.T1 = sd (STAI_T1, na.rm = TRUE)) #STAI_T1のSD
knitr::kable(stai_t1_stat, digits = 2) #出力

#PANASの記述統計量
panas_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(pa.mean.T1 = mean (PA_T1, na.rm = TRUE), #PA_T1の平均値
                   pa.sd.T1 = sd (PA_T1, na.rm = TRUE), #PA_T1のSD
                   na.mean.T1 = mean (NA_T1, na.rm = TRUE), #NA_T1の平均値
                   na.sd.T1 = sd (NA_T1, na.rm = TRUE)) #NA_T1のSD
knitr::kable(panas_t1_stat, digits = 2) #出力

#CES_Dの記述統計量
ces_t1_stat <- 
  data_s2 %>% 
  dplyr::summarise(ces.mean.T1 = mean (ces_T1, na.rm = TRUE), #CES_T1の平均値
                   ces.sd.T1 = sd (ces_T1, na.rm = TRUE)) #CES_T1のSD
knitr::kable(ces_t1_stat, digits = 2) #出力


#S2_4. 内的整合性 ----

library(psych)
library(GPArotation)

#HSPS-J10_T1の内的整合性
omega(data_s2[, c(23,24,26,27,30,16,17,20,31,33)],3,fm="ml") #alpha = .86

#EOE_T1の内的整合性
omega(data_s2[, c(23,24,26,27,30)],1,fm="ml") #alpha = .85

#LST_T1の内的整合性
omega(data_s2[, c(16,17,20)],1,fm="ml") #alpha = .76

#AES_T1の内的整合性
omega(data_s2[, c(31,33)],1,fm="ml") #alpha = .60

#BIS_T1の内的整合性
omega(data_s2[, c(36,41,45,48,50,54,56)],1,fm="ml") #alpha = .82

#BAS_T1の内的整合性
omega(data_s2[, c(37,38,39,40,42,43,44,46,47,49,51,52,55)],1,fm="ml") #alpha = .88

#E_T1の内的整合性
omega(data_s2[, c(57,64)],1,fm="ml") #alpha = .66
cor.test(data_s2$TIPI1_T1, data_s2$TIPI6_T1_t, na.rm=TRUE) #r=0.4893781 (p<.001) #2項目相関

#A_T1の内的整合性
omega(data_s2[, c(59,65)],1,fm="ml") #alpha = .28
cor.test(data_s2$TIPI2_T1_t, data_s2$TIPI7_T1, na.rm=TRUE) #r=0.163 (p<.001)　#2項目相関

#C_T1の内的整合性
omega(data_s2[, c(60,67)],1,fm="ml") #alpha = .52
cor.test(data_s2$TIPI3_T1, data_s2$TIPI8_T1_t, na.rm=TRUE) #r=0.349 (p<.001)　#2項目相関

#N_T1の内的整合性
omega(data_s2[, c(61,69)],1,fm="ml") #alpha = .42
cor.test(data_s2$TIPI4_T1, data_s2$TIPI9_T1_t, na.rm=TRUE) #r=0.267 (p<.001)　#2項目相関

#O_T1の内的整合性
omega(data_s2[, c(62,71)],1,fm="ml") #alpha = .48
cor.test(data_s2$TIPI5_T1, data_s2$TIPI10_T1_t, na.rm=TRUE) #r=0.3166 (p<.001)　#2項目相関

#STAI_T1の内的整合性
omega(data_s2[, c(73:77,79,81:83,85:87,90,91,93:95,97,98)],1,fm="ml") #alpha = .92

#PA_T1の内的整合性
omega(data_s2[, c(100,102,104,106,108,110,112,114)],1,fm="ml") #alpha = .89

#NA_T1の内的整合性
omega(data_s2[, c(99,101,103,105,107,109,111,113)],1,fm="ml") #alpha = .92

#CES_T1の内的整合性
omega(data_s2[, c(115:117,119:122,124:127,129:132,134:138)],1,fm="ml") #alpha = .90


#S2_5. 相関係数 ----

#Rだと面倒なのでHADで出力する
#csvで書き出し
write.csv(data_s2, file = "data_for_cor.csv", na = ".") #csvで書き出し

#結果の詳細はHADを参照


#S2_6. 確認的因子分析（3因子解） ----
library(lavaan)
library(semPlot)
library(semTools)

model_3f <-'
EOE =~ hsp8_T1 + hsp9_T1 + hsp11_T1 + hsp12_T1 + hsp15_T1
LST =~ hsp1_T1 + hsp2_T1 + hsp5_T1
AES =~ hsp16_T1 + hsp18_T1
'

fit_3f <- cfa(model_3f, data = data_s2, missing="fiml") 
summary(fit_3f, fit.measures = TRUE, standardized = TRUE)
plot_3f <- semPaths(fit_3f, "std", edge.label.cex=.8, fade = FALSE, gray = TRUE, mar=c(6,1,3,1), style="lisrel") #作図

#S2_7. 確認的因子分析（Bifactor解） ----

model_bi <-'
EOE =~ 1*hsp8_T1 + hsp9_T1 + hsp11_T1 + hsp12_T1 + hsp15_T1
LST =~ hsp1_T1 + hsp2_T1 + hsp5_T1
AES =~ hsp16_T1 + hsp18_T1
HSP =~ hsp8_T1 + hsp9_T1 + hsp11_T1 + hsp12_T1 + hsp15_T1 + hsp1_T1 + hsp2_T1 + hsp5_T1 + hsp16_T1 + hsp18_T1

#各因子間を0に固定
HSP ~~ 0*EOE
HSP ~~ 0*LST
HSP ~~ 0*AES
EOE ~~ 0*LST
EOE ~~ 0*AES
LST ~~ 0*AES

#各潜在変数の分散を1に固定
HSP ~~ 1*HSP
EOE ~~ 1*EOE
LST ~~ 1*LST
AES ~~ 1*AES
'

fit_bi <- cfa(model_bi, data = data_s2, missing = "fiml")
summary(fit_bi, fit.measures = TRUE, standardized = TRUE)
plot_bi <- semPaths(fit_bi, "std", edge.label.cex=.8, fade = FALSE, gray = TRUE, mar=c(6,1,3,1), style="lisrel") #作図


# Study3: Longitudinal CFA and Test-retest Validity and Convergent Validity ----

#S3_1. データハンドリング ----

data_s3 <- data_s2 #一応、Study3で使用するデータ名をdata_s3に変更

#S3_2. 度数分布とヒストグラム ----

#HSP-J10_T2平均得点の度数分布とヒストグラム
hsp <- dplyr::count(data_s3, hsp_T2)
knitr::kable(hsp) #テーブル化
ggplot(data = data_s3, mapping = aes(x = hsp_T2, fill = factor(hsp_T2))) + 
  geom_histogram(binwidth = 0.3) + guides(fill = "none") #視覚化

#EOE_T1平均得点の度数分布とヒストグラム
eoe <- dplyr::count(data_s3, eoe_T2)
knitr::kable(eoe) #テーブル化
ggplot(data = data_s3, mapping = aes(x = eoe_T2, fill = factor(eoe_T2))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#LST_T1平均得点の度数分布とヒストグラム
lst <- dplyr::count(data_s3, lst_T2)
knitr::kable(lst) #テーブル化
ggplot(data = data_s3, mapping = aes(x = lst_T2, fill = factor(lst_T2))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#AES_T1平均得点の度数分布とヒストグラム
aes <- dplyr::count(data_s3, aes_T2)
knitr::kable(aes) #テーブル化
ggplot(data = data_s3, mapping = aes(x = aes_T2, fill = factor(aes_T2))) + 
  geom_histogram(binwidth = 0.5) + guides(fill = "none") #視覚化

#S3_3. 記述統計量 ----

#HSPS_T2の記述統計量
hsp_t2_stat <- 
  data_s2 %>% 
  dplyr::summarise(hsp.mean.T2 = mean (hsp_T2, na.rm = TRUE), #HSPS-J10_T2の平均値
                   hsp.sd.T2 = sd (hsp_T2, na.rm = TRUE), #HSPS-J10_T2のSD
                   eoe.mean.T2 = mean (eoe_T2, na.rm = TRUE), #EOEの平均値
                   eoe.sd.T2 = sd (eoe_T2, na.rm = TRUE), #EOEのSD
                   lst.mean.T2 = mean (lst_T2, na.rm = TRUE), #LSTの平均値
                   lst.sd.T2 = sd (lst_T2, na.rm = TRUE), #LSTのSD
                   aes.mean.T2 = mean (aes_T2, na.rm = TRUE), #AESの平均値
                   aes.sd.T2 = sd (aes_T2, na.rm = TRUE)) #EOEのSD
knitr::kable(hsp_t2_stat, digits = 2) #出力


#S3_4. 内的整合性 ----

#HSPS-J10_T2の内的整合性
omega(data_s2[, c(146,147,149,150,153,139,140,143,154,156)],3,fm="ml") #alpha = .84

#EOE_T2の内的整合性
omega(data_s2[, c(146,147,149,150,153)],1,fm="ml") #alpha = .86

#LST_T2の内的整合性
omega(data_s2[, c(139,140,143)],1,fm="ml") #alpha = .73

#AES_T2の内的整合性
omega(data_s2[, c(154,156)],1,fm="ml") #alpha = .58

#S3_5. 再検査信頼性（級内相関係数ICC）-----

#irrパッケージ読み込み
library(irr)

#ICCに必要な変数だけのデータセットを作成
icc_hsp <- data_s3[,c(161,165)]
icc_eoe <- data_s3[,c(158,162)]
icc_lst <- data_s3[,c(159,163)]
icc_aes <- data_s3[,c(160,164)]

#HSCS得点のICC
icc(icc_hsp, "twoway", "agreement") #ICC(A,1) = 0.711, F(278,279) = 5.92 , p = 9.05e-45, 0.648 < ICC < 0.765
icc(icc_eoe, "twoway", "agreement") #ICC(A,1) = 0.705, F(278,279) = 5.79 , p = 7.44e-44, 0.641 < ICC < 0.759
icc(icc_lst, "twoway", "agreement") #ICC(A,1) = 0.628, F(278,279) = 4.38 , p = 1.47e-32, 0.552 < ICC < 0.694
icc(icc_aes, "twoway", "agreement") #ICC(A,1) = 0.651, F(278,278) = 4.75 , p = 1.22e-35, 0.578 < ICC < 0.714

#S3_6. 再検査信頼性（相関係数r）-----

cor.test(icc_hsp$hsp_T1, icc_hsp$hsp_T2, na.rm = TRUE) # r = 0.7109667, t = 16.827, df = 277, p-value < 2.2e-16, 95CI 0.6476035 0.7645649
cor.test(icc_eoe$eoe_T1, icc_eoe$eoe_T2, na.rm = TRUE) # r = 0.7056287, t = 16.574, df = 277, p-value < 2.2e-16, 95CI 0.6413406 0.7600779
cor.test(icc_lst$lst_T1, icc_lst$lst_T2, na.rm = TRUE) # r = 0.6285489, t = 13.450, df = 277, p-value < 2.2e-16, 95CI 0.5518502 0.6947034
cor.test(icc_aes$aes_T1, icc_aes$aes_T2, na.rm = TRUE) # r = 0.6541016, t = 14.392, df = 277, p-value < 2.2e-16, 65CI 0.5813226 0.7164975

#S3_7. 確認的因子分析(3因子解) ----

library(lavaan)
library(semPlot)
library(semTools)

model_3f <-'
EOE =~ hsp8_T2 + hsp9_T2 + hsp11_T2 + hsp12_T2 + hsp15_T2
LST =~ hsp1_T2 + hsp2_T2 + hsp5_T2
AES =~ hsp16_T2 + hsp18_T2
'

fit_3f <- cfa(model_3f, data = data_s3, missing = "fiml") 
summary(fit_3f, fit.measures = TRUE, standardized = TRUE)
plot_3f <- semPaths(fit_3f, "std", edge.label.cex=.8, fade = FALSE, gray = TRUE, mar=c(6,1,3,1), style="lisrel") #作図
#分散が負になり不適解

#S3_8. 確認的因子分析（Bifactor解） ----

model_bi <-'
EOE =~ 1*hsp8_T2 + hsp9_T2 + hsp11_T2 + hsp12_T2 + hsp15_T2
LST =~ hsp1_T2 + hsp2_T2 + hsp5_T2
AES =~ hsp16_T2 + hsp18_T2
HSP =~ hsp8_T2 + hsp9_T2 + hsp11_T2 + hsp12_T2 + hsp15_T2 + hsp1_T2 + hsp2_T2 + hsp5_T2 + hsp16_T2 + hsp18_T2

#各因子間を0に固定
HSP ~~ 0*EOE
HSP ~~ 0*LST
HSP ~~ 0*AES
EOE ~~ 0*LST
EOE ~~ 0*AES
LST ~~ 0*AES

#各潜在変数の分散を1に固定
HSP ~~ 1*HSP
EOE ~~ 1*EOE
LST ~~ 1*LST
AES ~~ 1*AES
'

fit_bi <- cfa(model_bi, data = data_s3, missing = "fiml")
summary(fit_bi, fit.measures = TRUE, standardized = TRUE)
plot_bi <- semPaths(fit_bi, "std", edge.label.cex=.8, fade = FALSE, gray = TRUE, mar=c(6,1,3,1), style="lisrel") #作図

#適合度すべて良好


#S3_9. 縦断的測定不変性 -----

#因子分析に必要な変数だけを抽出したデータセットを作成
data_cfa <- data_s3 %>% select(1,23,24,26,27,30,16,17,20,31,33,146,147,149,150,153,139,140,143,154,156) 
names(data_cfa)

### long型にデータ変換

longdata <- data_cfa %>% 
  gather(key = "time", value = "hsp8", hsp8_T1, hsp8_T2) %>% 
  gather(key = "jikan1", value = "hsp9", hsp9_T1, hsp9_T2) %>% 
  gather(key = "jikan2", value = "hsp11", hsp11_T1, hsp11_T2) %>% 
  gather(key = "jikan3", value = "hsp12", hsp12_T1, hsp12_T2) %>% 
  gather(key = "jikan4", value = "hsp15", hsp15_T1, hsp15_T2) %>% 
  gather(key = "jikan5", value = "hsp1", hsp1_T1, hsp1_T2) %>% 
  gather(key = "jikan7", value = "hsp2", hsp2_T1, hsp2_T2) %>% 
  gather(key = "jikan8", value = "hsp5", hsp5_T1, hsp5_T2) %>% 
  gather(key = "jikan9", value = "hsp16", hsp16_T1, hsp16_T2) %>% 
  gather(key = "jikan10", value = "hsp18", hsp18_T1, hsp18_T2)

#面倒だけど、重複するIDは手作業で削除して整形
#csvで書き出し
write.csv(longdata, file = "longdata.csv", na = ".") #csvで書き出し

#整形したデータを読み込み
longdata <- read_csv("longdata.csv", na = c(".", "")) #"longdata.csv"の読み込み
names(longdata)


###モデル記述（3因子モデル）
model_3f_long <-'
EOE =~ hsp8 + hsp9 + hsp11 + hsp12 + hsp15
LST =~ hsp1 + hsp2 + hsp5
AES =~ hsp16 + hsp18
'
measurementInvariance(model = model_3f_long, data = longdata, std.lv = TRUE, strict = TRUE, fit.measures = c("cfi", "rmsea", "aic"), group= "time", missing = "fiml")

config <- cfa(model_3f_long, data = longdata, group = "time", missing = "fiml")
summary(config, fit.measures = TRUE, standardized = TRUE) #列が長いから結構推定時間長くかかる


###モデル記述（bi-factorモデル）
model_bi_long <-'
EOE =~ 1*hsp8 + hsp9 + hsp11 + hsp12 + hsp15
LST =~ hsp1 + hsp2 + hsp5
AES =~ hsp16 + hsp18
HSP =~ hsp8 + hsp9 + hsp11 + hsp12 + hsp15 + hsp1 + hsp2 + hsp5 + hsp16 + hsp18

#各因子間を0に固定
HSP ~~ 0*EOE
HSP ~~ 0*LST
HSP ~~ 0*AES
EOE ~~ 0*LST
EOE ~~ 0*AES
LST ~~ 0*AES

#各潜在変数の分散を1に固定
HSP ~~ 1*HSP
EOE ~~ 1*EOE
LST ~~ 1*LST
AES ~~ 1*AES
'

measurementInvariance(model = model_bi_long, data = longdata, std.lv = TRUE, strict = TRUE, fit.measures = c("cfi", "rmsea", "aic"), group= "time", missing = "fiml")
#Δχ2検定、AIC、BICにもとづく結果、Model 5を支持
#Model 5: The factor loadings, intercepts, residual variances and means are constrained to be equal across groups.



#複数の因子がある場合はlongInvariance()が使用できないらしい。1因子構造しか提供していないらしい。https://stackoverflow.com/questions/41099149/lavaan-longitudinal-invariance-cfa-with-a-2-factor-model-in-r

#時点ごとに因子モデルを描く
#wb_cfa <- '
#wbt1 =~ health1_T1 + health2_T1 + health3_T1 + health4_T1 + health5_T1
#wbt2 =~ health1_T2 + health2_T2 + health3_T2 + health4_T2 + health5_T2
#wbt3 =~ health1_T3 + health2_T3 + health3_T3 + health4_T3 + health5_T3
#wbt4 =~ health1_T4 + health2_T4 + health3_T4 + health4_T4 + health5_T4
#'
##Create list of variables（時間ごとの変数名のリストをつくる）
#var1wb <- c("health1_T1","health2_T1", "health3_T1", "health4_T1", "health5_T1")#T1
#var2wb <- c("health1_T2","health2_T2", "health3_T2", "health4_T2", "health5_T2")#T2
#var3wb <- c("health1_T3","health2_T3", "health3_T3", "health4_T3", "health5_T3")#T3
#var4wb <- c("health1_T4","health2_T4", "health3_T4", "health4_T4", "health5_T4")#T4
#constrainedVarwb <- list(var1wb, var2wb, var3wb, var4wb)
#longInvariance(model = wb_cfa, auto = 1, data = DataSource, varList = constrainedVarwb, constrainAuto = TRUE, missing = "fiml")
