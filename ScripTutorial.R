#MINI TUTORIAL 

# ¿Como usar el paquete NMF para hacer un unsupervised clustering?

# 1) Establecemos el directorio de trabajo y cargamos el paquete correspondiente: 

directory <- setwd("Desktop/Bioinformatics/MiniTutoriales/")
library(NMF)
library(tidyverse)

set.seed(123456789) #Tambin establecemos un seed porque al ser un proceso 
                    #stocastico intentamos reducir la aleatoriedad

# 2) Cargamos los datos: la matriz con la que vamos a trabajar es una matriz de cuentas que se ha 
#obtenido de la TGCA. Son datos de RNA Seq de adenocarcinoma de pulmon en los que los pacientes 
# no presentan ninguna mutacion en genes oncogenicos comunes. 

# Recordemos que en esta matriz: 
  # 1) tiene que ser no negativa
  # 2) los genes tiene que estar en filas 
  # 3) los pacientes tienen que estar en columnas 

Matriz <- readRDS("Luad_Filtrado.rds") #Usamos RDS porque es mas facil para guardar y abrir trabajar 
                                       # con las funciones saveRDS y readRDS. 

dim(Matriz)
#[1] 17789    61

# 3) Nombramos variables que nos va a ser ultiles durante el proceso: 
     # Los valores de K con los que vamos a trabajar. 
     # El valor de K es el numero de grupos con los que queremos trabajar. 

NumeroCluster_inicial <- 2
NumeroCluster_final <- 10
clusteresTotales <- seq(NumeroCluster_inicial,NumeroCluster_final,1)

# 4) Opcion 1:  nuestra propia approach de validacion. 

# Primero calculamos el cophenetic value en cada valor de k en todas las muestras: 

i <- 1
#Donde vamos a guardar los resultados obtenidos: 

valorCorrelacionCofenetci_Total <- vector(mode="logical",length =length(clusteresTotales)) 
matricesConsenso <- list()
Resultados <- list()
MatrizH <- list()
MatrizW <- list()

#Como tenemos varios valores de K lo hacemos dentro de un bucle: 

for (i in seq_along(clusteresTotales)){
  cluster <- clusteresTotales[i] #Elegimos primero el 
  print(paste("Numero de Cluster",cluster,"en NMF para Lusc en Todas las Muestras"))
  Resultados[[i]] <-  nmf(Matriz,cluster,method="brunet",seed="random",nrun=100)
  valorCorrelacionCofenetci_Total[[i]] <- cophcor(Resultados[[i]])
  matricesConsenso[[i]] <- consensus(Resultados[[i]])
  MatrizH[[i]] <- coef(Resultados[[i]])
  MatrizW[[i]] <- basis(Resultados[[i]])
}

saveRDS(valorCorrelacionCofenetci_Total,"ValorCofenetic_Total_Luad.rds")
saveRDS(matricesConsenso,"MatricesConsenso_Total_Luad.rds")
saveRDS(MatrizW,"MatrizW_Total_Luad.rds")
saveRDS(MatrizH,"MatrizH_Total_Luad.rds")
saveRDS(Resultado,"Resultados_NMF_Total_Luad.rds")

ValorCofenetic_Total_Luad <- readRDS("ValorCofenetic_Total_Luad.rds")
matricesConsenso <- readRDS("MatricesConsenso_Total_Luad.rds")
MatrizW <- readRDS("MatrizW_Total_Luad.rds")
MatrizH <- readRDS("MatrizH_Total_Luad.rds")
ResultadosTotales <- readRDS("Resultados_NMF_Total_Luad.rds")


# Y ahora hacemos lo mismo que con diferentes subsamples para poder comprobar la robustes y la estabilidad del indice: 
  # Lo haemos con 100 subsamples para cada valor de k

i <- 1 #Son las filas
j <- 1 #Son las columnas
valorCorrelacionCofenetci_ParteMuestras <- matrix(0,nrow=length(clusteresTotales),ncol=100)
valorCorrelacionCofenetci_ParteMuestras <- data.frame(valorCorrelacionCofenetci_ParteMuestras)
rownames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("cluster",2:10)
colnames(valorCorrelacionCofenetci_ParteMuestras) <- paste0("NumeroVuelta",1:ncol(valorCorrelacionCofenetci_ParteMuestras))
Resultados <- list()

for (i in seq_along(clusteresTotales)){
  cluster <- clusteresTotales[i]
  print(paste("Cluster",cluster,"NMF, Lusc con Parte de las muestras"))
  for (j in seq(1,100,1)){
    subsample <- Matriz[sample(nrow(Matriz), round(nrow(Matriz) * 0.7) ), ]
    Resultados[[j]] <- nmf(subsample,cluster,method="brunet",seed="random",nrun=100)
    valorCorrelacionCofenetci_ParteMuestras[i,j] <- cophcor(Resultados[[j]]) #Guardamos en el cluster i en la posicion j
  }
}

saveRDS(valorCorrelacionCofenetci_ParteMuestras,"ValorCofenetic_Luad_SubSamples.rds")
saveRDS(Resultados,"Resultados_Lusc_SubSamples.rds")

valorCorrelacionCofenetci_ParteMuestras <- readRDS("ValorCofenetic_Parte_De_Las_Muestr_NMF_Luad.rds")
Resultados_Subsamples <- readRDS("Resultados_Parte_De_Las_Muestras_NMF_Luad.rds")

# 4.2) Realizamos los graficos correspondientes: 
NumeroCluster <- seq(2,7,1)
  # Los graficos realizados con todas las muestras: 

ValorCofenetic_Total_Luad <- round(ValorCofenetic_Total_Luad,digits=3)
ValorCofenetic_Total_Luad <- data.frame(NumeroCluster=NumeroCluster,Valores=ValorCofenetic_Total_Luad)

Grafico_Luad_Cofenetic_Total <- ggplot(data=ValorCofenetic_Total_Luad,aes(x=NumeroCluster,y=Valores))+
  geom_line(color="#D0E6A5",size=2,alpha=0.5)+
  geom_point(aes(x=NumeroCluster,y=Valores),color="#D0E6A5",size=1.5)+
  geom_text(aes(x=NumeroCluster,y=Valores,label=Valores))+
  scale_x_continuous(breaks =seq(2,10,1))+
  labs(title="Cophenetic Correlation Luad")
Grafico_Luad_Cofenetic_Total

  # Los graficos realizados con parte de las muestras: 
numeroClusteres <- factor(rep(2:7,each=100))
cofeneticValue_Luad <- readRDS("ValorCofenetic_Parte_De_Las_Muestr_NMF_Luad.rds")
cofeneticValue_Luad <- as.matrix(cofeneticValue_Luad)

cofenetic_filaUnica_Luad <- as.vector(cofeneticValue_Luad)
COFENETIC_LUAD <- data.frame(numeroClusteres=numeroClusteres,
                             valores=cofenetic_filaUnica_Luad)

GraficoCofenetic_Luad <- ggplot(data=COFENETIC_LUAD,aes(x=numeroClusteres,y=valores,col=numeroClusteres))+
  geom_violin()+
  geom_jitter(alpha=0.3)+
  ylab("Indice Cofenetic")+
  xlab("Numero de Clusteres")+
  labs(title="Indice Cofenetic para Luad, subsamples")
GraficoCofenetic_Luad

    # Dibujamos las matrices consenso: 
# M A T R I C E S  C O N S E N S O  L U A D: 

library(NMF)
resultadosNMFLuad <- readRDS("Resultados_NMF_Total_Luad.rds")

K2 <- resultadosNMFLuad[[1]]
K3 <- resultadosNMFLuad[[2]]
K4 <- resultadosNMFLuad[[3]]
K5 <- resultadosNMFLuad[[4]]
K6 <- resultadosNMFLuad[[5]]
K7 <- resultadosNMFLuad[[6]]

colnames(K2@consensus) <- sampleNames(K2)
rownames(K2@consensus) <- sampleNames(K2)
pdf("K2_LUAD.pdf")
consensusmap(K2, hclustfun="average")
dev.off()

colnames(K3@consensus) <- sampleNames(K3)
rownames(K3@consensus) <- sampleNames(K3)
pdf("K3_LUAD.pdf")
consensusmap(K3, hclustfun="average")
dev.off()

colnames(K4@consensus) <- sampleNames(K4)
rownames(K4@consensus) <- sampleNames(K4)
pdf("K4_LUAD.pdf")
consensusmap(K4, hclustfun="average")
dev.off()

colnames(K5@consensus) <- sampleNames(K5)
rownames(K5@consensus) <- sampleNames(K5)
pdf("K5_LUAD.pdf")
consensusmap(K5, hclustfun="average")
dev.off()

colnames(K6@consensus) <- sampleNames(K6)
rownames(K6@consensus) <- sampleNames(K6)
pdf("K6_LUAD.pdf")
consensusmap(K6, hclustfun="average")
dev.off()

colnames(K7@consensus) <- sampleNames(K7)
rownames(K7@consensus) <- sampleNames(K7)
pdf("K7_LUAD.pdf")
consensusmap(K7, hclustfun="average")
dev.off()


