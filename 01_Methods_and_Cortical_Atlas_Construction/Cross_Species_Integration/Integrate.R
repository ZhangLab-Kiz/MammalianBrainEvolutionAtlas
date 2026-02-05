library(dplyr)
library(Seurat)
library(SeuratObject)
library(anndata)
library(ggplot2)
library(reticulate)
library(harmony)
use_python("/bin/python")


spelist <- readLines("spe2.list")
Neuronlist <- list()
NonNeuronlist <- list()
for(i in 1:length(spelist)){
  name <- stringr::str_replace(spelist[i], "^\\d+\\.", "")
  if(dir.exists(paste0("../",spelist[i],"/",name,"_result"))){
    Neuronlist[[i]] <- readRDS(paste0("../",spelist[i],"/",name,"_result/","Neuron.RDS"))
    Neuronlist[[i]] <- RenameCells(Neuronlist[[i]], add.cell.id = spelist[i])
    NonNeuronlist[[i]] <- readRDS(paste0("../",spelist[i],"/",name,"_result/","NonNeuron.RDS"))
    NonNeuronlist[[i]] <- RenameCells(NonNeuronlist[[i]], add.cell.id = spelist[i])
    }
    else {
    warning("Directory does not exist: ", paste0("../",spelist[i], "/", name, "_result"))
    }
}

neuron <- merge(Neuronlist[[1]], Neuronlist[2:length(Neuronlist)])
nonneuron <- merge(NonNeuronlist[[1]], NonNeuronlist[2:length(NonNeuronlist)])
AllCells <- merge(neuron,nonneuron,merge.data = TRUE,merge.dr = TRUE)


saveRDS(neuron,"AllNeuron.RDS")
saveRDS(nonneuron,"AllNonNeuron.RDS")
saveRDS(AllCells,"AllCells.RDS")
