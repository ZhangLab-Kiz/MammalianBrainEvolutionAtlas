library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)


df_sign_gene <- read.csv("sign_gene.csv")
gene_id <- bitr(unique(df_sign_gene$gene.id), fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

disease_pathway <- enrichDO(gene          = gene_id$ENTREZID, 
              ont           = "HDO",   
              pvalueCutoff  = 1,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 1,
              readable      = TRUE)

write.csv(disease_pathway, "disease_pathway.csv")