##################################################################################################################
############################################## create phylogenetic trees and use them for analysis


#
library("msa")
library("geiger")
library("ape")
library("seqinr")
library("caper")

# algin Cyt B DNA sequences, create a tree and clean it up
names(mammal_subset_toren_cytb_dnaset2)<-mammal_toren_anage_pacif_df_c2$Species

msa_muscle_phylo_mammal_600sp <- msa(mammal_subset_toren_cytb_dnaset2, "ClustalOmega")
msa_muscle_phylo_mammal_600sp_conv <- msaConvert(msa_muscle_phylo_mammal_600sp, type="seqinr::alignment")
dm4<-dm <- dist.alignment(msa_muscle_phylo_mammal_600sp_conv , "identity") ## takes seuqinr alignment
mammal_tree_sp600m_clustal <- nj(dm4)
min_tip_length<-mammal_tree_sp600m_clustal$edge.length[mammal_tree_sp600m_clustal$edge.length>0]
min_tip_length<-min(min_tip_length)
mammal_tree_sp600m_clustal$edge.length[which(mammal_tree_sp600m_clustal$edge.length==0)]<-min_tip_length

## prepare data for PGLS
mammal_data<-data.frame(mammal_triplex,log(mammal_toren_anage_pacif_df_c2$MLS_anage),mammal_toren_anage_pacif_df_c2$Species,
                        mammal_subset_toren_dnaset_alph2$alpha_GC,mammal_subset_toren_dnaset_alph2$alpha_GC_skew,
                        mammal_subset_toren_dnaset_alph2$alpha_AT_skew,log(mammal_toren_anage_pacif_df_c2$adult_weight_anage))
names(mammal_data)<-c("trip","MLS_Log","Species","GC","GC_skew","AT_skew", "BM_Log")
rownames(mammal_data)<-mammal_toren_anage_pacif_df_c2$Species

## run PGLS
data_mammal_comparative<-comparative.data(phy=mammal_tree_sp600m_clustal,data=mammal_data,names.col=Species,vcv=TRUE,force.root=TRUE, warn.dropped=T)

pgls_mammal<-pgls(MLS_Log ~ trip+GC+GC_skew+AT_skew+BM_Log, data=data_mammal_comparative,lambda='ML')
summary(pgls_mammal)