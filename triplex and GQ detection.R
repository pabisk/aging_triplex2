################################################################################################################
############## GQ code
## count GQ motifs with pqsfinder

pqs_quads <- function(mtDNAstringsSet=s,mtDNAseq_sampled=seq_sampled, stringency=26, quietly=FALSE, return_starts=FALSE, return_midpoint=FALSE,
                      return_width=FALSE, return_mean_score=FALSE, return_seqs=FALSE) {
  ## stringency=26; default for pqsfinder
  ## return_seqs=T; gives you GQ sequences
  
  library("pqsfinder")
  # assumes input is +strand and counts GQ on both strands but this is not relevant
  # as GQ are not abundant on the other strand
  
  if(class(mtDNAstringsSet)!="DNAStringSet"){
    return(list(count=-2, sampled=-2))
  }
  
  if(return_midpoint ==TRUE){
    ## in this mode we do not return sampled dataset
    mtDNAseq_sampled<-mtDNAstringsSet
  }
  
  results_pqsfinder<-rep(0,length(mtDNAstringsSet))
  starts_pqsfinder<-list()
  mids_pqsfinder<-list()
  seqs_pqsfinder<-list()
  score_pqsfinder<-list()
  results_pqsfinder_sampled<-rep(0,length(mtDNAstringsSet))
  
  results_pqsfinder_cds<-rep(0,length(mtDNAstringsSet))
  cat("stringency is = ", stringency, "\n")
  for (i in c(1:length(mtDNAstringsSet))){
    cat("seaching spcies n :", i)
    cat("\n")
    string<-DNAString(paste(mtDNAstringsSet[i]))
    if(quietly==TRUE){
      result<-suppressMessages(pqsfinder(string, min_score=stringency))
    }else {
      result<-pqsfinder(string, min_score=stringency)
    }
    
    results_pqsfinder[i]<-length(result)
    starts_pqsfinder[[i]]<-start(result)
    mids_pqsfinder[[i]]<-start(result)+round((width(result)/2))
    seqs_pqsfinder[[i]]<-paste(result)
    score_pqsfinder[[i]]<-score(result)
    
    ## repeat for GC matched random sequence
    string<-DNAString(paste(mtDNAseq_sampled[i]))
    if(quietly==TRUE){
      results_pqsfinder_sampled[i]<-suppressMessages(length(pqsfinder(string, min_score=stringency)))
    }else{
      results_pqsfinder_sampled[i]<-length(pqsfinder(string, min_score=stringency))
    }
    
  }
  ret<-data.frame(results_pqsfinder,results_pqsfinder_sampled)
  ## can be plotted
  plot(ret$results_pqsfinder, ret$results_pqsfinder_sampled, xlim=c(0,120), ylim=c(0,120))
  abline(0,1)
  
  if(return_mean_score==TRUE){
    return(score_pqsfinder)
  }
  
  if(return_midpoint==TRUE){
    return(mids_pqsfinder)
  }
  
  if(return_seqs==TRUE){
    return(seqs_pqsfinder)
  }
  
  if(return_starts==TRUE){
    return(starts_pqsfinder)}
  
  return(list(count=results_pqsfinder, sampled=results_pqsfinder_sampled))
}



################################################################################################################
############## triplex code
## count triplex motfis with the triplex package

## quantify triplex motifs in dna stringset
triplex_search<-function(stringset=s,stringset_sampled=seq_sampled,pval=1,min_scor=15,return_both=F,
                         return_seq=T,return_full=F, type=0:7,strand="any", fast_mode=T){
  ## return_seq=T returns triplex sequences
  ## type=0:7 tests all known types of intramolecular triplex motifs
  ## strand=any default setting considers both + and - strabd
  ## pval=1 p value setting ignored
  ## min_scor=15 default; use this minimum score
  ## return_both=T; you may want to compare two stringsets at the same time
  library("triplex") ## from the bioconductor package
  if(class(stringset)!="DNAStringSet"){
    return(-2)
  }
  
  tripl_res<-vector()
  tripl_res_sampled<-vector()
  seqs<-list()
  full_data<-list()
  for(i in 1:length(stringset)){
    res<-triplex.search(stringset[[i]], min_score=min_scor, p_value=pval,type=type)
    
    if(strand=="+" || strand=="-" || strand=="plus"|| strand=="minus"){
      plus_strand<-strand(res)=="+"
      minus_strand<-strand(res)=="-"
      if(strand=="+" || strand=="plus") { res<-res[plus_strand] }
      if(strand=="-"|| strand=="minus") { res<-res[minus_strand] }
    }
    
    
    full_data[[i]]<-res
    tripl_res[i]<-length(res)
    if(fast_mode==F) {
      tripl_res_sampled[i]<-length(triplex.search(stringset_sampled[[i]], min_score=min_scor, p_value=pval,type=type))
    } else{
      tripl_res_sampled[i]<-0
    }
    if(return_seq==T){
      seqs[[i]]<-paste(res)
    }
  }
  if(return_full==T){ return(full_data)}
  if(return_seq==T){ return(seqs)}
  if(return_both==T){ return(list(tripl_res,tripl_res_sampled))}
  return(tripl_res)
}

## we can check different minimum score settings of triplex_search()
triplex_search_space<-function(dna_string_set=s,min=8,max=22,correl_df=log(df$Lifespan),return_single_sp_data=F,return_seq=F,strand="any"){
  singl_sp_df<-vector()
  sample_size<-vector()
  tripl_cor<-vector()
  for(i in min:max){
    message(i)
    
    if(return_single_sp_data==F){
      temp_tripl<-triplex_search(stringset=dna_string_set,stringset_sampled=dna_string_set,min_scor=i,return_seq=F,strand=strand) # 
      sample_size[i]<-mean(temp_tripl)
      print(temp_tripl)
      print(cor(temp_tripl,correl_df))
      
      tripl_cor[i]<-cor(temp_tripl,correl_df)
      if(is.na(cor(temp_tripl,correl_df)==TRUE)){tripl_cor[i]<-2}
      print(tripl_cor[i])
    }
    
    ## we can also just give one species as dnastringset and then search parameter space for this species
    if(return_single_sp_data==T){singl_sp_df[i]<-length(triplex.search(dna_string_set, min_score=i, p_value = 1))}
  }
  tripl_cor<-tripl_cor[!is.na(tripl_cor)]
  singl_sp_df<-singl_sp_df[!is.na(singl_sp_df)]
  sample_size<-sample_size[!is.na(sample_size)]
  
  if(return_single_sp_data==T){ return(singl_sp_df) }
  
  return(rbind(tripl_cor,sample_size))
}

