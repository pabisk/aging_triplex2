####################################################################################################################################################
#### breakpoint overlap functions

## example dataset from Lakshamanan et al. 2015
s<-readDNAStringSet(".../Lakshamanan 294 mammals.fasta",format="fasta")
hrep<-get_all_repeats_upto_kmer(min=6,kmer_length=15,unique_repeats = F,return_repeat_seq = T, remove_overlapping = F)
hrep<-sapply(hrep,as.character); hrep<-unlist(hrep) ## get all direct repeats of 6 to 15bp length

#
library("Biostrings")
library("stringr")

# generate controls
rep_bp_overlap2_multirun<-function(times=20,mtdna_seq=s[[250]],breakpoint_df=bps,kmer_min=6,kmer_length=15,window_delta=0,leftshift=0,detailed=F, type="direct",
                                   type2="single intersect", min = 5747, max = 16500,quietly=T, use_reb_bp1_method=F,repeat_df=hrep,
                                   rm_overlapping=T,sleepy=1,return_special_prime=0,mode="default"){
  multirun_df<-vector()
  print(head(breakpoint_df))
  Sys.sleep(sleepy)
  for(n in 1:times){
    Sys.sleep(sleepy/2)
    message(n)
    ctrl<-generate_reshuffled_bps_distribution_matched(breakpoints=breakpoint_df,min = min, max = max, mode=mode)
    head(ctrl)
    if(use_reb_bp1_method==T){
      ## the aeguments rm_overlapping and repeat_df are needed if we want to use rep_bp_overlap() method
      print("running old method")
      res<-rep_bp_overlap(mtdna_seq=mtdna_seq, breakpoint_df=ctrl, repeat_df=repeat_df, window_delta=window_delta,quietly=quietly,
                          leftshift=leftshift, detailed=detailed, type=type, type2=type2, rm_overlapping=rm_overlapping,return_special_prime=return_special_prime) 
    }else{
      if(quietly==T){
        res<-suppressMessages(rep_bp_overlap2(mtdna_seq=mtdna_seq, breakpoint_df=ctrl, kmer_min=kmer_min, kmer_length=kmer_length, window_delta=window_delta,
                                              leftshift=leftshift, detailed=detailed, type=type, type2=type2,return_special_prime=return_special_prime))
      }else{
        res<-rep_bp_overlap2(mtdna_seq=mtdna_seq, breakpoint_df=ctrl, kmer_min=kmer_min, kmer_length=kmer_length, window_delta=window_delta,
                             leftshift=leftshift, detailed=detailed, type=type, type2=type2,return_special_prime=return_special_prime) 
      }
    }
    
    
    multirun_df[n]<-length(res)
  }
  
  return(multirun_df)
}


# find all features near breakpoints
rep_bp_overlap<-function(mtdna_seq=s[[250]],breakpoint_df=bps,repeat_df=hrep,intersect=TRUE,type="direct",type2="single intersect",min_score=16,detailed=F, window_delta=0,
                         leftshift=0, rm_overlapping=F, quietly=T, return_special_prime=0){
  ## repeat_df has to be a df, preferably of UNIQUE repeat sequences
  ## check if ANY repeat from repeat_df is in the left breakpoint and return ALL such repeats
  ## check the same for the right breakpoint, then you can intersect the two and find repeats in both
  
  ## if you want to look at 5p or 3p associated features specifically
  ## arguments need to be detailed==T + return_special_prime=3 OR return_special_prime=5 then we can also specifically return 5' or 3' associated features in DIP mode
  
  library("stringr")
  library("stringi")
  
  ##########
  ### This allows you to adjust the size of the search window or shift it
  cat("ATTENTION window size (excluding central base) is.... ")
  wsz<-(-(breakpoint_df$n5pA-breakpoint_df$n5pB)[1])
  if(quietly==F) { cat(wsz) }
  print("")
  ## leftshift means a shift towards lower numbers, which is AWAY from OH, counter clockwise
  
  breakpoint_df$n5pA<-(breakpoint_df$n5pA-leftshift)
  breakpoint_df$n3pA<-(breakpoint_df$n3pA-leftshift)
  breakpoint_df$n5pB<-(breakpoint_df$n5pB-leftshift)
  breakpoint_df$n3pB<-(breakpoint_df$n3pB-leftshift)
  
  if(window_delta>0){
    print("changing the window size to fit your predefined windowsize")
    ## actual window size is larger than the desired window size -> decrease wsz
    
    xd<-(wsz-window_delta)/2
    breakpoint_df$n5pA<-breakpoint_df$n5pA+xd
    breakpoint_df$n3pA<-breakpoint_df$n3pA+xd
    breakpoint_df$n5pB<-breakpoint_df$n5pB-xd
    breakpoint_df$n3pB<-breakpoint_df$n3pB-xd
    
    cat("ATTENTION the NEW window size (excluding central base) is.... ")
    wsz<-(-(breakpoint_df$n5pA-breakpoint_df$n5pB)[1])
    cat(wsz)
    print("")
  }
  
  ##########
  
  names(breakpoint_df)<-c("n5pA","n5pB","n3pA","n3pB")
  
  repeat_count<-0
  repeat_pos_in_df<-vector()
  detailed_res_df<-data.frame(stringsAsFactors = F)
  detailed_res_df2<-data.frame(stringsAsFactors = F) ## store data when you do a double independent
  n_intersection2<-0
  
  for(i in 1:dim(breakpoint_df)[1]){

    ##########
    n5p_seq<-paste0(mtdna_seq[breakpoint_df$n5pA[i]:breakpoint_df$n5pB[i]])
    n3p_seq<-paste0(mtdna_seq[breakpoint_df$n3pA[i]:breakpoint_df$n3pB[i]])
    
    if(type=="direct"){
      ## given a vector with thousands of repeats, repeat_df, extract the repeats that are found within the breakpoints
      n5_extract<-str_extract(n5p_seq,repeat_df)
      n3_extract<-str_extract(n3p_seq,repeat_df)
    }
    if(type=="mirror" || type=="inverted" || type=="everted"){
      ## here we have to consider that the repeats are not located in the same direction and/or on the same strand
      
      if(type=="inverted"){
        n3_extract<-str_extract(paste0(reverseComplement(DNAStringSet(n3p_seq))),repeat_df)
        n5_extract2<-str_extract(paste0(reverseComplement(DNAStringSet(n5p_seq))),repeat_df)
      }
      if(type=="mirror"){
        n3_extract<-str_extract(reverse(n3p_seq),repeat_df)
        n5_extract2<-str_extract(reverse(n5p_seq),repeat_df)
      }
      if(type=="everted"){
        n3_extract<-str_extract(paste0(Biostrings::complement(DNAStringSet(n3p_seq))),repeat_df)
        n5_extract2<-str_extract(paste0(Biostrings::complement(DNAStringSet(n5p_seq))),repeat_df)
      }
      
      n5_extract<-str_extract(n5p_seq,repeat_df)
      n3_extract2<-str_extract(n3p_seq,repeat_df)
      
      n5_extract2<-n5_extract2[!is.na(n5_extract2)]
      n3_extract2<-n3_extract2[!is.na(n3_extract2)]
      
    }
    n5_extract<-n5_extract[!is.na(n5_extract)]
    n3_extract<-n3_extract[!is.na(n3_extract)]
    
    ################# default part
    intersection_length<-0
    n_intersection<-0
    if(intersect==TRUE){
      ## default=TRUE; requires the SAME feature to be found in both BPs
      ## if we intersect with GQ sequences we may prefer not to enforce the existence of the seq at BOTH breakpoints
      if(type=="direct"){
        n_intersection<-intersect(n5_extract,n3_extract)
        intersection_length<-length(n_intersection)
      }
      
      n_intersection2<-0
      if(type=="mirror" || type=="inverted" || type=="everted"){
        ## intersection has to be cross-wise between the given and the complemented/reversed/revcomp breakpoint region
        
        n_intersection<-intersect(n5_extract,n3_extract)
        intersection_length<-length(n_intersection)
        n_intersection2<-intersect(n5_extract2,n3_extract2)
        
        
        intersection_length2<-length(n_intersection2)
        intersection_length<-max(intersection_length,intersection_length2)
      }
      
    }else {
      ## this the case where we count ANY overlap with ANY of the breakpoint sequences
      ## useful for GQs
      intersection_length<-max(length(n5_extract),length(n3_extract))
    }
    #################
    
    ################# DIP part (for e.g. GQs, IRs, etc)
    if(type2=="double independent"){
      ## count only if both repeats have an element, but the element doesnt have to be the same
      ## below is for MR, IR and ER
      if(type=="mirror" || type=="inverted" || type=="everted"){
        dip_intersection1<-intersect(n5_extract,n5_extract2)
        dip_intersection2<-intersect(n3_extract,n3_extract2)
        dip_intersection_length<-max(length(dip_intersection1),length(dip_intersection2))
        
        ## the below code can be completely ignored for the standard GQ, Triplex cases as well as elements "spanning" both breakpoints
        if(rm_overlapping==T){
          ## Notes while for the search input, you want to include all potential repeats, also overlapping ones, during post-processing you have
          ## to remove them. OTherwise you will catch e.g.  inverted repeats that overlap with themselves.
          ## (You include overlapping XR initially because partly overlapping XR1 and XR2 may be located near the border of a BP, you could fail to detect these if you exclude the wrong one)
          if(length(dip_intersection1)>0){
            for(ab in 1:length(dip_intersection1)){
              locate1<-str_locate(n5p_seq,dip_intersection1[ab])
              if(type=="inverted"){ locate2<-str_locate(n5p_seq,paste(reverseComplement(DNAString(dip_intersection1[ab])))) }
              if(type=="everted"){ locate2<-str_locate(n5p_seq,paste(complement(DNAString(dip_intersection1[ab])))) }
              if(type=="mirror"){ locate2<-str_locate(n5p_seq,paste(reverse(DNAString(dip_intersection1[ab])))) }
              
              if((locate2[1,1]<=locate1[1,2]) &(locate2[1,1]>=locate1[1,1])){
                message("overlapping repeat, removing")
                dip_intersection1[ab]<-""
              }
              if((locate2[1,2]<=locate1[1,2]) &(locate2[1,2]>=locate1[1,1])){
                message("overlapping repeat, removing")
                dip_intersection1[ab]<-""
              }
            }
          }
          if(length(dip_intersection2)>0){
            for(ab in 1:length(dip_intersection2)){
              locate1<-str_locate(n3p_seq,dip_intersection2[ab])
              if(type=="inverted"){ locate2<-str_locate(n3p_seq,paste(reverseComplement(DNAString(dip_intersection2[ab])))) }
              if(type=="mirror"){ locate2<-str_locate(n3p_seq,paste(reverse(DNAString(dip_intersection2[ab])))) }
              if(type=="everted"){ locate2<-str_locate(n3p_seq,paste(Biostrings::complement(DNAString(dip_intersection2[ab])))) }
              
              if((locate2[1,1]<=locate1[1,2]) &(locate2[1,1]>=locate1[1,1])){
                message("overlapping repeat, removing")
                dip_intersection2[ab]<-""
              }
              if((locate2[1,2]<=locate1[1,2]) &(locate2[1,2]>=locate1[1,1])){
                message("overlapping repeat, removing")
                dip_intersection2[ab]<-""
              }
            }
          }
        }
        
        det_temp_df2<-data.frame(i,paste0(n5p_seq),paste0(n3p_seq),
                                 breakpoint_df$n5pA[i],breakpoint_df$n5pB[i],
                                 breakpoint_df$n3pA[i],breakpoint_df$n3pB[i],intersection_length,
                                 length(dip_intersection1),max(dip_intersection1),
                                 length(dip_intersection2), max(dip_intersection2),
                                 stringsAsFactors = F)
        
        detailed_res_df2<-rbind(detailed_res_df2,det_temp_df2)
        
      }else{
        ## this is for GQs, triplex etc
        n_intersection<-""
        n_intersection2<-""
        if(length(n5_extract)>0 || length(n3_extract)>0){
          repeat_count<-repeat_count+1
          repeat_pos_in_df<-c(repeat_pos_in_df,i)
          if(length(n5_extract)==0){n_intersection<-""
          } else {
            n_intersection<-n5_extract[which(nchar(n5_extract)==max(nchar(n5_extract)))[1]]
          }
          if(length(n3_extract)==0){
            n_intersection2 <-""
          } else{
            n_intersection2<-n3_extract[which(nchar(n3_extract)==max(nchar(n3_extract)))[1]]
          }
        }
        ## prepare a df for the detailed output
        det_temp_df2<-data.frame(i,paste0(n5p_seq),paste0(n3p_seq),
                                 breakpoint_df$n5pA[i],breakpoint_df$n5pB[i],
                                 breakpoint_df$n3pA[i],breakpoint_df$n3pB[i],intersection_length,
                                 nchar(max(n_intersection)),max(n_intersection),
                                 nchar(max(n_intersection2)),max(n_intersection2), 
                                 stringsAsFactors = F)
        detailed_res_df2<-rbind(detailed_res_df2,det_temp_df2)
      }
      
    }else { ## "single intersect" == default
      if(intersection_length>0){
        repeat_count<-repeat_count+1
        repeat_pos_in_df<-c(repeat_pos_in_df,i)
      }
      #################
      
      ## populate the data.frame with exact data in case you want detailed output
      det_temp_df<-data.frame(i,paste0(n5p_seq),paste0(n3p_seq),
                              breakpoint_df$n5pA[i],breakpoint_df$n5pB[i],
                              breakpoint_df$n3pA[i],breakpoint_df$n3pB[i],intersection_length,max(n_intersection),max(n_intersection2),nchar(max(n_intersection)),
                              nchar(max(n_intersection2)), stringsAsFactors = F)
      detailed_res_df<-rbind(detailed_res_df,det_temp_df)
    }
    
  }## end for loop
  
  ## clean up code, format output, etc.
  if(type2=="double independent"){
    names(detailed_res_df2)[2]<-"n5p_bp_region";names(detailed_res_df2)[3]<-"n3p_bp_region"
    names(detailed_res_df2)[4:7]<-c("n5pA","n5pB","n3pA","n3pB");names(detailed_res_df2)[10]<-c("dip_result1");names(detailed_res_df2)[12]<-c("dip_result2")
    if(detailed==F){
      
      detailed_res_df2$dip_result1[is.na(detailed_res_df2$dip_result1)==TRUE]=""; detailed_res_df2$dip_result2[is.na(detailed_res_df2$dip_result2)==TRUE]=""
      rez<-nchar(detailed_res_df2$dip_result1)+nchar(detailed_res_df2$dip_result2); if(quietly==F) { print(rez); }
      return(which(rez>0)) ##simplified return vector
      
    }
    if(detailed==T){ ## arguments need to be detailed==T + return_special_prime=3 OR return_special_prime=5 then we can also specifically return 5' or 3' associated features in DIP mode
      if(return_special_prime>0){
        ## we can also specifically return 5' or 3' associated features in DIP mode
        if(return_special_prime==5){ return(which(nchar(detailed_res_df2$dip_result1)>0)) }
        if(return_special_prime==3){ return(which(nchar(detailed_res_df2$dip_result2)>0)) }
      }
      return(detailed_res_df2)
    }
    
  }
  
  if(detailed==TRUE){
    names(detailed_res_df)[2]<-"n5p_bp_region";names(detailed_res_df)[3]<-"n3p_bp_region"
    return(detailed_res_df)
  }
  return(repeat_pos_in_df)
}


# helper function
generate_reshuffled_bps_distribution_matched<-function(min=5747,max=16500,window=50, breakpoints=breakpoints_human_mitobreak, quietly=T, mode="default"){
  ## this allows two shuffling modes that produce similar results
  ## mode default: keep the deletion size fixed when reshuffling
  ## mode new: use the breakpoint distribution to generate reshuffled breakpoints
  if(dim(breakpoints)[2]==4){
    message("transform breakpoint REGION to break POINT")
    ## in thise case we are dealing with breakpoint WINDOWS
    half_delta<-(breakpoints$n5pB-breakpoints$n5pA)/2
    n5p<-breakpoints$n5pA+((breakpoints$n5pB-breakpoints$n5pA)/2)
    n3p<-breakpoints$n3pA+((breakpoints$n3pB-breakpoints$n3pA)/2)
    if(quietly==F){print("dim = 4")}
  }else{
    n5p<-breakpoints$n5p_break
    n3p<-breakpoints$n3p_break
    if(quietly==F){print("dim = 2")}
  }
  bp_dist<-n3p-n5p
  
  if(quietly==F){print(head(bp_dist))}

  nu_reshuffled_bps<-data.frame()
  
  lengthn<-dim(breakpoints)[1]
  if(quietly==F){print(lengthn)}
  
  if(mode=="default"){
    
    for(i in 1:lengthn){
      a<-sample(min:max,1,replace=T)
      b<-a+bp_dist[i]
      if(quietly==F){message(i)}
      while(b>max){
        a<-sample(min:max,1,replace=T)
        b<-a+bp_dist[i]
      }
      nu_reshuffled_bps<-rbind(nu_reshuffled_bps,data.frame(a,b))
    }
    
  } else {
    ## mode="new"; generate matched BPs based on the distribution
    x<-density((bps_mitobreak$n5pA+bps_mitobreak$n5pB)/2,n=length(bps_mitobreak$n5pA))
    samp_n5<-sample(x$x, lengthn, replace=T, prob=x$y)
    
    y<-density((bps_mitobreak$n3pA+bps_mitobreak$n3pB)/2,n=length(bps_mitobreak$n3pA))
    samp_n3<-sample(x$x, lengthn,replace=T, prob=y$y)
    
    nu_reshuffled_bps<-cbind(round(samp_n5),round(samp_n3))
    nu_reshuffled_bps<-data.frame(nu_reshuffled_bps)
    names(nu_reshuffled_bps)<-c("a","b")
    print(nu_reshuffled_bps)
  }
  
  
  if(dim(breakpoints)[2]==4){
    ## return window if window was input (not POINT)
    nu_reshuffled_bps<-data.frame(nu_reshuffled_bps$a-half_delta, nu_reshuffled_bps$a+half_delta,nu_reshuffled_bps$b-half_delta, nu_reshuffled_bps$b+half_delta)
    names(nu_reshuffled_bps)<-c("n5pA","n5pB","n3pA","n3pB")
  }
  
  return(nu_reshuffled_bps)
}

