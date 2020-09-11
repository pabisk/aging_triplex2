##################################################################################################################
############################################## find all direct repeats

## example dataset
s<-readDNAStringSet(".../Lakshamanan 294 mammals.fasta",format="fasta")


#
library("Biostrings")
library("stringr")

# get repeats for a set of species
get_all_repeats_of_DNAStringset<-function(mtGenomes=s,min=10,kmer_length=13,unique_repeats=TRUE,repeats_type="direct",
                                          remove_overlapping=TRUE, ralgo="perfect",
                                          return_repeat_coord=FALSE, heuristic_exclude_DLoop=FALSE,return_repeat_seq=FALSE){
  library(stringr)
  library(stringi)
  print("calculating repeats for species: ")
  rcord<-list()

  
  if(repeats_type=="all"){ ## this will then return DR, MR, IR and ER in that order, species=row, columns=repeat types
    ## this will totally ignore the min=x value
    m <- matrix(0, ncol = 4, nrow = 0)## ncol is equal to 4 because there are 4 repeat types
    reps<-data.frame(m)
    
    for(a in 1:length(mtGenomes)){ 
      print(names(mtGenomes)[a])
      repeats<-get_all_repeat_types(mtGenomes[a],kmer_length=kmer_length,unique_repeats=unique_repeats,
                                    repeats_type=repeats_type,heuristic_exclude_DLoop=heuristic_exclude_DLoop)
      repeats<-t(data.frame(repeats)) #need to transpose because per default the created df has only 1 column & 4 rows
      
      reps<-rbind(reps,repeats)
      print(a)
    }
    names(reps)<-c("direct","mirror","inverted","everted")
  } else {
    ## this is the standard function that allows you to iterate between a min and max repeat length
    if(ralgo=="perfect"){
      m <- matrix(0, ncol = kmer_length-min+1, nrow = 0)
      reps<-data.frame(m)
      
      for(a in 1:length(mtGenomes)){
        print(names(mtGenomes)[a])
        repeats<-get_all_repeats_upto_kmer(mtGenomes[a],min=min,kmer_length=kmer_length,unique_repeats=unique_repeats,
                                           repeats_type=repeats_type,remove_overlapping=remove_overlapping,
                                           heuristic_exclude_DLoop=heuristic_exclude_DLoop,return_repeat_seq=return_repeat_seq)
        
        if(return_repeat_seq==TRUE){
          rcord[[a]]<-repeats
        }else{
          repeats<-t(data.frame(repeats)) # need to transpose because per default the created df has only 1 column & 4 rows
          
          reps<-rbind(reps,repeats) 
        }
        
        print(a)
      }
    }
    
  }
  if(return_repeat_seq==TRUE){return(rcord);            }
  return(reps)
}

# get repeats for a single species
get_all_repeats_upto_kmer<-function(mtGenome=s[250],min=10,kmer_length=13,unique_repeats=TRUE,return_repeat_seq=F,
                                    repeats_type="direct",ralgo="perfect",remove_overlapping=TRUE,
                                    return_repeat_coord=FALSE,heuristic_exclude_DLoop=FALSE){
  print("calculating repeats for kmer length of =")
  print(kmer_length)
  counter=1
  rep_seq<-vector()
  
  if(ralgo=="perfect"){
    repeat_results<-rep(1:((kmer_length-min)+1))
    for(xy in min:kmer_length){
      message("xy")
      if(return_repeat_seq==FALSE){
        repeat_results[counter]<-repeats_of_given_kmer(mtGenome,kmer_length=xy,return_repeat_seq=return_repeat_seq,  unique_repeats=unique_repeats,
                                                       repeats_type=repeats_type,remove_overlapping=remove_overlapping,
                                                       return_repeat_coord=return_repeat_coord,heuristic_exclude_DLoop=heuristic_exclude_DLoop)
      }
      counter<-counter+1
      if(return_repeat_seq==TRUE){
        res<-repeats_of_given_kmer(mtGenome,kmer_length=xy,return_repeat_seq=return_repeat_seq,  unique_repeats=unique_repeats,
                                   repeats_type=repeats_type,remove_overlapping=remove_overlapping,
                                   return_repeat_coord=return_repeat_coord,heuristic_exclude_DLoop=heuristic_exclude_DLoop)
        rep_seq<-c(rep_seq,res)
      }
      cat(xy)
    }
  }
  
  ## returns a vector for all repeats that have length between min and kmer_length
  if(return_repeat_seq==TRUE){ return(rep_seq) }
  return(repeat_results)
}

###############
# get repeats of fixed length for a single species
repeats_of_given_kmer<-function(mtGenome=s[250],kmer_length=13, clean=TRUE, return_repeat_seq=FALSE,
                                unique_repeats=TRUE,remove_overlapping=TRUE,return_repeat_coord=FALSE,
                                heuristic_exclude_DLoop=FALSE,repeats_type="direct"){
  mtg<-DNAString(paste(mtGenome))
  if(heuristic_exclude_DLoop==TRUE){
    message("excluding D loop by heuristic approach")
    mtg<-mtg[1000:15400]
  }
  
  if(clean==TRUE){
    mtg<-symbolClean(DNAStringSet(mtg))
    mtg<-paste(mtg)
  }else{  mtg<-paste(mtg) }
  
  if(repeats_type=="inverted"){ mtg_inverted<-reverseComplement(mtGenome) }
  if(repeats_type=="mirror"){ mtg_inverted<-reverse(mtGenome) }
  if(repeats_type=="everted"){ mtg_inverted<-Biostrings::complement(mtGenome) }
  
  ## generate truncated sequence; this should be fast
  ## for kmer=4 we will generate 4 mtDNAs truncated by one nt each
  ## this will alow us to split the mtDNA into chunks as long as the potential repeat
  truncated_strings<-rep(1:(kmer_length+1))
  for(i in 0:kmer_length){
    truncated_strings[i]<-substring(mtg,i,nchar(mtg)) ##start=i, end=nchar=mtDNA length=mtDNA end
  }

  ### get all the kmers for each of the truncated mtDNAs which should cover all unique repeats
  ### split the mtDNA every k bases, based on the kmer length
  substrings<-rep(1:kmer_length)
  for(n in 1:length(truncated_strings)){
    substrings[n]<-fixed_split(truncated_strings[n],kmer_length)##returns a list
  }
  
  ## exclude all duplicates and strings < kmer length
  ## in this step we are identifying all *potential* repeats
  substrings_unique<-unique(unlist(substrings))
  substrings_unique<-substrings_unique[nchar(substrings_unique)>=kmer_length]
  cat("Identified X unique substrings that could be potential repeats:")
  cat(length(substrings_unique))
  
  coords_of_repeatsv2<-vector()
  
  ## now count the repeats that actually match the mtDNA
  nr_repeats<-rep(1:length(substrings_unique))
  
  if(repeats_type=="direct"){
    print("checking direct repeats")
    for(x in 1:length(substrings_unique)){
      count_repeats<-str_count(mtg,substrings_unique[x]) ## can you find substrings_unique in mtg genome
      
      ## get the coordinate of the repeats (substrings_unique[x]) WITHIN the mtDNA (mtg)
      if(count_repeats>=2){
        interm<-str_locate_all(mtg, substrings_unique[x])[[1]][,1]
        coords_of_repeatsv2<-c(coords_of_repeatsv2,interm)
      }
    }
  }
  print("length(substrings_unique)")
  print(length(substrings_unique))
  
  if(repeats_type=="inverted" || repeats_type=="everted" || repeats_type=="mirror"){
    print("checking non direct repeats")
    for(x in 1:length(substrings_unique)){
      count_repeats<-str_count(mtg_inverted,substrings_unique[x]) ## can you find substrings_unique in mtg genome
      
      ## get the coordinate of the repeats (substrings_unique[x]) WITHIN the mtDNA (mtg)
      if(count_repeats>=1){
        interm<-str_locate_all(mtg, substrings_unique[x])[[1]][,1]
        coords_of_repeatsv2<-c(coords_of_repeatsv2,interm)
      }
    }
  }
  
  ## now handle zero repeats
  if(length(coords_of_repeatsv2)==0){
    print("no repeats found")
    if(return_repeat_seq==TRUE) return("")
    return(0)
  }
  
  ### the below code will exclude OVERLAPPING (but not redundant repeats thus double counting some)
  sr<-sort(coords_of_repeatsv2)
  if(remove_overlapping==TRUE){
    counter=1
    while(counter<=length(sr)){
      ## get coordinates for all repeats whose start point is equal or below the current one
      ## or whose start point does NOT overlap with startpoint of current one + repeat length
      sr<-sr[sr<=sr[counter] | sr>=(sr[counter]+kmer_length)]
      counter=counter+1
    }
  }
  
  ## get the sequences based on the positions
  return_vector<-DNAStringSet(rep("A",length(sr)))
  mtg<-DNAString(mtg)
  for(i in 1:length(sr)){
    # recover the sequences
    return_vector[[i]]<-mtg[sr[i]:(sr[i]+(kmer_length-1))]
    
  }
  ## now we also exclude redundant repeats
  if(unique_repeats==TRUE){
    return_vector<-levels(as.factor(return_vector))
  }
  if(return_repeat_seq==TRUE){return(return_vector)}
  
  
  return(length(return_vector))##returns the (cleaned) start coordinates for all repeats
}
