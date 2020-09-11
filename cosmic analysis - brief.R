################################################################################################################
############## COSMIC analysis
## all GRCh==38
cosmic<-read.delim(".../CosmicBreakpointsExport.tsv.txt") ## truncated file uploaded to github
cosmic2<-data.frame(c(cosmic$Chrom.From,cosmic$Chrom.To),
                    c(cosmic$Location.From.min,cosmic$Location.To.min))
names(cosmic2)<-c("chr","bp")

## exclude sex chromsomes and known problematic breakpoints
cosmic2_subset<-cosmic2
cosmic2_subset<-cosmic2_subset[-which(cosmic2_subset$chr==23),]
cosmic2_subset<-cosmic2_subset[-which(cosmic2_subset$chr==24),]
cosmic2_subset<-cosmic2_subset[-which(cosmic2_subset$bp==59118984),]

## analysis
cosmic_motif_test_full(input=cosmic2_subset,blocksize=2000, control=3000,around=250,sleeptimer=0.1,
                       type="triplex", error_handling=T,name_tag="full analysis")
cosmic_motif_test_full(input=cosmic2_subset,blocksize=2000, control=3000,around=250,sleeptimer=0.1,
                       type="GQ", error_handling=T,name_tag="full analysis")

################################################################################################################
############## COSMIC code
## requires code from: triplex and GQ detection.R

cosmic_motif_test_full<-function(input=cosmic2_subset,blocksize=500, control=3000,around=500,sleeptimer=0.1, type="triplex", error_handling=T,name_tag=""){
  ## measure number of GQ or Triplex motifs around COSMIC breakpoints in the human genome
  
  ## type = analyze triplex or GQ
  ## error_handling = if T ignore and discard the whole block of sequences if one of them cannot be found by get_sequence()
  ## around = consider X bp around the breakpoint
  ## control = shift the BP by +X bp which serves as a control
  ## blocksize = use this blocksize to get get_sequence() from human genome; that way we only need to hold one block of sequences in memory
  ## input = must be a df with chromosome and BP coordinates; incomaptible with X and Y chromsome numbering as chr23 and chr24
  library("BSgenome.Hsapiens.NCBI.GRCh38")
  xlen<-(floor(dim(input)[1]/blocksize))
  rez<-vector();  ctrl<-vector()
  rand_name<-sample(1:100000000,1)
  
  for(n in 1:xlen){
    bps<-cosmic2_subset[(((n-1)*blocksize)+1):(n*blocksize),]
    #print((((n-1)*blocksize)+1):(n*blocksize))
    
    if(error_handling==T){
      ## in rare instaces the coordinates are wrong, one mistake around every 50 000 BPs, probably due to +/- adding of coords
      seqs = tryCatch({
        get_sequence(chromosome=bps$chr,breakpoint=bps$bp,around=around, control=0)
      }, error = function(error_condition) {
        seqs<-0
      })
      seqs_ctrl = tryCatch({
        get_sequence(chromosome=bps$chr,breakpoint=bps$bp,around=around,control=control)
      }, error = function(error_condition) {
        seqs_ctrl<-0
      })
    }else {
      seqs<-get_sequence(chromosome=bps$chr,breakpoint=bps$bp,around=around, control=0)
      seqs_ctrl<-get_sequence(chromosome=bps$chr,breakpoint=bps$bp,around=around,control=control)
    }
    
    if(type=="triplex"){
      result_file<-paste0(name_tag," rez-data-trip",rand_name,".txt")
      ctrl_file<-paste0(name_tag," ctrl-data-trip",rand_name,".txt")
      rez<-triplex_search(stringset=seqs,stringset_sampled=seqs,return_seq=F,min_scor=15)
      ctrl<-triplex_search(stringset=seqs_ctrl,stringset_sampled=seqs_ctrl,return_seq=F,min_scor=15)
      
    }
    
    if(type=="GQ"){
      result_file<-paste0(name_tag," rez-data-GQ",rand_name,".txt")
      ctrl_file<-paste0(name_tag," ctrl-data-GQ",rand_name,".txt")
      rez<-pqs_quads(mtDNAstringsSet=seqs,mtDNAseq_sampled=seqs)$count
      ctrl<-pqs_quads(mtDNAstringsSet=seqs_ctrl,mtDNAseq_sampled=seqs_ctrl)$count
    }
    ctrl<-cbind(ctrl,bps$chr,bps$bp)
    rez<-cbind(rez,bps$chr,bps$bp)
    write.table(rez,file=result_file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
    write.table(ctrl,file=ctrl_file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
    
    Sys.sleep(sleeptimer)
  }
  return(list(rez,ctrl))
}
