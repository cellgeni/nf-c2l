#!/usr/bin/env Rscript
# somehow in Rscript libPath order is wrong, so I need to overwrite it
.libPaths("/nfs/cellgeni/R/x86_64-pc-linux-gnu-library/4.3")
library(spacexr)
library(Seurat)
library(argparse)

print(.libPaths())

get_parser = function(){
  parser <- ArgumentParser()
  
  parser$add_argument("--ref_h5ad", type = "character",
                      help = "Path to reference h5ad")
  
  parser$add_argument("--vis_h5ads", type = "character",nargs='+',
                      help = "Path to visium h5ad[s]. Will concatenate with spatial offetting (to avoid overlaps) if mutliple samples are provided. Though it seems that RCTD is not using coordinates for deconvolution, so offsetting is not really essential.")
  
  parser$add_argument("--outbase", type = "character",
                      help = "base name for output files")
  
  parser$add_argument("--doublet_mode", type = "character",
                      help = "RCTD doublet_mode: full or multi (can be also doublet)")
  
  parser$add_argument("--celltype_col", type = "character",
                      help = "name of celltype filed in adata.obs. Precomputed c2l reference will be used if the parameter is not specified.")
  
  parser$add_argument("--nthreads", type = "integer",
                      help = "number of parallel threads")
  
  parser
}

load_reference = function(ref_h5ad,celltype_col){
  ref = schard::h5ad2list(ref_h5ad)
  ctnames = ref$obs[,celltype_col]
  ctnames = gsub('/','_',ctnames)
  ctnames = factor(ctnames)
  names(ctnames) = colnames(ref$X)
  reference = Reference(ref$X,ctnames)
  reference
}

load_c2l_reference = function(ref_h5ad){
  ref = schard::h5ad2data.frame(ref_h5ad,'varm/means_per_cluster_mu_fg')
  colnames(ref) = sub('means_per_cluster_mu_fg_','',colnames(ref))
  colnames(ref) = gsub('/','_',colnames(ref))
  ref$`_index` = NULL
  ref = as.matrix(ref)
  ref
}

load_visium = function(vis_h5ads){
  # in case of multiple samples it will offset coordinates by imagerow to make sure they do not overlap
  counts = NULL
  coords = NULL
  shift = 0
  for(p in vis_h5ads){
    v = schard::h5ad2list(p,load.obsm = TRUE)
    v$X = v$X[v$var$feature_types=='Gene Expression',]
    rownames(v$obsm$spatial) = colnames(v$X)
    v$obsm$spatial[,1] = v$obsm$spatial[,1] + shift
    shift = max(v$obsm$spatial[,1]) + (max(v$obsm$spatial[,1])-min(v$obsm$spatial[,1]))/10
    counts = cbind(counts,v$X)
    coords = rbind(coords,v$obsm$spatial)  
  }
  visium = SpatialRNA(as.data.frame(coords), counts)
  visium
}

run_rctd = function(reference,visium,ncpu,doublet_mode,precomp_ref){
  if(precomp_ref){
    ctnames = colnames(reference)
    ctnames = factor(ctnames)
    names(ctnames) = colnames(reference)
    reference_mock = Reference(round(reference),ctnames)
    
    myRCTD =  create.RCTD(visium, reference = reference_mock,
                          max_cores = ncpu,CELL_MIN_INSTANCE=1,cell_type_profiles = reference)
  }else{
    myRCTD =  create.RCTD(visium, reference, max_cores = ncpu,CELL_MIN_INSTANCE=1)
  }
  myRCTD = run.RCTD(myRCTD, doublet_mode = doublet_mode) 
  myRCTD
}

main = function(){
  parser = get_parser()
  args = parser$parse_args()
  #args = parser$parse_args(cmdargs)

  print(jsonlite::toJSON(args))
  
  vis = load_visium(args$vis_h5ads)
  precomp_ref = is.null(args$celltype_col)
  
  if(precomp_ref){
    print('Using precomputed c2l reference signatures')
    ref = load_c2l_reference(args$ref_h5ad) 
  }else{
    print('Using RCTD to compute reference signatures')
    ref = load_reference(args$ref_h5ad,args$celltype_col)
  }
  
  rctd = run_rctd(ref,vis,args$nthreads,args$doublet_mode,precomp_ref=precomp_ref)
  
  saveRDS(rctd,paste0(args$outbase,'.rds'))
  write.csv(rctd@results$weights,paste0(args$outbase,'_weights.csv'))
}

main()

# cmdargs = c("--ref_h5ad","/nfs/team283/GBM_LEAP/GBM_LEAP_annotations/cell2location/reference_signatures/run_sc_AT10.h5ad",
#   #"--celltype_col","TME_GBM_granular",
#   "--vis_h5ads", "/nfs/team283/GBM_LEAP/GBM_LEAP_annotations/web_atlas/anndata/AT10-BRA-5-FO-1_1.h5ad",
#   "--outbase", "rctd_out/A10",
#   "--doublet_mode","full",
#   "--nthreads","16")

