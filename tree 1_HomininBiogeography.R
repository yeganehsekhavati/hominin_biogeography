library(paleotree)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(phytools)

wd="~/3. R analyses_final"
setwd(wd)
getwd()

trfn = np(paste(addslash(wd), "mbl1.newick", sep="")) # tree path
geogfn = np(paste(addslash(wd), "Range.txt", sep="")) #range path

tr = read.tree(trfn)
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

max_range_size = 3 #maximum number of areas


####################################################
# Run DEC
####################################################
# Intitialize DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# tree and geographic files
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn

BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = FALSE  # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50  
BioGeoBEARS_run_object$speedup = TRUE          
BioGeoBEARS_run_object$use_optimx = "GenSA"    
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) #
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE 

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "DEC"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}


####################################################
# Run DEC+j
####################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = FALSE  # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50    
BioGeoBEARS_run_object$speedup = TRUE        
BioGeoBEARS_run_object$use_optimx = "GenSA"   
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) 
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE   

# Set up DEC+J model
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Starting values for anagenetic parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DECj"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}


####################################################
# Run DIVALIKE
####################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001   
BioGeoBEARS_run_object$include_null_range = FALSE   # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50   
BioGeoBEARS_run_object$speedup = TRUE        
BioGeoBEARS_run_object$use_optimx = "GenSA"   
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE 

# Set up DIVALIKE model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "DIVALIKE"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}


####################################################
# Run DIVALIKE+j
####################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = FALSE    # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50    
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"  
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) 
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE   

# Set up DIVALIKE+J model
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DIVALIKEj"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}


#######################################################
# Run BAYAREALIKE
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = FALSE   # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50    
BioGeoBEARS_run_object$speedup = TRUE        
BioGeoBEARS_run_object$use_optimx = "GenSA"  
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object) 
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE 

# Set up BAYAREALIKE model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "BAYAREALIKE"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}


#######################################################
# Run BAYAREALIKE+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = FALSE   # set TRUE for models in which null range is included
BioGeoBEARS_run_object$on_NaN_error = -1e50  
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"    
BioGeoBEARS_run_object$num_cores_to_use = 4
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE   

dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Setup BAYAREALIKE+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Setting the maximum of "j" tp 1, as required for BAYAREALIKE+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "BAYAREALIKEj"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resBAYAREALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKEj = res
}


#######################################################
# Model comparison: DEC vs. DEC+J
#######################################################
restable = NULL
teststable = NULL

LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)


#######################################################
# Model comparison: DIVALIKE vs. DIVALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Model comparison: BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
conditional_format_table(stats)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

# Get the results tables 
teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable = put_jcol_after_ecol(restable)
restable

# Save
save(restable, file="restable_v1.Rdata")
load(file="restable_v1.Rdata")

save(teststable, file="teststable_v1.Rdata")
load(file="teststable_v1.Rdata")

write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

restable2 = restable

# AIC
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
conditional_format_table(restable_AIC_rellike)

#  AICcs
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
conditional_format_table(restable_AICc_rellike)


#######################################################
# Plot ancestral states - BAYAREALIKE
#######################################################
pdffn = "Hominin Bayarealike.pdf"
pdf(pdffn, width=30, height=20)

analysis_titletxt ="BAYAREALIKE on Hominins"

results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="text", label.offset=0.3, tipcex=2, statecex=1.6, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)
plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="pie", label.offset=0.3, tipcex=2, statecex=1.6, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

dev.off() 
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
pdffn = "Hominin Bayarealike+j tree1.pdf"
pdf(pdffn, width=6, height=6)

analysis_titletxt ="BAYAREALIKE+j on Hominins"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=FALSE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Plot legend
#######################################################
pdffn = "colors_legend.pdf"
pdf(pdffn, width=12, height=12)

areanames = names(tipranges@df)
areanames

include_null_range = FALSE

states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range)

statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range, split_ABC=FALSE)
statenames

relprobs_matrix = resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
MLprobs = get_ML_probs(relprobs_matrix)
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, returnwhat="states", if_ties="takefirst")

colors_matrix = get_colors_for_numareas(length(areanames))
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=include_null_range)
colors_list_for_states

possible_ranges_list_txt = areas_list_to_states_list_new(areas=areanames,  maxareas=max_range_size, split_ABC=FALSE, include_null_range=include_null_range)
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, colors_list_for_states, MLstates)

legend_ncol=NULL
legend_cex=0.4
colors_legend(possible_ranges_list_txt, colors_list_for_states, legend_ncol=legend_ncol, legend_cex=legend_cex)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#######################################################
#Get the probability for each node 
#######################################################
relprobs_matrix = resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node
get_ML_probs(relprobs_matrix, unlist_TF = TRUE)
resBAYAREALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node


#######################################################
# stochastic mapping
#######################################################
model_name = "BAYAREALIKE"
res = resBAYAREALIKE

BSM_inputs_fn = "BSM_inputs_file.Rdata"
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  load(BSM_inputs_fn)
} 

names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=100, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  
  load(file="RES_clado_events_tables.Rdata")
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} 

clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = FALSE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 3

############################################
# Painting one stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]


############################################
# Saving the result as a PDF
############################################
pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
pdf(pdffn, width=12, height=12)

master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)
plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


#######################################################
# Save all stochastic maps 
#######################################################
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range = FALSE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = stratified

pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v.pdf")
pdf(pdffn, width=30, height=20)

nummaps_goal = 100
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"),  plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range,label.offset=0.2, tipcex=1.8, statecex=1.6, splitcex=0.6,plotsplits=FALSE)  #do not show splits
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(model_name, "_", tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } 
} 
cat("...done.\n")

#######################################################
# Check that ancestral state
#######################################################
library(MultinomialCI)   
check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)


## Adopted from http://phylo.wikidot.com/biogeobears
