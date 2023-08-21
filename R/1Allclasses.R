
#'@exportClass sample
setClass("sample",
         slots = list(sample_ID="character",
                      variant_site_table="data.frame"))


#'@exportClass transmission_pair
setClass("transmission_pair",
         slots = list(transmission_pair_ID="character",
                      Donor="sample",
                      Recipient="sample"))


#'@exportClass ob_for_wf
setClass("ob_for_wf",
         slots=list(variant_ids="character",
                    proportion_table="data.frame"))


