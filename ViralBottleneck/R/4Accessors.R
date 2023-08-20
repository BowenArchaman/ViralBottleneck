
#######################################################################################################
setGeneric("get_shared_variant_site_number",function(x){#
  standardGeneric("get_shared_variant_site_number")
})

setGeneric("get_transmission_pair_ID",function(x){#
  standardGeneric("get_transmission_pair_ID")
})

setGeneric("get_sample_frequency_table",function(x){#
  standardGeneric("get_sample_frequency_table")
})


setGeneric("get_sample_from_transmission_pairs",function(x){
  standardGeneric("get_sample_from_transmission_pairs")
})

#######################################################################################################
setMethod("get_transmission_pair_ID",c(x = "transmission_pair"),
          function(x){
            id=x@transmission_pair_ID
            return(id)
          })

setMethod("get_shared_variant_site_number",c(x = "transmission_pair"),function(x){
  Donor_id=x@Donor@variant_site_table$variant_ID
  Recipient_id=x@Recipient@variant_site_table$variant_ID
  number=length(intersect(Donor_id,Recipient_id))
  return(number)
})

setMethod("get_sample_frequency_table",c(x="sample"),function(x){
  variant_table=x@variant_site_table
  variant_table = data.frame(cbind(variant_table[,8],variant_table[,3:7]))
  return(variant_table)
})

#######################################################################################################
get_transmission_pairs<-function(transmission_ob){
 ids=as.character(lapply(transmission_ob, get_transmission_pair_ID))
 return(ids)
}

get_shared_site_number<-function(transmission_ob){
 numbers=as.numeric(lapply(transmission_ob, get_shared_variant_site_number))
 return(numbers)
}

