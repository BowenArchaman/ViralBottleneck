#########################################################################################################################
#create a transmission object
#using one transmission pair to create transmission object
#'@importFrom methods new
#'@importFrom pbapply pbapply
Create_object_using_one_pair <- function(row){
  Donor=read.table(paste0(row[1],".csv"),header = TRUE,,sep = ",")
  Donor=produce_variant_id(Donor)
  Recipient=read.table(paste0(row[2],".csv"),header = TRUE,,sep = ",")
  Recipient=produce_variant_id(Recipient)
  Recipient_name=row[2]
  sample_donor=new("sample",sample_ID=row[1],variant_site_table=Donor)
  sample_recipient=new("sample",sample_ID=row[2],variant_site_table=Recipient)
  transmission_pair_sample = new("transmission_pair",transmission_pair_ID=row[3],
                                 Donor=sample_donor,Recipient=sample_recipient)
  return(transmission_pair_sample)
}

#using transmission pair table to create transmission object
Create_object_using_table <- function(transmission_pairs_table){
  transmission_pairs_table=produce_transmission_pairs(transmission_pairs_table)
  transmission_object=list(pbapply(transmission_pairs_table,1,Create_object_using_one_pair))
  return(transmission_object)
}

#########################################################################################################################
#Create transmission object containing above every thing
#' Create Transmission Object
#'
#' This function could check the input files and catch the files in working directory according to transmission
#'  pairs table to create a transmission object.
#'
#'
#' @param transmission_pairs_table is a dataframe only contatining 2 columns, Donor and Recipient
#'
#' @return Transmission_object, a environment containing all the infomation of transmission pairs
#' @examples
#' mytransmission_object = CreateTransmissionObject(transmission_pairs_table)
#'
#' @export
CreateTransmissionObject <- function(transmission_pairs_table){
  if(Check_everything(transmission_pairs_table)){
    transmisson_object=Create_object_using_table(transmission_pairs_table)
    transmission_pairs_number= length(transmisson_object[[1]])
    sample_number=length(union(transmission_pairs_table[,1],transmission_pairs_table[,2]))
    message(paste("There are",transmission_pairs_number,"transmission pairs"))
    message(paste("Total",sample_number,"samples take part in creating transmission object"))
    transmisson_object=transmisson_object[[1]]
    return(transmisson_object)
  }
}


