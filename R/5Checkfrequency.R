Create_shared_variant_site_table <- function(transmission_pair,NonSyn_or_Syn){
  donor_table=get_sample_frequency_table(transmission_pair@Donor)
  recipient_table=get_sample_frequency_table(transmission_pair@Recipient)
  if(NonSyn_or_Syn=="Synonymous"){
    donor_table=donor_table[donor_table[,6]=="Syn",]
    recipient_table=recipient_table[recipient_table[,6]=="Syn",]
    donor_table=donor_table[,-6]
    recipient_table=recipient_table[,-6]
    }
  if(NonSyn_or_Syn=="Non-Synonymous"){
    donor_table=donor_table[donor_table[,6]=="Non",]
    recipient_table=recipient_table[recipient_table[,6]=="Non",]
    donor_table=donor_table[,-6]
    recipient_table=recipient_table[,-6]
  }
  if(NonSyn_or_Syn=="All"){
    donor_table=donor_table
    recipient_table=recipient_table
    donor_table=donor_table[,-6]
    recipient_table=recipient_table[,-6]
  }
  shared_table=merge(donor_table,recipient_table, by.x=1,by.y=1)
  names(shared_table)=c("variant_ID",names(donor_table)[2:5],names(recipient_table)[2:5])
  row.names(shared_table)=shared_table[,1]
  shared_table=shared_table[,-1]
  return(shared_table)
}

n_used_get <- function(sample_table,threshold){
  sample_table$sum=rowSums(sample_table)
  n = nrow(sample_table[sample_table[,5]>=threshold,])
  return(n)
}


n_unused_get <- function(sample_table,threshold){
  sample_table$sum=rowSums(sample_table)
  n = nrow(sample_table[sample_table[,5]<threshold,])
  return(n)
}

create_log_sample <- function(n_used,n_not_used){
  log_list=list(n_used,n_not_used)
  return(log_list)
}

create_log_transmission_pair <- function(log_table,donor_list,recipient_list,donor_id,recipirnt_id){
  row=data.frame(list(donor_id,recipirnt_id,unlist(donor_list),unlist(recipient_list)))
  return(row)
}

tidy_up_sample_table <- function(sample_table,threshold){
  sample_table$sum=rowSums(sample_table)
  sample_table=sample_table[sample_table[,5]>threshold,]
  for(i in 1:4){
    sample_table[,i][sample_table[,i]==0]=0
  }
  return(sample_table)
}

create_frequency_table <- function(sample_sum_table,error_calling){
  for(i in 1:4){
    sample_sum_table[,i]=sample_sum_table[,i]/sample_sum_table[,5]
    sample_sum_table[,i][sample_sum_table[,i]<error_calling]=0
  }
  return(sample_sum_table)
}

tidy_up_shared_sites_table <- function(shared_table,donor_depth_threshold, recipient_depth_threshold,error_calling){
  Donor=shared_table[,1:4]
  Donor=tidy_up_sample_table(Donor,threshold = donor_depth_threshold)
  recipient=shared_table[,5:8]
  recipient=tidy_up_sample_table(recipient,threshold = recipient_depth_threshold)
  Donor=create_frequency_table(Donor,error_calling = error_calling)
  recipient=create_frequency_table(recipient,error_calling = error_calling)
  tidy_table=merge(Donor,recipient,by.x=0,by.y=0)
  tidy_table=tidy_table[,-11]
  tidy_table=tidy_table[,-6]
  row.names(tidy_table)=tidy_table[,1]
  tidy_table=tidy_table[,-1]
  names(tidy_table)=c(paste0("do.",names(Donor)[1:4]),paste0("re.",names(recipient)[1:4]))
  return(tidy_table)
}




