#' @importFrom pbapply pblapply
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom rmutil dbetabinom
#' @importFrom rmutil pbetabinom
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme_bw

#Tidy table, check coverage, error_calling, depth_threshold(donor,recipient) one transmission pair
Before_calculation_table<-function(one_pair,donor_depth_threshold, recipient_depth_threshold,error_calling,NonSyn_or_Syn){
  table=Create_shared_variant_site_table(transmission_pair = one_pair,NonSyn_or_Syn=NonSyn_or_Syn)
  tidy_table=tidy_up_shared_sites_table(table,donor_depth_threshold = donor_depth_threshold,recipient_depth_threshold = recipient_depth_threshold,error_calling = error_calling)
  return(tidy_table)
}


create_max_f <- function(shared_site,tidy_table){
  mix=merge(tidy_table,shared_site,by.x=0,by.y=0)
  row.names(mix)=mix[,1]
  mix=mix[,-1]
  donor=mix[,1:4]
  donor$d_max=apply(donor,1,max)
  max_table=cbind.data.frame(donor,mix[,5:16])
  max_table=max_table[,-(10:13)]
  return(max_table)
}

find_dominant_not_mutated<-function(row){
  if(which.max(row[1:4])==which.max(row[5:8])){
    return(TRUE)
  }
  else{return(FALSE)}
}

find_dominant_in_recipient <- function(row){
  index=which.max(row[1:4])
  dominant=row[index+4]
  return(dominant)
}

find_variant_in_recipient <- function(row){
  index=which(row[1:4]==row[10])
  variant=row[index+4]
  return(variant)
}

Create_matrix_for_biallelic <- function(shared_table,tidy_shared_table,variant_calling){
  mix=create_max_f(shared_table,tidy_shared_table)
  donor=mix[,1:4]
  sort=t(apply(donor,1,sort,decreasing=TRUE))#sort to find dominant and variant
  sort=sort[sort[,2]>variant_calling,]
  sort=sort[sort[,2]!=sort[,3],]#check the same variants
  sort=sort[sort[,3]<=variant_calling,]#filter non-biallelic
  res=merge(mix,sort,by.x=0,by.y=0)
  row.names(res)=res[,1]
  res=res[,2:18]
  res[,"V3"]=apply(cbind.data.frame(res[,1:4],res[,6:9],res[,14:15]),1,find_dominant_in_recipient)
  res[,"V4"]=apply(cbind.data.frame(res[,1:4],res[,6:9],res[,14:15]),1,find_variant_in_recipient)
  res$V5=apply(cbind.data.frame(res[,1:4],res[,10:13],res[,14:15]),1,find_dominant_in_recipient)
  res$V6=apply(cbind.data.frame(res[,1:4],res[,10:13],res[,14:15]),1,find_variant_in_recipient)
  return(res)
}



#filter the variants which would not be distinguished from errors.
Prepared_matrix_for_methods <- function(shared_sites_table,tidy_sites_table,variant_calling){
  res=Create_matrix_for_biallelic(shared_sites_table,tidy_sites_table,variant_calling)
  res_info=res[,14:19]
  names(res_info)=c("do.dom","do.subdom","re.dom","re.subdom","re.dom_reads","re.subdom_reads")
  prepared_matrix=cbind.data.frame(res[,1:13],res_info)
  return(prepared_matrix)
}


Convert_to_Approxmate_method_matrix <- function(prepared_matrix){
  App_matrix=cbind.data.frame(prepared_matrix[,15],prepared_matrix[,17])
  return(App_matrix)
}

Convert_to_Exact_method_matrix<- function(prepared_matrix){
  matrix=cbind.data.frame(prepared_matrix[,15],prepared_matrix[,18:19])
  matrix$sum=rowSums(matrix[,2:3])
  Exact_matix=cbind.data.frame(matrix[,1],matrix[,3:4])
  return(Exact_matix)
}

Create_variant_identificatin_forKL <- function(shared_table,tidy_shared_table,variant_calling){
  mix=create_max_f(shared_table,tidy_shared_table)
  donor=mix[,1:4]
  #if error filtering is delete it should be delete
  for(i in 1:8){
    tidy_shared_table[,i][tidy_shared_table[,i]<variant_calling]=0
  }
  sort=t(apply(donor,1,sort,decreasing=TRUE))#sort to find dominant and variant
  sort=sort[sort[,2]>variant_calling,] #filtered the no-variation sites
  var_sites=merge(sort[,1:2],tidy_shared_table,by="row.names")
  row.names(var_sites)=var_sites[,1]
  var_sites=var_sites[,-(1:3)]
  return(var_sites)
}



find_fixed_variant_app_onerow <- function(row,variant_calling){
  if(row[2]<=(1-variant_calling)){
    row=row
    
  }
  else{
    row[1]=1-row[1]
    row[2]=1-row[2]
  }
  return(row)
}

find_fixed_variant_app<- function(table,variant_calling){
  table=as.data.frame(t(apply(table,1,find_fixed_variant_app_onerow,variant_calling=variant_calling)))
  return(table)
}



find_fixed_variant_exact_onerow <- function(row,variant_calling){
  if(row[2]<=(1-variant_calling)*row[3]){
    row=row
    
  }
  else{
    row[1]=1-row[1]
    row[2]=row[3]-row[2]
  }
  return(row)
}
find_fixed_variant_exact<- function(table,variant_calling){
  table=as.data.frame(t(apply(table,1,find_fixed_variant_exact_onerow,variant_calling=variant_calling)))
  return(table)
}


find_confidence_interval <- function(final_vector,Nbmin){
  final_vector=as.numeric(final_vector[1,])
  l_OK = is.finite(final_vector)
  if(FALSE %in% l_OK){
    lowest_index=Nbmin
  }
  else{lowest_index=Nbmin-1}
  final_vector=final_vector[l_OK]
  max_value=max(final_vector)
  max_index=which.max(final_vector)+lowest_index
  height=max_value- qchisq(0.95,df=1)/2
  CI=final_vector[final_vector>=height]
  CI_low=which(final_vector==CI[1]) + lowest_index
  CI_low=min(CI_low)
  CI_high=which(final_vector==CI[length(CI)])+lowest_index
  CI_high=max(CI_high)
  if(is.na(CI[1])){
    CI_low=max_index
    CI_high=max_index+1
  }
  CI_list=list(CI_low,CI_high,max_index)
  return(CI_list)
}
######################################method 
#KL 
KL_shared_table <- function(shared_site_table){
  sum=0
  for(i in 1:4){
    q=shared_site_table[,i]
    p=shared_site_table[,i+4]
    p[p == 0] <- 1e-9
    q[q == 0] <- 1e-9
    ratio=p/q
    add=p*log2(ratio)
    sum=sum+add
  }
  return(sum)
}

Effective_bottleneck_size_KL <- function(k,shared_site_table){
  n_s=nrow(shared_site_table)
  KL=sum(KL_shared_table(shared_site_table))
  one_val_likelihood=-(k*KL)+(n_s/2)*log(k)
  return(one_val_likelihood)
}

Range_function_KL<-function(shared_site_table,Nbmin,Nbmax){
  likelihood_vector=pblapply(Nbmin:Nbmax,Effective_bottleneck_size_KL,shared_site_table=shared_site_table)
  print(paste(nrow(shared_site_table),"sites were used in calculation"))
  final_likelihood_vector=data.frame(likelihood_vector)
  name_v=Nbmin:Nbmax
  class(name_v)="character"
  names(final_likelihood_vector)=name_v
  return(final_likelihood_vector)
}

#Beta-binomial Approximate verison
one_Nbval_function_Approximate<- function(k,table,variant_calling){
  table=find_fixed_variant_app(table,variant_calling)
  table=subset(table,table[,1]>=variant_calling)
  present=table[table[,2]>=variant_calling,]
  absent=table[table[,2]<variant_calling,]
  likelihood_vector_present=numeric(nrow(present))
  likelihood_vector_absent=numeric(nrow(absent))
  for(i in 0:k){
    P1=dbeta(present[,2], i, k-i)
    pbinVd1=dbinom(i, size=k, prob=present[,1])
    add1=P1*pbinVd1
    #print(is.infinite(add1))
    likelihood_vector_present=likelihood_vector_present+add1
    #print(likelihood_vector_present)
  }
  
  for(j in 0:k){
    P2=pbeta(variant_calling, j, k-j)
    pbinVd2=dbinom(j, size=k, prob=absent[,1])
    add2=P2*pbinVd2
    #print(is.infinite(add2))
    likelihood_vector_absent=likelihood_vector_absent+add2
    #print(likelihood_vector_absent)
  }
  sum=sum(log(likelihood_vector_present))+sum(log(likelihood_vector_absent))
  return(sum)
}

Range_function_Approximate<-function(variant_calling,table,Nbmin,Nbmax){
  res_list=pblapply(Nbmin:Nbmax,one_Nbval_function_Approximate,variant_calling=variant_calling,table=table)
  print(paste(nrow(table),"sites were used in calculation"))
  final_likelihood_vector=data.frame(res_list)
  name_v=Nbmin:Nbmax
  class(name_v)="character"
  names(final_likelihood_vector)=name_v
  return(final_likelihood_vector)
}

#Beta-binomial exact verison
one_Nbval_function_Exact<- function(k,table,variant_calling){
  table=find_fixed_variant_exact(table,variant_calling)
  table=subset(table,table[,1]>=variant_calling)
  present=table[table[,2]>=variant_calling*table[,3],]
  absent=table[table[,2]<variant_calling*table[,3],]
  likelihood_vector_present=numeric(nrow(present))
  likelihood_vector_absent=numeric(nrow(absent))
  if(nrow(present) != 0){
    for(i in 0:k){
      alpha=i
      Beta=(k-i)
      if(alpha==0){alpha=10^-9}
      if(Beta==0){Beta=10^-9}
      m=alpha / (alpha+Beta)
      s=(alpha+Beta)
      pbinVd1=dbinom(i, size=k, prob=present[,1])
      add1=dbetabinom(present[,2],present[,3], m, s, log = FALSE)*pbinVd1
      likelihood_vector_present=likelihood_vector_present+add1
    }
  }
  if(nrow(absent) != 0){
    for(j in 0:k){
      alpha=j
      Beta=(k-j)
      if(alpha==0){alpha=10^-9}
      if(Beta==0){Beta=10^-9}
      m=alpha / (alpha+Beta)
      s=(alpha+Beta)
      pbinVd2=dbinom(j, size=k, prob=absent[,1])
      add2=pbetabinom(floor(variant_calling*absent[,3]),absent[,3], m, s)*pbinVd2
      likelihood_vector_absent=likelihood_vector_absent+add2
    }
  }
  sum=sum(log(likelihood_vector_present))+sum(log(likelihood_vector_absent))
  return(sum)
}

Range_function_Exact<-function(variant_calling,table,Nbmin,Nbmax){
  res_list=pblapply(Nbmin:Nbmax,one_Nbval_function_Exact,variant_calling=variant_calling,table=table)
  print(paste(nrow(table),"sites were used in calculation"))
  final_likelihood_vector=data.frame(res_list)
  name_v=Nbmin:Nbmax
  class(name_v)="character"
  names(final_likelihood_vector)=name_v
  return(final_likelihood_vector)
}

#Presence-Absence
one_Nbval_function_preOrabsent <- function(k,table,variant_calling){
  present=table[table[,2]>=variant_calling,]
  absent=table[table[,2]<variant_calling,]
  likelihood_vector_present=numeric(nrow(present))
  likelihood_vector_absent=numeric(nrow(absent))
  likelihood_vector_present=log(1-(1-present[,1])^k)
  likelihood_vector_absent=k*log(1-absent[,1])
  one_val_likelihood=sum(likelihood_vector_absent)+sum(likelihood_vector_present)
  return(one_val_likelihood)
}

Range_function_preOrabsent <- function(variant_calling,table,Nbmin,Nbmax){
  likelihood_list=pblapply(Nbmin:Nbmax,one_Nbval_function_preOrabsent,table=table,variant_calling=variant_calling)
  print(paste(nrow(table),"sites were used in calculation"))
  final_likelihood_vector=data.frame(likelihood_list)
  name_v=Nbmin:Nbmax
  class(name_v)="character"
  names(final_likelihood_vector)=name_v
  return(final_likelihood_vector)
}

#binomial method
one_Nbval_function_binomial <- function(k,table,variant_calling){
  table=find_fixed_variant_exact(table,variant_calling)
  present=table[table[,2]>=variant_calling*table[,3],]
  absent=table[table[,2]<variant_calling*table[,3],]
  likelihood_vector_present=numeric(nrow(present))
  likelihood_vector_absent=numeric(nrow(absent))
  for(i in 0:k){
    pbinVd1=dbinom(i, size=k, prob=present[,1])
    add1=dbinom(present[,2],size = present[,3],prob = i/k)*pbinVd1
    likelihood_vector_present=likelihood_vector_present+add1
  }
  for(j in 0:k){
    pbinVd2=dbinom(j,size = k,prob = absent[,1])
    add2=pbinom(floor(variant_calling*absent[,3]),size = absent[,3],prob = j/k)*pbinVd2
    likelihood_vector_absent=likelihood_vector_absent+add2
  }
  sum=sum(log(likelihood_vector_present))+sum(log(likelihood_vector_absent))
  return(sum)
}

Range_function_binomial <- function(variant_calling,table,Nbmin,Nbmax){
  likelihood_list=pblapply(Nbmin:Nbmax,one_Nbval_function_binomial,table=table,variant_calling=variant_calling)
  print(paste(nrow(table),"sites were used in calculation"))
  final_likelihood_vector=data.frame(likelihood_list)
  name_v=Nbmin:Nbmax
  class(name_v)="character"
  names(final_likelihood_vector)=name_v
  return(final_likelihood_vector)
}

################################################################################################################
#Wrighter-Fisher method
Get_variant_ids_from_one_pair <- function(one_ob,donor_depth_threshold, recipient_depth_threshold,error_calling,NonSyn_or_Syn){
  s_t=Create_shared_variant_site_table(one_ob,NonSyn_or_Syn = NonSyn_or_Syn)
  t_t=tidy_up_shared_sites_table(s_t,donor_depth_threshold, recipient_depth_threshold,error_calling)
  f_t=Prepared_matrix_for_methods(s_t,t_t,error_calling)
  ids=row.names(f_t)
  f_t=f_t[,-6]
  f_t=f_t[,-5]
  f_t=f_t[,-3]
  before_wf_ob=new("ob_for_wf",variant_ids=ids,proportion_table=f_t)
  return(before_wf_ob)
}

Get_variant_ids_from_objects <- function(transmission_ob,donor_depth_threshold, recipient_depth_threshold,error_calling,NonSyn_or_Syn){
  processed_ob=pblapply(transmission_ob,Get_variant_ids_from_one_pair,donor_depth_threshold=donor_depth_threshold,recipient_depth_threshold=recipient_depth_threshold,error_calling=error_calling,NonSyn_or_Syn=NonSyn_or_Syn)
  return(processed_ob)
}

Get_shared_sites_ID_cross_samples <- function(processed_ob){
  intersect_id=processed_ob[[1]]@variant_ids
  for(i in 2:length(processed_ob)){
    intersect_id=intersect(intersect_id,processed_ob[[i]]@variant_ids)
  }
  return(intersect_id)
}

Apply_sd_table <- function(table){
  sd_vector=apply(table,1,sd)
  return(sd_vector)
}

one_ob_wf_process <- function(wf_ob,intersect_id){
  id_df=data.frame(intersect_id)
  merge_table=merge(id_df,wf_ob@proportion_table,by.x=1,by.y=0)
  row.names(merge_table)=merge_table[,1]
  proportion_matrix=merge_table[,-1]
  proportion_matrix$var=(proportion_matrix[,2]-proportion_matrix[,3])
  new_ob=new("ob_for_wf",variant_ids=intersect_id,proportion_table=proportion_matrix)
  return(new_ob)
}

obs_wf_process <- function(wf_obs,intersect_ids){
  new_obs=lapply(wf_obs,one_ob_wf_process,intersect_id=intersect_ids)
  return(new_obs)
}

return_value_table <- function(new_obs,n){
   starter=new_obs[[1]]@proportion_table
   table=starter[,n]
   for(i in 2:length(new_obs)){
     value=new_obs[[i]]@proportion_table
     table=cbind.data.frame(table,value[,n])
   }
   return(table)
}

return_var_with_variance <- function(P_table,P1_table,signal){
  sd_p=Apply_sd_table(P_table)
  sd_p1=Apply_sd_table(P1_table)
  if(signal=="+"){
    P_table=P_table+sd_p
    P1_table=P1_table+sd_p1
    var_table=P_table-P1_table
    var_v=(Apply_sd_table(var_table))^2
    return(var_v)
  }
  if(signal=="-"){
    P_table=P_table-sd_p
    P1_table=P1_table-sd_p1
    var_table=P_table-P1_table
    var_v=(Apply_sd_table(var_table))^2
    return(var_v)
  }
}

Calculate_Bottleneck_size_wrighter_fisher <- function(p,q,var){
  p_q_var_table=cbind.data.frame(p,q,var)
  nb_v=(p_q_var_table[,1]*p_q_var_table[,2]) / (2*p_q_var_table[,3])
  return(nb_v)
}

Collect_values_in_table_for_wf <- function(transmission_obs,donor_depth_threshold, recipient_depth_threshold,error_calling,NonSyn_or_Syn=NonSyn_or_Syn){
  new_obs=Get_variant_ids_from_objects(transmission_obs,donor_depth_threshold, recipient_depth_threshold,error_calling,NonSyn_or_Syn)
  shared_ids=Get_shared_sites_ID_cross_samples(new_obs)
  ids_df=data.frame(shared_ids)
  new_ob_processed=obs_wf_process(new_obs,shared_ids)
  Q_t=return_value_table(new_ob_processed,1)
  P_t=return_value_table(new_ob_processed,2)
  P1_t=return_value_table(new_ob_processed,3)
  var_t=return_value_table(new_ob_processed,4)
  Q_mean=rowMeans(Q_t)
  P_mean=rowMeans(P_t)
  P1_mean=rowMeans(P1_t)
  #var_mean=rowMeans(var_t)
  var_mean=(Apply_sd_table(var_t))^2
  sd_Q=Apply_sd_table(Q_t)
  sd_P=Apply_sd_table(P_t)
  Q_minus=rowMeans(Q_t-sd_Q)
  Q_add=rowMeans(Q_t+sd_Q)
  P_minus=rowMeans(P_t-sd_P)
  P_add=rowMeans(P_t+sd_P)
  var_minus=return_var_with_variance(P_t,P1_t,"-")
  var_add=return_var_with_variance(P_t,P1_t,"+")
  N_t=round(Calculate_Bottleneck_size_wrighter_fisher(P_mean,Q_mean,var_mean))
  N_t_add=round(Calculate_Bottleneck_size_wrighter_fisher(P_add,Q_minus,var_add))
  N_t_minus=round(Calculate_Bottleneck_size_wrighter_fisher(P_minus,Q_add,var_minus))
  final_table=cbind.data.frame(shared_ids,P_mean,Q_mean,var_mean,N_t,P_add,Q_minus,var_add,N_t_add,P_minus,Q_add,var_minus,N_t_minus)
  ids=t(data.frame(strsplit(final_table[,1],"-")))
  final_table=cbind.data.frame(ids[,1],ids[,2],final_table[,2:ncol(final_table)])
  names(final_table)=c("segments","position","P_mean","Q_mean","var","Nb","P+eps","Q-eps","var+eps","Nb+eps","P-eps","Q+eps","var-eps","Nb-eps")
  row.names(final_table)=1:nrow(final_table)
  final_table[final_table[,14]<=0,14]=NA
  return(final_table)
}

#####################################################################################################################

#' Summary_ob
#'
#' @param transmission_ob is a object created by `CreateTransmissionObject` function
#' @param save is a logical value to determine whether the table need to be saved in path or not.
#' @param file_name is the name of file. It only would be used when `save=TRUE`.
#'
#' @return would return a csv file containing the transmission pairs: donors, recipients and shared variant sites.

#' @examples
#' #The example can not be run directly
#' summary_table = Summary_ob(transmission_object)

#' @export
Summary_ob <- function(transmission_ob,save=FALSE,file_name=NULL){
  transmission_info=c("donors","recipients","number of shared sites")
  content=data.frame(cbind(get_transmission_pairs(transmission_ob),get_shared_site_number(transmission_ob)))
  ids=t(data.frame(strsplit(content[,1],"-")))
  summary_table=cbind.data.frame(ids[,1],ids[,2],content[,2:ncol(content)])
  row.names(summary_table)=1:nrow(summary_table)
  names(summary_table)=transmission_info
  if(save==TRUE && is.null(file_name)==TRUE){stop("File name is null, please input the file name.")}
  if(save==TRUE && is.null(file_name)==FALSE){
    check_file(file_name)
    write.csv(summary_table,file = file_name)
    }
  return(summary_table)
}


#############################################################################################
find_index<-function(ele,r_f){
  index=which(r_f==ele)
  return(index)
}
get_from_ob<- function(index,get_l){
  if(length(index)==0){
    return(NA)
  }
  else{
    ele=get_l[[as.numeric(index)]]
    return(ele)
  }
}

tidy_get <- function(user_ids,ob_list,ob){
  res=lapply(user_ids,find_index,r_f=ob_list)
  get_res=lapply(res,get_from_ob,get_l=ob)
  get_res=get_res[!is.na(get_res)]
  return(get_res)
}


create_difference <- function(user_ids,ob_ids){
  diff=user_ids[!(user_ids %in% ob_ids)]
  return(diff)
}
#############################################################################################
message_diff_ele <- function(ele){
  message(paste(ele,"is missing in object, please check"))
}

message_diff <- function(diff_list){
  lapply(diff_list, message_diff_ele)
}
# subset transmission_pair object according to transmission pairs which is user input
#############################################################################################
get_from_in_object <- function(transmission_ob,transmission_pairs){
  ob_ids=get_transmission_pairs(transmission_ob)
  pro_table=produce_transmission_pairs(transmission_pairs)
  user_ids=union(pro_table[,3],pro_table[,3])
  get=tidy_get(user_ids,ob_ids,transmission_ob)
  diff_l=create_difference(user_ids,ob_ids)
  message_diff(diff_l)
  return(get)
}



plot_likelihood_function <- function(id,final_likelihood_vector,confidence_res,method){
  plot_name=paste0(id,"_",method,"_","plot.png")
  table_final_likelihood_vector=rbind.data.frame(final_likelihood_vector,names(final_likelihood_vector))
  x=t(table_final_likelihood_vector)
  x=data.frame(x)
  png(plot_name)
  ggp<-ggplot(x,aes(x=as.numeric(X2),y=as.numeric(X1)))+
    xlab("Nb")+ylab("log-likelihood")+geom_line(col="darkorange2",linewidth=1.5)+
    geom_vline(xintercept=confidence_res[[3]],col="red",linewidth=0.5) +
    annotate("rect",xmin=confidence_res[[1]], xmax=confidence_res[[2]], ymin=-Inf, ymax=Inf,alpha=.2,fill="yellow")+
    theme_bw()
  print(ggp)
  dev.off()
}
#############################################################################################
#check directory exist
check_dir <- function(name_dir){
  if(file.exists(name_dir)){
    setwd(file.path(name_dir))
  } else {
    dir.create(name_dir)
    setwd(name_dir)
  }
}

check_file <- function(name_file){
  if(file.exists(name_file)){
    file.remove(name_file)
  }
  else{
  }
}

log_one_pair<- function(one_ob,donor_threshold,recipient_threshold,NonSyn_or_Syn,error_calling,variant_calling){
  transmisson_id = one_ob@transmission_pair_ID
  shared_table=Create_shared_variant_site_table(one_ob,NonSyn_or_Syn=NonSyn_or_Syn)
  tidy_table = Before_calculation_table(one_ob,donor_depth_threshold = donor_threshold, recipient_depth_threshold=recipient_threshold,error_calling=error_calling,NonSyn_or_Syn=NonSyn_or_Syn)
  table=Prepared_matrix_for_methods(shared_table,tidy_table,variant_calling)
  donor_used=nrow(table)
  donor_unused=nrow(one_ob@Donor@variant_site_table)-nrow(table)
  recipient_used=nrow(table)
  recipient_unused=nrow(one_ob@Recipient@variant_site_table)-nrow(table)
  devide=strsplit(transmisson_id,split = "-")
  row=list(devide[[1]][1],devide[[1]][2],donor_used,donor_unused,recipient_used,recipient_unused)
  return(row)
}

log_main_function<-function(ob,donor_threshold,recipient_threshold,NonSyn_or_Syn,error_calling,variant_calling){
  main_t=lapply(ob,log_one_pair,donor_threshold=donor_threshold,recipient_threshold=recipient_threshold,NonSyn_or_Syn=NonSyn_or_Syn,error_calling=error_calling,variant_calling=variant_calling)
  main_t=do.call(rbind.data.frame,main_t)
  names(main_t)=c("donor","recipient","donor_used","donor_unused","recipient_used","recipient_unused")
  return(main_t)
}

# One transmission pair
#############################################################################################


one_transmission_pair_process <- function(one_pair,method,donor_depth_threshold,recipient_depth_threshold,error_calling,plot,log,log_dir_name,log_table,variant_calling,Nbmin,Nbmax,NonSyn_or_Syn){
  con=NA
  transmisson_id = one_pair@transmission_pair_ID
  shared_table=Create_shared_variant_site_table(one_pair,NonSyn_or_Syn=NonSyn_or_Syn)
  tidy_table = Before_calculation_table(one_pair,donor_depth_threshold = donor_depth_threshold, recipient_depth_threshold=recipient_depth_threshold,error_calling=error_calling,NonSyn_or_Syn=NonSyn_or_Syn)
  table=Prepared_matrix_for_methods(shared_table,tidy_table,variant_calling)
  if(log==TRUE){
   write.csv(table,paste0(transmisson_id,"_log.csv"),row.names = FALSE) 
  }
  
  if(method == "KL"){
    KL_tidy_table=Create_variant_identificatin_forKL(shared_table, tidy_table,variant_calling)
    v=Range_function_KL(KL_tidy_table,Nbmin,Nbmax)
    res=find_confidence_interval(v,Nbmin=Nbmin)
    Nb=list(res[[3]],res[[1]],res[[2]])
    if(plot==TRUE){
      dir_name=paste0(transmisson_id,"_plot")
      check_dir(dir_name)
      plot_likelihood_function(transmisson_id,v,res,method)
      setwd('..')
    }
    return(Nb)
  }

  if(method=="Beta_binomial_Approximate"){
    matrix_app=Convert_to_Approxmate_method_matrix(table)
    v=Range_function_Approximate(variant_calling=variant_calling ,table=matrix_app,Nbmin = Nbmin,Nbmax=Nbmax)
    res=find_confidence_interval(v,Nbmin=Nbmin)
    Nb=list(res[[3]],res[[1]],res[[2]])
    if(plot==TRUE){
      dir_name=paste0(transmisson_id,"_plot")
      check_dir(dir_name)
      plot_likelihood_function(transmisson_id,v,res,method)
      setwd('..')
    }
    return(Nb)
  }

  if(method=="Beta_binomial_Exact"){
    matrix_exact=Convert_to_Exact_method_matrix(table)
    v=Range_function_Exact(table = matrix_exact,variant_calling=variant_calling,Nbmin = Nbmin,Nbmax = Nbmax)
    res=find_confidence_interval(v,Nbmin=Nbmin)
    Nb=list(res[[3]],res[[1]],res[[2]])
    if(plot==TRUE){
      dir_name=paste0(transmisson_id,"_plot")
      check_dir(dir_name)
      plot_likelihood_function(transmisson_id,v,res,method)
      setwd('..')
    }
    return(Nb)
  }
  if(method=="Presence-Absence"){
    matrix_PA=Convert_to_Approxmate_method_matrix(table)
    v=Range_function_preOrabsent(variant_calling = variant_calling,table = matrix_PA,Nbmin = Nbmin,Nbmax = Nbmax)
    res=find_confidence_interval(v,Nbmin=Nbmin)
    Nb=list(res[[3]],res[[1]],res[[2]])
    if(plot==TRUE){
      dir_name=paste0(transmisson_id,"_plot")
      check_dir(dir_name)
      plot_likelihood_function(transmisson_id,v,res,method)
      setwd('..')
    }
    return(Nb)
  }
  if(method=="Binomial"){
    matrix=Convert_to_Exact_method_matrix(table)
    v=Range_function_binomial(variant_calling=variant_calling,table = matrix,Nbmin = Nbmin,Nbmax=Nbmax)
    res=find_confidence_interval(v,Nbmin=Nbmin)
    Nb=list(res[[3]],res[[1]],res[[2]])
    if(plot==TRUE){
      dir_name=paste0(transmisson_id,"_plot")
      check_dir(dir_name)
      plot_likelihood_function(transmisson_id,v,res,method)
      setwd('..')
    }
    return(Nb)
  }
}
################################################################################################




#transmission bottleneck size calculation
#' Calculate transmission bottleneck size
#'
#' @param transmission_ob is a transmission object created by,`CreateTransmissionObject`.
#' @param method is characters containing `KL` method, `Presence-Absence` method, `Binomial` method, `Beta_binomial_Approximate` method, `Beta_binomial_Exact` method and `Wright-Fisher` method
#' @param plot is a logical value to determine plot the likelihood or not. Each plot would be stored in an individual folder named by transmission id.
#' @param show_table is a logical value to determine output a result table as `csv` format.
#' @param transmission_pairs is a dataframe which is the subset of the transmission pairs table during the object creation.
#' @param donor_depth_threshold is a integer to filter variants in donor with the lower sequencing coverage
#' @param recipient_depth_threshold is a integer to filter variants in recipient with the lower sequencing coverage
#' @param error_filtering is a numeric and filter. The proportion of the variant sites lower than error calling means the variant is difficult to distinguish from PCR errors. Then, the sites with lower proportion than error calling would be filtered.
#' @param log is a logical value to show the log or not. The log would be stored in an individual folder.
#' @param variant_calling is a parameter used for `Presence-absence`method, `Binomial` method, `Beta_binomial` method to determine the variants is absent in recipient or not. Proportion of variant site lower than variant calling means this site is absent in recipient.
#' @param Nbmin is a integer representing the minimum number in likelihood range.
#' @param Nbmax is a integer representing the maximum number in likelihood range.
#' @param NonSyn_or_Syn is a filter to keep the variant sites which meet the requirements. If user need the mutation of variant sites are synonymous mutation, then `NonSyn_or_Syn="Synonymous"`. If user need the mutation of variant sites are non-synonymous mutation, then `NonSyn_or_Syn="Non-Synonymous"`. The default is "All"


#' @return Bottleneck size table containing the transmission pairs, bottleneck size, the lower confidence interval and higher confidence interval

#' @examples
#' #The example can not be run directly
#' #Use different methods for estimation
#' #use KL method
#'KL_table = Bottleneck_size_calculation(transmission_object, method="KL")
#' #use presence-absence method
#'PA_table = Bottleneck_size_calculation(transmission_object, method="Presence-Absence")
#' #use Binomial method
#'Bi_table = Bottleneck_size_calculation(transmission_object, method="Binomial")
#' #use Beta-binomial approximate method
#'BB_App_table = Bottleneck_size_calculation(transmission_object, method="Beta_binomial_Approximate")
#' #use Beta-binomial approximate method
#'BB_Exact_table = Bottleneck_size_calculation(transmission_object, method="Beta_binomial_Exact")
#' #use wright-fisher method
#'BB_WF_table = Bottleneck_size_calculation(transmission_object, method="Wright-Fisher")


#' @references
#' Emmett, K. J., Lee, A., Khiabanian, H., & Rabadan, R. (2015) 
#' High-resolution genomic surveillance of 2014 Ebolavirus using shared 
#' subclonal variants. \emph{PLOS Currents Outbreaks} \bold{7}, 
#' ecurrents.outbreaks.\cr 
#' Sacristán, S., Malpica, J. M., Fraile, A., & García-Arenal, F. (2003)
#' Estimation of population bottlenecks during systemic movement of
#' tobacco mosaic virus in tobacco plants. \emph{Journal of Virology}
#' \bold{77}(18), 9906–9911.\cr
#' Poon, L. L. M., Song, T., Rosenfeld, R., Lin, X., Rogers, M. B.,
#' Zhou, B., Sebra, R., Halpin, R., Guan, Y., Twaddle, A., DePasse, J.,
#' Stockwell, T., Wentworth, D., Holmes, E., Greenbaum, B., Peiris, J. S. M.,
#' Cowling, B. J., & Ghedin, E. (2016) Quantifying influenza virus
#' diversity and transmission in humans. \emph{Nature Genetics}
#' \bold{48}(2), 195--200.\cr
#' Sobel Leonard, A., Weissman, D. B., Greenbaum, B., Ghedin, E.,
#' & Koelle, K. (2017) Transmission bottleneck size estimation from
#' pathogen deep-sequencing data, with an application to human
#' influenza A virus. \emph{Journal of Virology} \bold{91}(14),
#' e00171-17.\cr
#' 
#' @export
Bottleneck_size_Calculation <- function(transmission_ob,method="Beta_binomial_Approximate",plot=FALSE,show_table=FALSE,transmission_pairs=NULL,donor_depth_threshold=10, recipient_depth_threshold=10,error_filtering=0.001,log=FALSE,variant_calling=0.03,Nbmin=1,Nbmax=1000,NonSyn_or_Syn="All"){
  method_list=c("KL","Presence-Absence","Binomial","Beta_binomial_Approximate","Beta_binomial_Exact","Wright-Fisher")
  if(method%in%method_list==FALSE){stop("Please choose valid methods (KL,Presence-Absence,Binomial,Beta_binomial_Approximate,Beta_binomial_Exact,Wright-Fisher)!")}
  if(is.null(transmission_pairs)){
    ob=transmission_ob
    }
  else{
    ob=get_from_in_object(transmission_ob = transmission_ob,transmission_pairs = transmission_pairs)
  }
  id_number=paste0(donor_depth_threshold,recipient_depth_threshold,error_filtering)
  log_dir=paste0(substitute(transmission_ob),"_",method,id_number,"_log")
  log_table=log_main_function(transmission_ob,donor_threshold = donor_depth_threshold,recipient_threshold = recipient_depth_threshold,NonSyn_or_Syn=NonSyn_or_Syn,error_calling=error_filtering,variant_calling=variant_calling)
  if(log==TRUE && method != "Wright-Fisher"){
    check_file(paste0(log_dir,".csv"))
    write.csv(log_table,paste0(log_dir,".csv"))
    setwd('..')
  }
  if(method != "Wright-Fisher"){
    N_bs=lapply(ob, one_transmission_pair_process,method=method,donor_depth_threshold=donor_depth_threshold, recipient_depth_threshold=recipient_depth_threshold,error_calling=error_filtering,plot=plot,log=log,log_dir_name=log_dir,variant_calling=variant_calling,Nbmin=Nbmin,Nbmax=Nbmax,NonSyn_or_Syn=NonSyn_or_Syn)
    transmission_info=c("transmission pairs","transmission_bottleneck_size")
    N_bs=do.call(rbind.data.frame,N_bs)
    N_table=cbind.data.frame(get_transmission_pairs(ob),N_bs)
    ids=t(data.frame(strsplit(N_table[,1],"-")))
    N_table=cbind.data.frame(ids[,1],ids[,2],N_table[,2:ncol(N_table)])
    row.names(N_table)=1:nrow(N_table)
    if(ncol(N_table)==5){
      names(N_table)=c("donor","recipient","transmission_bottleneck_size","CI_low","CI_high")
    }
    row.names(N_table)=c(1:nrow(N_table))
    if(show_table==TRUE){
      name=paste0(substitute(transmission_ob),"_",method,".csv")
      check_file(name)
      write.csv(N_table,name)
    }
    return(N_table)
  }
  if(method == "Wright-Fisher"){
    N_table=Collect_values_in_table_for_wf(transmission_ob,donor_depth_threshold = donor_depth_threshold,recipient_depth_threshold = recipient_depth_threshold,error_calling = error_filtering,NonSyn_or_Syn=NonSyn_or_Syn)
    names(N_table)=c("segments","position","P_mean","Q_mean","var","Nb","P_add_eps","Q_minus_eps","var_add_eps","Nb_add_eps","P_minus_eps","Q_add_eps","var_minus_eps","Nb_minus_eps")
    if(show_table==TRUE){
      name=paste0(substitute(transmission_ob),"_",method,".csv")
      check_file(name)
      write.csv(N_table,name)
    }
    return(N_table)
  }
}




