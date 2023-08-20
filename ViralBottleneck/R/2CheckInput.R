#'@importFrom methods setRefClass
#'@importFrom utils read.table
########################################################################################################################
#########################################################################################################################
#make sure the column don't contain the repeated info
check_not_repeated <- function(table_column){
  list = data.frame(table(table_column))[,2]
  repeated = list[list>1]
  if (length(repeated)==0){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#produce repeated message
repeated_message <- function(table_column){
  message=data.frame(table(table_column))
  for(i in 1:nrow(message)){
    row = message[i,]
    if(row$Freq > 1){
      message(paste(row$table_column,"is deplicated",row$Freq,"times."))
    }
  }
}

#########################################################################################################################
#create transmission pairs id
get_trans_ID <- function(row){
  ID=paste0(row[1],"-",row[2])
  return(ID)
}

produce_transmission_pairs <- function(transmission_pairs_table){

  transmission_id = apply(transmission_pairs_table,1, get_trans_ID)
  transmission_pairs_table$transmission_ID = transmission_id
  return(transmission_pairs_table)
}

#########################################################################################################################
#create variant id
get_variant_ID <- function(row){
  ID=paste0(row[2],"-",row[1])
  return(ID)
}

produce_variant_id <- function(variant_site_table){
  variant_id = apply(variant_site_table,1, get_variant_ID)
  variant_site_table$variant_ID = variant_id
  return(variant_site_table)
}

#########################################################################################################################
#make sure the coloumn do not have missing value
check_not_missing <- function(table_coloumn){
  if (NA %in% table_coloumn){
    return(FALSE)
  }
  else {
    return(TRUE)
  }
}

#tell the users the info of missing value
Print_missing_inTable <- function(table,table_name){
  for (i in 1:ncol(table)){
    if (check_not_missing(table[,i])==FALSE){message(paste("There are missing values in column", i,"in",table_name))}
  }
}

check_not_missing_table <- function(table){
  if(FALSE %in% lapply(t(table),check_not_missing)){
    return(FALSE)
    stop(paste(table,"has missing value"))
  }
  else{return(TRUE)}
}

#########################################################################################################################
#make sure the number of columns are correct
Check_table_dim <-function(table, valided_number_column){
  if(ncol(table)==valided_number_column){return(TRUE)}
  else{
    return(FALSE)
    message(paste("The number of column is not", valided_number_column, "in", table))
  }
}

#make sure the transmission pairs not duplicated
Check_transmission_pairs_not_duplicated <- function(transmission_pairs_table){
  if(check_not_repeated(transmission_pairs_table[,3])){
    return(TRUE)
  }
  else{
    return(FALSE)
    message("There are transmission pairs duplicated:")
    repeated_message(transmission_pairs_table$transmission_ID)
  }
}

#make sure the variant not duplicated
Check_variant_not_duplicated <- function(variant_site_table){
  if(check_not_repeated(variant_site_table[,8])){
    return(TRUE)
  }
  else{
    return(FALSE)
    message("There are variants duplicated:")
    repeated_message(variant_site_table$variant_ID)
  }
}

#make sure each column containing the valid info
Check_column_type <- function(column){
  if(is.integer(column)){return(TRUE)}
  else{
    message(paste(column,"is not valid, it is not integer, please check"))
    return(FALSE)
  }
}

Check_each_column_variant_sites_table <- function(variant_site_table){
  if(FALSE %in% lapply(t(variant_site_table[,3:6]),Check_column_type)){
    message(paste("The frequency part in",variant_site_table,"is not filled by integers,please check"))
    return(FALSE)
  }
  else{return(TRUE)}
}

#########################################################################################################################
#convert the name to the file name
Produce_sample_file_name <- function(sample_name){
  file_name=paste0(sample_name,".csv")
  return(file_name)
}

#convert a name list to a file name list in order to compare with the file names in working directory.
Produce_sample_file_name_list<- function(sample_name_list){
  file_name_list = lapply(sample_name_list, Produce_sample_file_name)
  return(file_name_list)
}

#find the file which is not exist in directory but exist in transmission pairs
Print_difference_message<-function(element){
  message(paste("Can not find",element,"please check"))
}

Print_difference_message_inlist <- function(difference_list){
  lapply(difference_list, Print_difference_message)
}
#make sure the sample file exsit in working directory
Check_sample_file_exist <- function(transmission_pairs_table){
  files_in_working_directory=list.files(full.names = FALSE, recursive = FALSE)
  transmission_samples=union(transmission_pairs_table[,1],transmission_pairs_table[,2])
  transmission_samples = Produce_sample_file_name_list(transmission_samples)
  check=intersect(transmission_samples,files_in_working_directory)
  if(setequal(check,transmission_samples)){
    return(TRUE)
  }
  else{
    difference=transmission_samples[!(transmission_samples %in% files_in_working_directory)]
    Print_difference_message_inlist(difference)
    stop("Some samples files are missing, please check directory!")
    return(FALSE)
  }
}

#########################################################################################################################
#make sure the transmission pairs table is valid
Check_transmission_pair_table <- function(transmission_pairs_table){
  if(is.data.frame(transmission_pairs_table)){
    if(Check_table_dim(transmission_pairs_table,2)){
      if(check_not_missing_table(transmission_pairs_table)){
        transmission_pairs_table=produce_transmission_pairs(transmission_pairs_table)
        if(Check_transmission_pairs_not_duplicated(transmission_pairs_table)){
          if(check_not_repeated(transmission_pairs_table[,1])&&check_not_repeated(transmission_pairs_table[,2])){
            return(TRUE)
          }
          else{
            message("In donors:")
            repeated_message(transmission_pairs_table[,1])
            message("In recipients:")
            repeated_message(transmission_pairs_table[,2])
            return(TRUE)
          }
        }
        else{return(FALSE)}
      }
      else{return(FALSE)}
    }
    else{return(FALSE)}}
  else{return(FALSE)}
}

#########################################################################################################################
#make sure the variant site table is valid
Check_variant_site_table <-function(variant_site_table){
  name_table=variant_site_table
  variant_site_table=read.table(paste0(name_table,".csv"),header = TRUE,sep = ",")
  if(Check_table_dim(variant_site_table,7)){
    Print_missing_inTable(variant_site_table,name_table)
    if(check_not_missing_table(variant_site_table)){
      if(Check_each_column_variant_sites_table(variant_site_table)){
        variant_site_table=produce_variant_id(variant_site_table)
        if(Check_variant_not_duplicated(variant_site_table)){
          return(TRUE)
        }
        else{return(FALSE)}
      }
      else{return(FALSE)}
    }
    else{return(FALSE)}
  }
  else{return(FALSE)}
}


Check_variant_tables_list<-function(sample_list){
  if(FALSE %in% lapply(sample_list, Check_variant_site_table)){
    return(FALSE)
    message("There are some invalided inputs, please check!")
  }
  else{return(TRUE)}
}

#########################################################################################################################
#check every thing is valid before creating transmission object
Check_everything <- function(transmission_pair_table){
  if(Check_transmission_pair_table(transmission_pair_table)){
    if(Check_sample_file_exist(transmission_pair_table)){
      sample_list=union(transmission_pair_table[,1],transmission_pair_table[,2])
      if(Check_variant_tables_list(sample_list)){
        return(TRUE)
      }
      else{return(FALSE)}
    }
    else{return(FALSE)}
  }
  else{
    return(FALSE)
  }
}


