if (!requireNamespace("ViralBottleneck", quietly = TRUE)) {
  # devtools is heavier; remotes works too
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")          # only once
  devtools::install_github("BowenArchaman/ViralBottleneck",
                           dependencies = TRUE,   # pull CRAN deps
                           upgrade = "never")     # don’t re-upgrade installed deps
}
library(ViralBottleneck)
setwd("/home3/2748683z/R_pub/First_simulation/BottleneckRDataset_noerror/R_input")
transmission_pairs = read.table("BottleneckR_transmission_pairs.csv",header = TRUE,sep = ",")
H1N1_ob = ViralBottleneck::CreateTransmissionObject(transmission_pairs)
KL_output = ViralBottleneck::Bottleneck_size_Calculation(H1N1_ob, method = "KL",variant_calling = 0.01,Nbmin=1,Nbmax=400,donor_depth_threshold =0,recipient_depth_threshold = 0,error_filtering = 0,show_table =TRUE)
PA_output= ViralBottleneck::Bottleneck_size_Calculation(H1N1_ob, method = "Presence-Absence",variant_calling = 0.01,Nbmin=1,Nbmax=400,donor_depth_threshold =0,recipient_depth_threshold = 0,error_filtering = 0,show_table =TRUE)
App=ViralBottleneck::Bottleneck_size_Calculation(H1N1_ob, method = "Beta_binomial_Approximate",variant_calling = 0.01,Nbmin=1,Nbmax=400,donor_depth_threshold =0,recipient_depth_threshold = 0,error_filtering = 0,show_table =TRUE)
Exact=ViralBottleneck::Bottleneck_size_Calculation(H1N1_ob, method = "Beta_binomial_Exact",variant_calling = 0.01,Nbmin=1,Nbmax=400,donor_depth_threshold =0,recipient_depth_threshold = 0,error_filtering = 0,show_table =TRUE)
binomial=ViralBottleneck::Bottleneck_size_Calculation(H1N1_ob, method = "Binomial",variant_calling = 0.01,Nbmin=1,Nbmax=400,donor_depth_threshold =0,recipient_depth_threshold = 0,error_filtering = 0,show_table =TRUE)
