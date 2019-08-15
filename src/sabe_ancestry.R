library(data.table)

sabe_ancestry <- fread('/scratch/genevol/sabe/SABE_Global_ancestry.txt')

s_sabe <- sabe_ancestry[NAM < 0.05 & AFR > 0.1 & AFR < 0.35]

sabe_ancestry[,.(mean(NAM), mean(AFR), mean(EUR))]

s_sabe[,.(.N, mean(NAM), mean(AFR), mean(EUR))]

#file_name <- 
#  s_sabe[,  paste0('sample_list_sabe_N_', .N, 
#            '_NAM_', round(mean(NAM), digits = 2), 
#            '_AFR_' , round(mean(AFR), digits = 2),
#            '_EUR_' , round(mean(EUR), digits = 2),
#            '.txt')]

file_path <- s_sabe[,  paste0('../data/BRL', .N, '/')]
file_name <- 'BRL_samples.txt' 

system(paste0('mkdir -p ', file_path))
fwrite(s_sabe[,.(IND)] , paste0(file_path, file_name), col.name = FALSE) 
