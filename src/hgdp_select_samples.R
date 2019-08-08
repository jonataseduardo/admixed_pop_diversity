library(data.table)


hgdp_samples <- fread('/scratch/genevol/hgdp3.0/hgdp_wgs.20190516.metadata.txt')

regions <- hgdp_samples[,unique(region)]
populations <- hgdp_samples[,unique(population)]

for(r in regions){
  fwrite(hgdp_samples[region == r, .(sample)], paste0('/scratch/genevol/hgdp3.0/', r, '.txt'), col.names = FALSE)
}

for(p in populations){
  fwrite(hgdp_samples[population == p, .(sample)], paste0('/scratch/genevol/hgdp3.0/', p, '.txt'), col.names = FALSE)
}
