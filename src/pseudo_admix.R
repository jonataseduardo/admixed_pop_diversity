library(data.table)
library(diverdt)
library(ggplot2)
library(devtools)
library(pkgload)
library(hexbin)
library(ggpmisc)

devtools::document('~/diverdt')
devtools::load_all('~/diverdt')


merge_pop_data <-
  function(pop_freqs){
    fst_dt <- hudson_fst(pop_freqs[1:2])

    on_cols = c("CHROM", "POS", "ID", "REF", "ALT")
    scols = c("AF", "HTZ", "OBS_CT", "POP", on_cols)

    pop_dt1 <- merge(pop_freqs[[1]][, ..scols], 
                    pop_freqs[[2]][, ..scols], 
                    by = on_cols,
                    all = TRUE)
    pop_dt2 <- merge(pop_dt1, fst_dt, by=on_cols, all=TRUE)
    pop_dt3 <- merge(pop_dt2, pop_freqs[[3]][,..scols], by=on_cols, all=TRUE)
    return(pop_dt3)
  }

plot_wd_stats <- 
  function(data_wd, xformula, yformula, xlabel, ylabel, yintercept = 1, figname = 'teste.png'){
    ggplot(data_wd, aes_string(x = xformula, y = yformula)) + 
      theme_bw() + 
      xlab(xlabel) + 
      ylab(ylabel) +
      labs(fill = expression(log[10]*'(n windows)')) +
      theme(legend.position = "top") +
      geom_hex(aes(fill = log10(..count..))) +
      scale_fill_viridis_c(direction = -1) +
      geom_hline(yintercept = yintercept, color = 'black') +
      geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'red') + 
      stat_poly_eq(formula = y ~ x, parse = TRUE) 
    ggsave(paste0('../figures/', figname))
  }

pops <- c('YRI', 'CEU', 'ASW')
l_chr <- 10:16

pop_freqs <- 
    lapply(
      pops, 
      function(pop_name){
        rbindlist(
          lapply(l_chr, 
            function(chr){
              fname <- paste0('../data/', pop_name, '61/', pop_name, '_chr', chr)
              freq_data <- load_pvar_afreq(fname, pop_name = pop_name, read_info = FALSE)
              freq_data[, HTZ := 2 * AF * (1 - AF)]
              return(freq_data[FILTER == 'PASS'])
            })
        )
      }
    )

 
tg_data <- merge_pop_data(pop_freqs)

width = 2e4
S_wd <- 
  tg_data[, .(S_adx= length(.I[which(!is.na(AF))]), 
              S_s1= length(.I[which(!is.na(AF.x))]), 
              S_s2= length(.I[which(!is.na(AF.y))]), 
              HTZ_adx = mean(HTZ, na.rm = TRUE), 
              HTZ_s1 = mean(HTZ.x, na.rm = TRUE), 
              HTZ_s2 = mean(HTZ.y, na.rm = TRUE), 
              HTZn_adx = sum(HTZ, na.rm = TRUE) / width, 
              HTZn_s1 = sum(HTZ.x, na.rm = TRUE) / width, 
              HTZn_s2 = sum(HTZ.y, na.rm = TRUE) / width, 
              FST = sum(T1, na.rm = TRUE)/sum(T2, na.rm = TRUE) 
             )
          ,by = .(CHROM,floor(POS / width))]

min_S <- 15
data_wd <- S_wd[!is.na(FST) & S_adx > min_S & S_s1 > min_S & S_s1 > min_S]

plot_wd_stats(data_wd = data_wd, 
              xformula = 'FST',
              yformula = 'HTZn_adx / HTZn_s1',
              xlabel = expression('F'['st']*'(CEU,YRI)'), 
              ylabel = expression(pi['ASW']/pi['YRI']), 
              figname = 'teste.png')
