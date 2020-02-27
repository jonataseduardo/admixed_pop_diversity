library(data.table)
library(ggplot2)
library(hexbin)
library(ggpmisc)


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

serial_fst_plots <-
  function(data_wd, pop_names){

    stats <- c('HTZn_', 'S_') 
    stat_labels <- c('\\pi', 'S')
    for(i in 1:2){
      for(j in 1:2){
        sl <- stat_labels[j]
        sv <- stats[j]
        xformula <- 'FST'
        yformula <- paste0(sv, 'adx / ', sv, 's', i)
        xlabel <- TeX(paste0("F_{st}(",pop_names[1], ",", pop_names[2], ")"))
        ylabel <- TeX(paste0("$", sl, "_{", pop_names[3], "}","/", sl, "_{", pop_names[i], "}","$"))
        plot_wd_stats(data_wd = data_wd, 
                      xformula = xformula,
                      yformula = yformula,
                      xlabel = xlabel,
                      ylabel = ylabel,
                      figname = fname)
      }
    }
  }


asw_dt <- fread('../data/results/asw_results.csv')
brl_dt <- fread('../data/results/brl_results.csv')

serial_fst_plots(asw_dt, c('YRI', 'CEU', 'ASW'))
serial_fst_plots(brl_dt, c('AFR', 'EUR', 'BRL'))
