library(data.table)
library(ggplot2)
library(hexbin)
library(ggpmisc)
library(latex2exp)
library(ggpubr)
library(scales)

plot_wd_stats <- 
  function(data_wd, xformula, yformula, xlabel, ylabel, gtitle, yl, yintercept = 1, figname = NULL){
    fig <- 
      ggplot(data_wd, aes_string(x = xformula, y = yformula)) + 
        theme_pubr(base_size=18) + 
        xlab(xlabel) + 
        ylab(NULL) +
        labs(title = gtitle, subtitle = ylabel) + 
        labs(fill = TeX("$\\log_{10}$ ($N^{\\underline{0}}$  windows)")) +
        theme(legend.position = "bottom") +
        geom_hex(aes(fill = log10(..count..))) +
        scale_fill_viridis_c(direction = -1, breaks=pretty_breaks()) +
        geom_hline(yintercept = yintercept, color = 'black') +
        geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'red') + 
        ylim(yl) + 
        stat_poly_eq(formula = y ~ x, parse = TRUE, size = 6) 
      if(!is.null(figname)){
        ggsave(paste0('../figures/', figname))
      }
      fig
  }

serial_fst_plots <-
  function(data_wd, pop_names){
    stats <- c('HTZn_', 'S_') 
    stat_labels <- c('\\pi', 'S')
    stat_title <- c('Nucleotide diversity ratio', 'Number of segragating sites ratio')
    for(i in 1:2){
      for(j in 1:2){
        sl <- stat_labels[j]
        st <- stat_title[j]
        sv <- stats[j]
        xformula <- 'FST'
        yformula <- paste0(sv, 'adx / ', sv, 's', i)
        xlabel <- TeX(paste0("$F_{st}$(",pop_names[1], ", ", pop_names[2], ")"))
        #ylabel <- TeX(paste0("$", sl, "(", pop_names[3], ")","/", sl, "(", pop_names[i], ")","$"))
        ylabel <- paste0(pop_names[3], " / ", pop_names[i])
        fname <- paste0(sv, pop_names[3], '-', pop_names[i], '.png') 
        yl <- data_wd[,quantile(eval(parse(text=yformula)),c(0.001, 0.999))]
        plot_wd_stats(data_wd = data_wd, 
                      xformula = xformula,
                      yformula = yformula,
                      xlabel = xlabel,
                      ylabel = ylabel,
                      gtitle = st,
                      yl = yl,
                      #figname = NULL )
                      figname = fname)
      }
    }
  }

asw_dt <- fread('../data/results/asw_results.csv')
brl_dt <- fread('../data/results/brl_results.csv')

serial_fst_plots(asw_dt, c('YRI', 'CEU', 'ASW'))
serial_fst_plots(brl_dt, c('AFR', 'EUR', 'BRL'))
