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

brl_dt[,`:=`(POP_adx = "BRL", POP_s1 = "AFR", POP_s2="EUR")]
asw_dt[,`:=`(POP_adx = "ASW", POP_s1 = "YRI", POP_s2="CEU")]

cols <- intersect(names(asw_dt), names(brl_dt))
div_dt <- rbindlist(list(asw_dt[,..cols], brl_dt[,..cols]))





ratio_summary <- rbindlist(list(
  asw_dt[, .(POP = 'AWS', STAT = 'S', MEAN = mean(S_adx/S_s1), t(quantile(S_adx/S_s1, c(0.025, 0.975))))],
  asw_dt[, .(POP = 'AWS', STAT = 'HTZ', MEAN = mean(HTZn_adx/HTZn_s1), t(quantile(HTZn_adx/HTZn_s1, c(0.025, 0.975))))],
  brl_dt[, .(POP = 'BRL', STAT = 'S', MEAN = mean(S_adx/S_s1), t(quantile(S_adx/S_s1, c(0.025, 0.975))))],
  brl_dt[, .(POP = 'BRL', STAT = 'HTZ', MEAN = mean(HTZn_adx/HTZn_s1), t(quantile(HTZn_adx/HTZn_s1, c(0.025, 0.975))))]
))

ratio_summary

serial_fst_plots(asw_dt, c('YRI', 'CEU', 'ASW'))
serial_fst_plots(brl_dt, c('AFR', 'EUR', 'BRL'))


div_dt <- 
  rbindlist(list(  
       brl_dt[, .(stat_key = 'Nuc. diversity ratio', POP_adx = 'BRL', POP_s1 = 'AFR', POP_s2 = 'EUR', FST, value = HTZn_adx/ HTZn_s1)],
       brl_dt[, .(stat_key = 'Num. seg. sites ratio', POP_adx = 'BRL', POP_s1 = 'AFR', POP_s2 = 'EUR', FST, value = S_adx/ S_s1)],
       asw_dt[, .(stat_key = 'Nuc. diversity ratio', POP_adx = 'ASW', POP_s1 = 'YRI', POP_s2 = 'CEU', FST, value = HTZn_adx/ HTZn_s1)],
       asw_dt[, .(stat_key = 'Num. seg. sites ratio', POP_adx = 'ASW', POP_s1 = 'YRI', POP_s2 = 'CEU', FST, value = S_adx/ S_s1)]
  ))

div_dt <- div_dt[value < 2]
div_dt[,POP_ratio := paste0(POP_adx, " / ", POP_s1)]
plot_keys <- div_dt[,.GRP, keyby = .(stat_key, POP_ratio)]
plot_keys[,`:=`(label = paste0("(", letters[1:4], ")"), x = rep(div_dt[,min(FST)],4), y = rep(div_dt[,max(value)],4))]

  {
ggplot(div_dt[value < 2], aes(y=value, x=FST)) +
    geom_hex(aes(fill = log10(..count..))) + 
    theme_pubr(base_size=16) + 
    xlab(TeX("$F_{st}$ between ancestrals")) + 
    ylab(NULL) +
    labs(fill = TeX("$\\log_{10}$ ($N^{\\underline{0}}$  windows)")) +
    theme(legend.position = "bottom",
          legend.title = element_text(size=12),
          legend.text = element_text(size=10),
          strip.text = element_text(size=16),
          strip.background = element_blank(),
          strip.placement = "outside") +
    geom_hex(aes(fill = log10(..count..))) +
    geom_text(data = plot_keys, aes(x=x,y=y, label=label), size = 6)+
    scale_fill_viridis_c(direction = -1, breaks=pretty_breaks()) +
    geom_hline(yintercept = 1, color = 'black') +
    geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = 'red') + 
    stat_poly_eq(formula = y ~ x, parse = TRUE, size = 5, coef.digits=3, label.x = 'right', label.y = 'top') +
    facet_grid(stat_key~POP_ratio, 
               scales = "free_y",
               switch = "y"
    )
    ggsave("../figures/stat_ratio_fst.pdf")
  }

