
margin_value<-0.1

library(tidyverse)

threshold_mat <- read.csv("threshold_mat.csv")
results_mat <- read.csv("results_mat.csv")



table(apply(results_mat[results_mat$true_es<results_mat$margin,c("equiv", "hdi", "bf", "tost")],1, which.max))
table(apply(results_mat[results_mat$true_es>results_mat$margin,c("equiv", "hdi", "bf", "tost")],1, which.min))





pr_all <- results_mat %>% gather(key = "method",
                                 value = "pr_equiv",
                                 equiv, hdi, bf, tost) %>%
  mutate(method = factor(x = .$method),
         scale = scale) %>%
  distinct(margin,
           n,
           true_es,
           scale,
           method,
           .keep_all = TRUE) %>%
  mutate(sr = as.character(round(x = .$scale,
                                 digits = 1))) %>%
  mutate(id = ifelse(test = method %in% c("equiv", "tost"),
                     yes = as.character(method),
                     no = paste(.$method, .$sr, sep = "_"))) %>%
  mutate(id = factor(
    x = id,
    levels = c("equiv",
               "hdi_0.4", "hdi_0.7", "hdi_1.4",               
               "bf_0.4", "bf_0.7", "bf_1.4","tost")
  ))


pr_all<- pr_all[pr_all[,"sr"]%in%c("0.4", "0.7", "1.4"),]
###  We will only look at one-third of the frequentist 
###  results since we have three results (one for every 
###  scale parameter value)

freq_ids<- apply(pr_all[pr_all[,"method"]%in%c("TOST", "equiv"), 
                        !colnames(pr_all)%in%c("jjj", "sr", "scale", "pr_equiv")],1,
                 function(x) paste(x, collapse = ""))

pr_all_freq <- pr_all[match(unique(freq_ids), freq_ids),]

dim(pr_all[pr_all[,"method"]%in%c("TOST", "equiv"), ])

# removing all freq methods because there are triplicates (due to scale)
pr_all <- pr_all[!pr_all[,"method"]%in%c("TOST", "equiv"), ]
dim(pr_all)
# 
pr_all <- rbind(pr_all,pr_all_freq)
dim(pr_all)
head(pr_all)
tail(pr_all)
####






### making the plot ###


plot_int03 <- ggplot(data = pr_all[pr_all$margin == margin_value, ],
                     mapping = aes(x = true_es)) +
  geom_hline(yintercept = 0.05,
             size = 0.5,
             colour = "grey80",
             linetype = "dashed") +
  
  geom_segment(mapping=aes(x=margin_value, xend=margin_value, y=0, yend=0.7),
               size = 0.5,
               colour = "grey80",
               linetype = "dashed") +
  geom_line(aes(y = pr_equiv, 
                colour = id,
                linetype = id),
            lwd = 0.75) +
  facet_wrap(facets = ~n,
             labeller = label_bquote(italic("n") == .(n))) +
  ggtitle(label = expression(italic("m") == 0.1))  +
  scale_x_continuous(name = expression(delta),
                     breaks = seq(from = 0,
                                  to = 0.5,
                                  by = 0.1),
                     labels = c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                     limits = c(0, 0.5)) +
  scale_y_continuous(name = "Proportion predict equivalence",
                     breaks = seq(from = 0,
                                  to = 1,
                                  by = 0.2),
                     labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
                     limits = c(0, 1)) +
  scale_colour_manual(name = NULL,
                      values = c(
                        "purple",
                        "blue",
                        "blue",
                        "blue",
                        "green",
                        "green",
                        "green", "darkred"
                      ),
                      labels = c(
                        "Freq. `optimal test`",
                        bquote("HDI(" ~Z1~ "), " ~ italic("r") == 0.5 / sqrt(2)),
                        bquote("HDI(" ~Z2~ "), " ~ italic("r") == 1 / sqrt(2)),
                        bquote("HDI(" ~Z3~ "), " ~ italic("r") == 2 / sqrt(2)),
                        bquote("BF > " ~ B1 ~ ", " ~ italic("r") == 0.5 / sqrt(2)),
                        bquote("BF > " ~ B2 ~ ", " ~ italic("r") == 1 / sqrt(2)),
                        bquote("BF > " ~ B3 ~ ", " ~ italic("r") == 2 / sqrt(2)),
                        "TOST"
                      )) +
  scale_linetype_manual(name = NULL,
                        values = c(
                          "solid",
                          "dashed",
                          "solid",
                          "dotted",
                          "dashed",
                          "solid",
                          "dotted", "solid"
                        ),
                        labels =  c(
                          "Freq. `optimal test`",
                          bquote("HDI(" ~Z1~ "), " ~ italic("r") == 0.5 / sqrt(2)),
                          bquote("HDI(" ~Z2~ "), " ~ italic("r") == 1 / sqrt(2)),
                          bquote("HDI(" ~Z3~ "), " ~ italic("r") == 2 / sqrt(2)),
                          bquote("BF > " ~ B1 ~ ", " ~ italic("r") == 0.5 / sqrt(2)),
                          bquote("BF > " ~ B2 ~ ", " ~ italic("r") == 1 / sqrt(2)),
                          bquote("BF > " ~ B3 ~ ", " ~ italic("r") == 2 / sqrt(2)),
                          "TOST"
                        )) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12),
        strip.background = element_rect(fill = "white",
                                        colour = "white"),
        plot.title = element_text(hjust = 0.5),
        panel.spacing = unit(x = 0.25, units = "in"),
        legend.key.height = unit(x = 0.35, units = "in"),
        legend.key.width = unit(x = 0.35, units = "in"),
        legend.text.align = 0,
        legend.margin = margin(t = 0,
                               r = 10,
                               b = 0,
                               l = 0),
        legend.box.margin = margin(t = 0,
                                   r = 0,
                                   b = 0,
                                   l = -5),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.title.x = element_text(margin = margin(t = 5)))  
plot_int03

BF_thres <- round(t(
  matrix(  threshold_mat[threshold_mat[,"margin"]==margin_value, "BFthres_sim"],4,)))

dat_text <- data.frame(
  label = c(
    paste0(c("B1 = ", "B2 = ","B3 = "), BF_thres[,1], collapse=", "), 
    paste0(c("B1 = ", "B2 = ","B3 = "), BF_thres[,2], collapse=", "),
    paste0(c("B1 = ", "B2 = ","B3 = "), BF_thres[,3], collapse=", "), 
    paste0(c("B1 = ", "B2 = ","B3 = "), BF_thres[,4], collapse=", ")),
  n   = c(50, 100, 250, 500)
)
HDI_thres <- round(t(
  matrix(  threshold_mat[threshold_mat[,"margin"]==margin_value, "HDIprob_sim"],4,)),2)


hdi_text <- data.frame(
  label = c(paste0(c("Z1 = ", "Z2 = ","Z3 = "), HDI_thres[,1], collapse=", "), 
            paste0(c("Z1 = ", "Z2 = ","Z3 = "), HDI_thres[,2], collapse=", "),
            paste0(c("Z1 = ", "Z2 = ","Z3 = "), HDI_thres[,3], collapse=", "), 
            paste0(c("Z1 = ", "Z2 = ","Z3 = "), HDI_thres[,4], collapse=", ")),
  n   = c(50, 100, 250, 500)
)

plot_int03 +
  geom_text(
    data    = dat_text,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.05,
    vjust   = 2, cex=2.75)+
  
  geom_text(
    data    = hdi_text,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1.05,
    vjust   = 4, cex=2.75)

