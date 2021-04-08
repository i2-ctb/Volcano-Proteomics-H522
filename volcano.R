################################################################################
# Install and load libraries
################################################################################

# Uncomment the lines below to install any missing packages
#install.packages("tidyverse")
#install.packages("ggthemes")
#install.packages("ggrepel")
#install.packages("here")

# Load
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(here)

################################################################################
# Read in and process our data
################################################################################

data <- read_tsv(here("example data", "proteomics.tsv"))


################################################################################
# Make the plot
################################################################################

proteomics.log.fc.threshold <- log2(2)
proteomics.log.stats.threshold <- -log10(0.05)

data$logFC <- -data$logFC
min.x <- min(data$logFC)
max.x <- max(data$logFC)
x.range <- max.x - min.x
x.limit <- max(abs(min.x), abs(max.x))

data$log.adj.P.Val <- -log10(data$adj.P.Val)
y.range <- max(data$log.adj.P.Val) - min(data$log.adj.P.Val)

data <- data %>% mutate(significance = case_when(SARS_CoV_2 == T ~ "SARS-CoV-2",
                                                 abs(logFC) >= proteomics.log.fc.threshold & log.adj.P.Val >= proteomics.log.stats.threshold ~ "Sig. FC and q-value",
                                                 #abs(logFC) >= proteomics.log.fc.threshold ~ "Sig. fold-change",
                                                 #log.adj.P.Val >= proteomics.log.stats.threshold ~ "Sig. statistic",
                                                 TRUE ~ "Not significant"),
                        Label = case_when(SARS_CoV_2 == T ~ str_replace(`Gene names`, "SARS_CoV_2", ""),
                                          TRUE ~ ""))

################################################################################
# Volcano plot 96 hrs
################################################################################
COV2.color <- "#fb8072"

panel <- (ggplot(data=data, aes(x=logFC, y=log.adj.P.Val, label=Label, color=significance))
        + geom_point(size=0.5)

        + geom_vline(xintercept = -proteomics.log.fc.threshold, linetype="dashed", color=grey(0.7))
        + geom_vline(xintercept = proteomics.log.fc.threshold, linetype="dashed", color=grey(0.7))
        + geom_hline(yintercept = proteomics.log.stats.threshold, linetype="dashed", color=grey(0.7))

        + geom_text_repel(size=2, color=COV2.color)

        + scale_color_manual(name="Significance",
                             values=c("Not significant"="#AAAAAA55",
                                      "Sig. fold-change"="#AAAAAA55",
                                      "Sig. statistic"="#AAAAAA55",
                                      "Sig. FC and q-value"="#1f78b455",
                                      "SARS-CoV-2"="#fb8072"))

        + scale_x_continuous(limits=c(min.x - x.range/20, x.limit + x.range/20),
                             breaks = seq(-3,6,1))

        + ylab("-log10(q-value)")
        + xlab("log2(fold-change)")

        + ggtitle("96 hpi vs 96 hr mock",
                  subtitle = str_c("n =",
                                   prettyNum(nrow(data), big.mark=","),
                                   "proteins",
                                   sep=" "))

        + theme_clean(base_size = 7)
        + theme(panel.grid.major.y = element_blank(),
                legend.background = element_rect(fill=grey(0.95), color=NA),
                legend.title = element_text(size=6),
                legend.text = element_text(size=6),
                legend.key = element_blank(),
                legend.key.size = unit(3.0, 'mm'),
                legend.position = c(0.85, 0.13),
                legend.margin = margin(3,3,3,3),
                plot.title = element_text(hjust=0.5),
                plot.subtitle = element_text(hjust=0.5, face="italic"))
        )

print(panel)

ggsave(here("Volcano.pdf"), width=3.5, height=3, compress=F, panel)
