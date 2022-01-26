#############################################################################
## CA OAT CEA Model
## CE Plane - PSA Results
## Feb 27, 2017
#############################################################################
#library(XLConnect)
library(ggplot2)
#library(xlsx)
library(readxl)
library(grid)
library(gtable)
#library(hexbin)
#library(SIBER)

setwd(dirname("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE"))

#Input data for total costs and QALYs CE plot
input.totals <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\CA_OAT_CEA_totals.xlsx", sheet = "Sheet1", col_names = TRUE)
#Subset for societal and health perspective
total.societal <- subset(input.totals, select = c(QALY, total.costs, intervention))
total.health <- subset(input.totals, select = c(QALY, health.total, intervention))

#Input data
input.results <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\PSA_SecondDraft_RESULTS-v4_RR_BE.xlsx", sheet = "EBvOP CEA Results", col_names = TRUE)
input.results2 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\PSA_SecondDraft_RESULTS-v4_RR_BE.xlsx", sheet = "EBvCG CEA Results", col_names = TRUE)
input.results3 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\PSA_SecondDraft_RESULTS-v4_RR_BE.xlsx", sheet = "OPvCG CEA Results", col_names = TRUE)

#PSA - OWSA
input.results4 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\R&R PSA OWSA\\PSA_OWSA_HRU_REL_LO_2000nSim_RR_BE.xlsx", sheet = "CEA Results", col_names = TRUE)
input.results5 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\R&R PSA OWSA\\PSA_OWSA_Crime_REL_LO_2000nSim_RR_BE.xlsx", sheet = "CEA Results", col_names = TRUE)
input.results6 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\R&R PSA OWSA\\PSA_OWSA_REL_ABS_2000nSim_RR_BE.xlsx", sheet = "CEA Results", col_names = TRUE)
input.results7 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\R&R PSA OWSA\\PSA_OWSA_All_to_REL_2000nSim_RR_BE.xlsx", sheet = "CEA Results", col_names = TRUE)
input.results8 <- read_excel("Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\R&R PSA OWSA\\PSA_OWSA_DET-to-ABS_2000nSim_RR_BE.xlsx", sheet = "CEA Results", col_names = TRUE)

#Societal perspective
dat <- subset(input.results, select = c(d.QALYs, d.total.costs)) #Evidence-based vs. observed practice
dat2 <- subset(input.results2, select = c(d.QALYs, d.total.costs)) #Evidence-based vs. clinical guidelines
dat3 <- subset(input.results3, select = c(d.QALYs, d.total.costs)) #Observed practice vs. clinical guidelines

dat4 <- subset(input.results4, select = c(d.QALYs, d.total.costs)) #PSA-OWSA for HRU costs
dat5 <- subset(input.results5, select = c(d.QALYs, d.total.costs)) #PSA-OWSA for crime costs
dat6 <- subset(input.results6, select = c(d.QALYs, d.total.costs)) #PSA-OWSA for REL to ABS
dat7 <- subset(input.results7, select = c(d.QALYs, d.total.costs)) #PSA-OWSA for All to REL
dat8 <- subset(input.results8, select = c(d.QALYs, d.total.costs)) #PSA-OWSA for DET to ABS

#Health costs
health.dat <- subset(input.results, select = c(d.QALYs, d.health.total)) #Evidence-based vs. observed practice
health.dat2 <- subset(input.results2, select = c(d.QALYs, d.health.total)) #Evidence-based vs. clinical guidelines
health.dat3 <- subset(input.results3, select = c(d.QALYs, d.health.total)) #Observed practice vs. clinical guidelines

health.dat4 <- subset(input.results4, select = c(d.QALYs, d.health.total)) #PSA-OWSA for HRU costs
health.dat5 <- subset(input.results5, select = c(d.QALYs, d.health.total)) #PSA-OWSA for crime costs
health.dat6 <- subset(input.results6, select = c(d.QALYs, d.health.total)) #PSA-OWSA for REL to ABS
health.dat7 <- subset(input.results7, select = c(d.QALYs, d.health.total)) #PSA-OWSA for All to REL
health.dat8 <- subset(input.results8, select = c(d.QALYs, d.health.total)) #PSA-OWSA for DET to ABS

t1 <- ggplot(total.societal, aes(QALY, total.costs, shape=intervention, colour=intervention))
t2 <- ggplot(total.health, aes(QALY, health.total, shape=intervention, colour=intervention))

p1 <- ggplot(dat, aes(d.QALYs, d.total.costs))
p2 <- ggplot(dat2, aes(d.QALYs, d.total.costs))
p3 <- ggplot(dat3, aes(d.QALYs, d.total.costs))

p4 <- ggplot(dat4, aes(d.QALYs, d.total.costs))
p5 <- ggplot(dat5, aes(d.QALYs, d.total.costs))
p6 <- ggplot(dat6, aes(d.QALYs, d.total.costs))
p7 <- ggplot(dat7, aes(d.QALYs, d.total.costs))
p8 <- ggplot(dat8, aes(d.QALYs, d.total.costs))

health.p1 <- ggplot(health.dat, aes(d.QALYs, d.health.total))
health.p2 <- ggplot(health.dat2, aes(d.QALYs, d.health.total))
health.p3 <- ggplot(health.dat3, aes(d.QALYs, d.health.total))

health.p4 <- ggplot(health.dat4, aes(d.QALYs, d.health.total))
health.p5 <- ggplot(health.dat5, aes(d.QALYs, d.health.total))
health.p6 <- ggplot(health.dat6, aes(d.QALYs, d.health.total))
health.p7 <- ggplot(health.dat7, aes(d.QALYs, d.health.total))
health.p8 <- ggplot(health.dat8, aes(d.QALYs, d.health.total))

#Scatter Plot for all PSA points (full display)
#Total costs and QALYs (societal)
psa_plot_total_societal <- t1 +
  geom_point() +
  labs(shape="Treatment Strategy", colour="Treatment Strategy") +
  scale_x_continuous(breaks = c(10,12,14,16), labels = c("10","12","14","16"), limits = c(10,16)) +
  scale_y_continuous(breaks = c(500000, 750000, 1000000, 1250000, 1500000, 1750000), labels = c("$500,000", "$750,000", "$1,000,000", "$1,250,000", "$1,500,000", "$1,750,000"), limits = c(470000, 1800000)) +
  xlab("Total QALYs") + ylab("Total Costs (2016 USD)")

#Total costs and QALYs (health sector)
psa_plot_total_health <- t2 +
  geom_point() +
  labs(shape="Treatment Strategy", colour="Treatment Strategy") +
  scale_x_continuous(breaks = c(10,12,14,16), labels = c("10","12","14","16"), limits = c(10,16)) +
  scale_y_continuous(breaks = c(150000, 400000, 650000, 900000, 1150000, 1400000), labels = c("$150,000", "$400,000", "$650,000", "$900,000", "$1,150,000", "$1,400,000"), limits = c(180000, 1300000)) +
  xlab("Total QALYs") + ylab("Total Costs (2016 USD)")

#Incremental plots
#Societal perspective
#Evidence based v. observed practice
psa_plot_points <- p1 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point() +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 99.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 99.9%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.9%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)") 

#Evidence based v. observed practice (zoom)
psa_plot_points_zoom <- p1 + 
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_point() +
  annotate("text", x = 0.75, y = 30000, label = "ICER = $50,000") +
  #Add point at average
  #geom_point(data = mean, aes(x,y), colour = "red", size = 4) +
  #stat_ellipse(data = dat, level = 0.95) +
  #stat_ellipse(data = dat, level = 0.99) +
  scale_x_continuous(breaks = c(-0.5,0,1), labels = c("-0.5","0","1"), limits = c(-0.5,1)) +
  scale_y_continuous(breaks = c(50000, 0, -100000), labels = c("$50,000", "$0", "-$100,000"), limits = c(-100000, 50000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)") 

#Evidence based v. clinical guidelines
psa_plot_points2 <- p2 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "blue") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 100%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 100%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 100%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)") 

#Observed practice v. clinical guidelines
psa_plot_points3 <- p3 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "red") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 100%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 100%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 100%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA HRU Costs)
psa_plot_points4 <- p4 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "forestgreen") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 90.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 98.2%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.5%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA Crime Costs)
psa_plot_points5 <- p5 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "orange") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 98.1%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 99.5%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.9%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA REL-ABS)
psa_plot_points6 <- p6 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "darkorchid4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 97.7%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 98.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.1%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA All to REL)
psa_plot_points7 <- p7 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "green4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 96.5%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 96.7%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 96.8%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA DET to ABS)
psa_plot_points8 <- p8 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "brown4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 89.5%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 90.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 90.9%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#Healthcare perspective
#Evidence based v. observed practice
health_psa_plot_points <- health.p1 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point() +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 94.1%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 98.5%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.6%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)") 

#Evidence based v. clinical guidelines
health_psa_plot_points2 <- health.p2 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "blue") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 96.3%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 99.3%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.8%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)") 

#Observed practice v. clinical guidelines
health_psa_plot_points3 <- health.p3 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "red") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 97.3%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 99.5%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.8%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA HRU Costs)
health_psa_plot_points4 <- health.p4 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "forestgreen") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 1.0%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 33.4%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 81.9%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA Crime Costs)
health_psa_plot_points5 <- health.p5 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "orange") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 94.1%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 98.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 99.6%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA REL-ABS)
health_psa_plot_points6 <- health.p6 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "darkorchid4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 88.3%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 94.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 97.2%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA All-to-REL)
health_psa_plot_points7 <- health.p7 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "green4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 90.4%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 95.1%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 96.3%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

#EBvOP (PSA-OWSA DET-ABS)
health_psa_plot_points8 <- health.p8 +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_abline(mapping = NULL, data = NULL, slope = 50000, intercept = 0, na.rm = FALSE, show.legend = "ICER = $50,000") + # Add linear ICER line
  geom_abline(mapping = NULL, data = NULL, slope = 100000, intercept = 0, na.rm = FALSE, color = "navyblue", show.legend = "ICER = $100,000") + # Add linear ICER line
  geom_point(color = "brown4") +
  annotate("text", x = 2.50, y = 100000, label = "ICER = $50,000/QALY") +
  annotate("text", x = 0.45, y = 100000, color = "navyblue", label = "ICER = $100,000/QALY") +
  annotate("text", x = 2.5, y = -450000, label = "Proportion of simulations:") +
  annotate("text", x = 1.8, y = -470000, label = "Cost-saving = 76.6%", hjust = 0) +
  annotate("text", x = 1.8, y = -485000, label = "Cost-effective ($50,000/QALY) = 83.7%", hjust = 0) +
  annotate("text", x = 1.8, y = -500000, label = "Cost-effective ($100,000/QALY) = 86.8%", hjust = 0) +
  scale_x_continuous(breaks = c(-1,0,1,2,3), labels = c("-1","0","1","2","3"), limits = c(-1,3)) +
  scale_y_continuous(breaks = c(100000, 0, -100000, -200000, -300000, -400000, -500000), labels = c("$100,000", "$0", "-$100,000", "-$200,000", "-$300,000", "-$400,000", "-$500,000"), limits = c(-500000, 100000)) +
  xlab("Incremental QALYs") + ylab("Incremental Costs (2016 USD)")

ggsave(psa_plot_total_societal, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_total_societal_RR2.png", dpi=500)
ggsave(psa_plot_total_health, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_total_health_RR2.png", dpi=500)

ggsave(psa_plot_points, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_EBvOP_RR2.png", dpi=500)
ggsave(psa_plot_points2, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_EBvCG_RR2.png", dpi=500)
ggsave(psa_plot_points3, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_OPvCG_RR2.png", dpi=500)
#ggsave(psa_plot_points4, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_HRU_PSA_OWSA_RR2.png", dpi=500)
ggsave(psa_plot_points5, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_crime_PSA_OWSA_RR2.png", dpi=500)
ggsave(psa_plot_points6, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_REL-ABS_PSA_OWSA_RR2.png", dpi=500)
ggsave(psa_plot_points7, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_All-to-REL_PSA_OWSA_RR2.png", dpi=500)
ggsave(psa_plot_points8, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_DET-to-ABS_PSA_OWSA_RR2.png", dpi=500)

ggsave(health_psa_plot_points, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_EBvOP_RR2.png", dpi=500)
ggsave(health_psa_plot_points2, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_EBvCG_RR2.png", dpi=500)
ggsave(health_psa_plot_points3, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_OPvCG_RR2.png", dpi=500)
#ggsave(health_psa_plot_points4, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_HRU_PSA_OWSA_RR2.png", dpi=500)
ggsave(health_psa_plot_points5, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_crime_PSA_OWSA_RR2.png", dpi=500)
ggsave(health_psa_plot_points6, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_REL-ABS_PSA_OWSA_RR2.png", dpi=500)
ggsave(health_psa_plot_points7, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_All-to-REL_PSA_OWSA_RR2.png", dpi=500)
ggsave(health_psa_plot_points8, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_health_DET-to-ABS_PSA_OWSA_RR2.png", dpi=500)

ggsave(psa_plot_points_zoom, width=10, height=10, filename="Y:\\Projects\\Xnational project\\CEA paper\\FilesForBen\\CODE\\ce_plane_points_zoom.png", dpi=500)