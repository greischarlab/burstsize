
########################################
### OPTIMAL TRANSMISSION INVESTMENT #####
########################################

FITNESS_ALL_OPTIMAL_CV <- NULL
for (k in seq(1, nrow(FITNESS_MAX))) {
  tmp <- FITNESS_MAX[k, ]
  tmp_fitness <- FITNESS_ALL_SPLIT[[k]]
  
  optimal_cv <- subset(tmp_fitness, tmp_fitness$C_V == tmp$C_V &
                         tmp_fitness$status == "success")
  
  FITNESS_ALL_OPTIMAL_CV[[k]] <- cbind.data.frame(
    B_V = optimal_cv$B_V,
    C_V = optimal_cv$C_V,
    change = optimal_cv$change,
    grouping = optimal_cv$grouping,
    endfitness = optimal_cv$end_fitness
  )
}


FITNESS_ALL_OPTIMAL_CV_F <- do.call(rbind, FITNESS_ALL_OPTIMAL_CV)

###THIS IS THE Data-frame that does not include the original sim
FITNESS_ALL_OPTIMAL_CV_1 <- subset(
  FITNESS_ALL_OPTIMAL_CV_F,
  FITNESS_ALL_OPTIMAL_CV_F$change != 1
)


FITNESS_ALL_OPTIMAL_CV_2 <- subset(FITNESS_ALL_OPTIMAL_CV_F, FITNESS_ALL_OPTIMAL_CV_F$C_V == 0.75 &
                                     FITNESS_ALL_OPTIMAL_CV_F$grouping == "Original")[, c("B_V", "endfitness")]


OPT_CV_GG <- ggplot(
  FITNESS_ALL_OPTIMAL_CV_1,
  aes(
    x = B_V,
    y = endfitness,
    group =as.factor(change),
    color = as.factor(change)
  )
) +
  geom_line(size = 1.2) +
  facet_wrap(~grouping,
             ncol = 1,
             labeller = as_labeller(Parameter_Names)
  ) +
  annotate("line",
           x = FITNESS_ALL_OPTIMAL_CV_2$B_V,
           y = FITNESS_ALL_OPTIMAL_CV_2$endfitness,
           size = 1.2
  ) +
  scale_color_manual(values = c(
    "#118ab2", "#ef233c",
    "#118ab2", "#ef233c",
    "#118ab2", "#ef233c"
  )) +
  xlab("Burst size") +
  ylab("Cumulative transmission potential") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15)
  )

########################################
### Optimal burst size              ####
########################################
FITNESS_ALL_OPTIMAL_BV <- NULL
for (k in seq(1, nrow(FITNESS_MAX))) {
  tmp <- FITNESS_MAX[k, ]
  tmp_fitness <- FITNESS_ALL_SPLIT[[k]]
  
  optimal_bv <- subset(tmp_fitness, tmp_fitness$B_V == tmp$B_V &
                         tmp_fitness$status == "success")
  
  
  
  FITNESS_ALL_OPTIMAL_BV[[k]] <- cbind.data.frame(
    B_V = optimal_bv$B_V,
    C_V = optimal_bv$C_V,
    change = optimal_bv$change,
    grouping = optimal_bv$grouping,
    endfitness = optimal_bv$end_fitness
  )
}

FITNESS_ALL_OPTIMAL_BV_F <- do.call(rbind, FITNESS_ALL_OPTIMAL_BV)

FITNESS_ALL_OPTIMAL_BV_1 <- subset(
  FITNESS_ALL_OPTIMAL_BV_F,
  FITNESS_ALL_OPTIMAL_BV_F$change != 1
)


FITNESS_ALL_OPTIMAL_BV_2 <- subset(
  FITNESS_ALL_OPTIMAL_BV_F,
  round(FITNESS_ALL_OPTIMAL_BV_F$B_V,2) == 14.5 & 
    FITNESS_ALL_OPTIMAL_BV_F$grouping == 'Original')[, c("C_V", "endfitness")]




OPT_BV_GG <- ggplot(
  FITNESS_ALL_OPTIMAL_BV_1,
  aes(
    x = C_V,
    y = endfitness,
    group = as.factor(change),
    color = as.factor(change)
  )
) +
  geom_line(size = 1.2) +
  facet_wrap(~grouping,
             ncol = 1,
             labeller = as_labeller(Parameter_Names)
  ) +
  annotate("line",
           x = FITNESS_ALL_OPTIMAL_BV_2$C_V,
           y = FITNESS_ALL_OPTIMAL_BV_2$endfitness,
           size = 1.2
  ) +
  scale_color_manual(values = c(
    "#118ab2", "#ef233c",
    "#118ab2", "#ef233c",
    "#118ab2", "#ef233c"
  )) +
  xlab("Transmission investment") +
  ylab("Cumulative transmission potential") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "off",
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 15)
  )


OPT_CV_GG + OPT_BV_GG


ggsave(here("Figures", "Raw", "Supp_Optimal_BV_CV_SA.pdf"),
       height = 8, width = 9, units = "in"
)
