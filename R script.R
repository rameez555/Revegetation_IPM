
## R scripts for reproducing the analysis presented in the paper


#### PACKAGES ####
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for data wrangling
library(metafor) #for meta-analysis
library(orchaRd) #for orchard plots
library(ggplot2) #for plotting
library(patchwork) # for combining plots
library(metaAidR) # for variance-covariance matric

#### DATA ####
# getting the data into R emvironment
full_data<-read.csv("data_all.csv")



#remove cases where control SD = 0,  that usually results in having vi=NA 
data_all<-full_data %>% filter(sdC != 0)

# Effect size calculation
# Calculating lnRR:

data.all.lnRR.sc<-escalc(measure="ROM", 
                         n1i=nT, n2i=nC, 
                         m1i=mT, m2i=mC, 
                         sd1i=sdT, sd2i=sdC,
                         data=data_all,
                         var.names=c("lnRR","lnRR.sv"), add.measure=FALSE,
                         append=FALSE)

#combined effect sizes with original data frame

data_es <-bind_cols(data_all, data.all.lnRR.sc)

# summarize studies by studied variables

# nativity status
data_es %>% group_by(Status_parameter1) %>% summarise(`Number of observations` = n(),
                                                      `Number of studies` = n_distinct(PR_ID))    

#  biodiversity metrics
data_es %>% group_by(Parameter_Category) %>% summarise(`Number of observations` = n(),
                                                       `Number of studies` = n_distinct(PR_ID))
# continent
data_es %>% group_by(Continent) %>% summarise(`Number of observations` = n(),
                                              `Number of studies` = n_distinct(PR_ID))

# ecosystem
data_es %>% group_by(Ecosystem) %>% summarise(`Number of observations` = n(),
                                              `Number of studies` = n_distinct(PR_ID))

# re-vegetation method
data_es %>% group_by(Revegetation.method) %>% summarise(`Number of observations` = n(),
                                                        `Number of studies` = n_distinct(PR_ID))

# functional group
data_es %>% group_by(Method_GF_Fun) %>% summarise(`Number of observations` = n(),
                                                        `Number of studies` = n_distinct(PR_ID))
# taxonomic composition
data_es %>% group_by(revegetation_type) %>% summarise(`Number of observations` = n(),
                                                      `Number of studies` = n_distinct(PR_ID))
# growth form
data_es %>% group_by(Method_GF_Gen) %>% summarise(`Number of observations` = n(),
                                                      `Number of studies` = n_distinct(PR_ID))


# make vcv matrix with 0.5 correlation
data_es_rr <- data_es[order(data_es$CHT_ID), ]
vcv_RR_0.5 <- make_VCV_matrix(data = data_es_rr, V = "lnRR.sv",
                              cluster = "CHT_ID", obs = "ES_ID",
                              type = "vcv", rho = 0.5)


# overall
lnRR_full <- rma.mv(lnRR, vcv_RR_0.5, random = list(~1 | PR_ID, ~1 | ES_ID, ~1 | CHT_ID), method = "REML", data = data_es_rr)



summary(lnRR_full, digits = 3) #summarise model results
round(i2_ml(lnRR_full), 5) 

# CHT_ID as random effect not explaining any heterogeneity
# re-run model while dropping CHT_ID from random factors

lnRR_full <- rma.mv(lnRR, vcv_RR_0.5, random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)

summary(lnRR_full, digits = 3) #summarise model results
round(i2_ml(lnRR_full), 2) 


full<-mod_results(lnRR_full, mod = "1", data = data_es_rr, group = "PR_ID")


orchard_plot(full, 
             mod="lnRR_full", 
             xlab = "log(Response ratio) (lnRR)",
             alpha = 0.75
) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position= c(0.02, 0.99), 
    legend.justification = c(0,1), 
    legend.key.size = unit(1, "mm"),
    legend.direction="horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9),
    panel.grid = element_blank())


ggsave("orchard_overall.jpg", height=3.5, width=7, units='in', dpi=1200)


# ORCHARD PLOT FUNCTION ----
# this function modifies aesthetics in subsequent orchard plots
orchard_plot <- function(object, mod = "1", group, xlab, N = NULL,
                         alpha = 0.5, angle = 90, cb = TRUE, k = TRUE, g = TRUE,
                         tree.order = NULL, trunk.size = 3, branch.size = 0.5, twig.size = 0.5,
                         transfm = c("none", "tanh"), condition.lab = "Condition",
                         legend.pos = c("bottom.right", "bottom.left",
                                        "top.right", "top.left",
                                        "top.out", "bottom.out",
                                        "none"), # "none" - no legends
                         k.pos = c("right", "left", "none"),
                         colour = FALSE,
                         fill = TRUE,
                         weights = "prop", by = NULL, at = NULL, upper = TRUE, flip = TRUE)
{
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  
  if(any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))){
    
    if(mod != "1"){
      results <-  orchaRd::mod_results(object, mod, group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    } else {
      results <-  orchaRd::mod_results(object, mod = "1", group,  N,
                                       by = by, at = at, weights = weights, upper = upper)
    }
  }
  
  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }
  
  mod_table <- results$mod_table
  
  data_trim <- results$data
  # making sure factor names match
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)
  
  data_trim$scale <- (1/sqrt(data_trim[,"vi"]))
  legend <- "Precision (1/SE)"
  
  #if tree.order isn't equal to NULL, and length of tree order does not match number of categories in categorical moderator, then stop function and throw an error
  if(!is.null(tree.order)&length(tree.order)!=nlevels(data_trim[,'moderator'])){
    stop("Length of 'tree.order' does not equal number of categories in moderator")
  }
  
  #if tree.order isn't equal to NULL but passes above check, then reorder mod table according to custom order if there is one
  if (!is.null(tree.order)){
    data_trim$moderator<-factor(data_trim$moderator, levels = tree.order, labels = tree.order)
    mod_table <- mod_table %>% dplyr::arrange(factor(name, levels = tree.order))
  }
  
  if(is.null(N) == FALSE){
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
    #latex2exp::TeX()
  }
  
  if(transfm == "tanh"){
    cols <- sapply(mod_table, is.numeric)
    mod_table[,cols] <- Zr_to_r(mod_table[,cols])
    data_trim$yi <- Zr_to_r(data_trim$yi)
    label <- xlab
  }else{
    label <- xlab
  }
  
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[,"moderator"], function(x) length(x[,"yi"])))
  
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[,2])
  
  # the number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))
  
  #data_no <- nrow(data)
  
  cbpl <- c("#999999", "#999999", "#999999", 
            "#999999", "#999999", "#999999", "#999999")
  
  # setting fruit colour
  if(colour == TRUE){
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  }else{
    color <- data_trim$mod
    color2 <- mod_table$name
  }
  
  # whether we fill fruit or not
  if(fill == TRUE){
    fill <- color
  }else{
    fill <- NULL
  }
  
  # whether marginal
  if(names(mod_table)[2] == "condition"){
    
    # the number of levels in the condition
    condition_no <- length(unique(mod_table[, "condition"]))
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = moderator), shape = 21, fill = alpha("white", 0.7), groupOnX = FALSE) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", size = 0.3) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
                              size = branch.size, position = ggplot2::position_dodge2(width = 0.3)) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name, ymin = lowerCL, ymax = upperCL,  shape = as.factor(condition), fill = color2), size = twig.size, position = ggplot2::position_dodge2(width = 0.3), fatten = trunk.size) +
      # this will only work for up to 5 different conditions
      # flipping things around (I guess we could do use the same geoms but the below is the original so we should not change)
      ggplot2::scale_shape_manual(values =  20 + (1:condition_no))  +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::labs(shape = condition.lab) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle))
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
    
  } else {
    
    plot <- ggplot2::ggplot() +
      # pieces of fruit (bee-swarm and bubbles)
      ggbeeswarm::geom_quasirandom(data = data_trim, ggplot2::aes(y = yi, x = moderator, size = scale, colour = moderator), shape = 21, fill = alpha("white", 0.7), groupOnX = FALSE) +
      
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "black", size = 0.3) +
      # creating CI
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL), size = branch.size) +
      # drowning point estimate and PI
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, x = name,  ymin = lowerCL, ymax = upperCL, fill = color2), size = twig.size, fatten = trunk.size, shape = 21) +
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none", colour = "none") +
      ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) +
      ggplot2::theme(legend.direction="horizontal") +
      ggplot2::theme(legend.background = ggplot2::element_blank()) +
      ggplot2::labs(y = label, x = "", size = latex2exp::TeX(legend)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, colour ="black",
                                                         hjust = 0.5,
                                                         angle = angle)) #+
    #ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
    if(flip){
      plot <- plot + ggplot2::coord_flip()
    }
  }
  
  # adding legend
  if(legend.pos == "bottom.right"){
    plot <- plot + ggplot2::theme(legend.position= c(1, 0), legend.justification = c(1, 0))
  } else if ( legend.pos == "bottom.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 0), legend.justification = c(0, 0))
  } else if ( legend.pos == "top.right") {
    plot <- plot + ggplot2::theme(legend.position= c(1, 1), legend.justification = c(1, 1))
  } else if (legend.pos == "top.left") {
    plot <- plot + ggplot2::theme(legend.position= c(0, 1), legend.justification = c(0, 1))
  } else if (legend.pos == "top.out") {
    plot <- plot + ggplot2::theme(legend.position="top")
  } else if (legend.pos == "bottom.out") {
    plot <- plot + ggplot2::theme(legend.position="bottom")
  } else if (legend.pos == "none") {
    plot <- plot + ggplot2::theme(legend.position="none")
  }
  
  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values = cbpl) +
      ggplot2::scale_colour_manual(values = cbpl)
  }
  
  # putting k and g in
  if(k == TRUE && g == FALSE && k.pos == "right"){
    plot <- plot +
      ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                        label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "right", size = 3.5)
  } else if(k == TRUE && g == FALSE && k.pos == "left") {
    plot <- plot +  ggplot2::annotate('text', y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                      label= paste("italic(k)==", mod_table$K[1:group_no]), parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "right"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text', y = (max(data_trim$yi) + (max(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "right", size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos == "left"){
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate('text',  y = (min(data_trim$yi) + (min(data_trim$yi)*0.10)), x = (seq(1, group_no, 1)+0.3),
                                     label= paste("italic(k)==", mod_table$K[1:group_no], "~","(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, hjust = "left", size = 3.5)
  } else if (k == TRUE && g == FALSE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no]),
                                     parse = TRUE, size = 3.5)
  } else if (k == TRUE && g == TRUE && k.pos%in%c('right','left','none')==FALSE) {
    # get group numbers for moderator
    plot <- plot + ggplot2::annotate("text", y = k.pos, x = (seq(1, group_no,
                                                                 1) + 0.3), label = paste("italic(k)==", mod_table$K[1:group_no],
                                                                                          "~", "(", mod_table$g[1:group_no], ")"),
                                     parse = TRUE, size = 3.5)
  }
  return(plot)
}


#plant nativity status
lnRR_full_Status_parameter <- rma.mv(yi = lnRR, V = vcv_RR_0.5, mods = ~Status_parameter1 -1, 
                                     random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)
lnRR_Status_parameter<-mod_results(lnRR_full_Status_parameter, mod = "Status_parameter1", data = data_es_rr, group = "PR_ID")

round(r2_ml(lnRR_full_Status_parameter)*100, 2)

# drawing plot
orchard_plot(lnRR_Status_parameter, mod = "Status_parameter1", xlab = "log(Response ratio) (lnRR)", angle = 45, cb = TRUE, alpha = 0.1, group = "PR_ID", legend.pos = "bottom.right") + 
  theme(
    legend.direction="horizontal",
    panel.grid = element_blank())
ggsave("orchard_nativity.jpg", height=3.5, width=5, units='in', dpi=1200)


## Moderators: Functional group vs Status parameter1

lnRR_full_Method_GF_Fun_cat <- rma.mv(yi = lnRR, V = vcv_RR_0.5, mods = ~Method_GF_Fun -1 + Status_parameter1, 
                                       random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)
lnRR_Method_GF_Fun_cat<- mod_results(lnRR_full_Method_GF_Fun_cat, group = "PR_ID", mod = "Method_GF_Fun", by = "Status_parameter1", weights = "prop" )

# drawing plots
orchard_plot(lnRR_Method_GF_Fun_cat, mod = "Method_GF_Fun", xlab = "log(Response ratio) (lnRR)", angle = 45, cb = TRUE, alpha = 0.1, group = "PR_ID", g = TRUE, legend.pos = "top.left", k.pos = "right", condition.lab = "Status") + 
  theme(
    legend.direction="vertical",
    panel.grid = element_blank())
ggsave("Method_GF_Fun_cat.jpg", height=4, width=6, units='in', dpi=1200)

(exp(0.44) - 1) * 100
(exp(0.33) - 1) * 100
(exp(-0.38) - 1) * 100

round(r2_ml(lnRR_full_Method_GF_Fun_cat)*100, 2)

## Moderators: re-vegetation type (taxonomic) vs Status parameter1

lnRR_full_revegetation_type <- rma.mv(yi = lnRR, V = vcv_RR_0.5, mods = ~revegetation_type -1 + Status_parameter1, 
                                      random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)
lnRR_Method_revegetation_type<- mod_results(lnRR_full_revegetation_type, group = "PR_ID", mod = "revegetation_type", by = "Status_parameter1", weights = "prop" )

# drawing plots
orchard_plot(lnRR_Method_revegetation_type, mod = "revegetation_type", xlab = "log(Response ratio) (lnRR)", angle = 45, cb = TRUE, alpha = 0.1, group = "PR_ID", g = TRUE, legend.pos = "top.left", k.pos = "right", condition.lab = "Status") + 
  theme(
    legend.direction="vertical",
    panel.grid = element_blank())
ggsave("Method_revegetation_type.jpg", height=4, width=6, units='in', dpi=1200)

round(r2_ml(lnRR_full_revegetation_type)*100, 2)

## Moderators: growth form vs Status parameter1
data_es_rr$Method_GF_Gen = as.factor(data_es_rr$Method_GF_Gen)

data_es_rr$Method_GF_Gen_reordered <- factor(data_es_rr$Method_GF_Gen, c("Herbaceous and Woody", "Woody", "Herbaceous"))


lnRR_full_Method_GF_Gen <- rma.mv(yi = lnRR, V = vcv_RR_0.5, mods = ~Method_GF_Gen_reordered -1 + Status_parameter1, 
                                  random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)
lnRR_Method_Method_GF_Gen<- mod_results(lnRR_full_Method_GF_Gen, group = "PR_ID", mod = "Method_GF_Gen_reordered", by = "Status_parameter1", weights = "prop" )

# drawing plots
orchard_plot(lnRR_Method_Method_GF_Gen, mod = "Method_GF_Gen_reordered", xlab = "log(Response ratio) (lnRR)", angle = 45, cb = TRUE, alpha = 0.1, group = "PR_ID", g = TRUE, legend.pos = "top.left", k.pos = "right", condition.lab = "Status") + 
  theme(
    legend.direction="vertical",
    panel.grid = element_blank())
ggsave("Method_GF_Gen.jpg", height=4, width=6, units='in', dpi=1200)

round(r2_ml(lnRR_full_Method_GF_Gen)*100, 2)


uni_mod_plot_ns<-function(m, df, log_ratio, response, variance){
  p <- predict.rma(m)
  df %>% mutate(ymin = p$ci.lb, 
                ymax = p$ci.ub, ymin2 = p$cr.lb, 
                ymax2 = p$cr.ub, pred = p$pred) %>% 
    ggplot(aes(x = response, y = log_ratio, size = sqrt(1/variance))) + geom_point(shape = 21, alpha= 0.2,
                                                                                   fill = "grey90") + 
    geom_hline(yintercept = 0, size = .5, colour = "gray70")+
    geom_smooth(aes(y = ymin2), method = "lm", se = FALSE, lty = "dashed", lwd = 0.75, 
                colour = "#0072B2") + geom_smooth(aes(y = ymax2), method = "lm", se = FALSE, 
                                                  lty = "dashed", lwd = 0.75, colour = "#0072B2") + geom_smooth(aes(y = ymin), 
                                                                                                                method = "lm", se = FALSE, lty = "dashed", lwd = 0.75, colour = "#D55E00") + 
    geom_smooth(aes(y = ymax), method = "lm", se = FALSE, lty = "dashed", lwd = 0.75, 
                colour = "#D55E00") + geom_smooth(aes(y = pred), method = "lm", se = FALSE, 
                                                  lty = "dashed", lwd = 1, colour = "black") + 
    theme_classic() + theme(legend.position = c(0, 1), legend.justification = c(0, 1)) + theme(legend.direction = "horizontal") + 
    theme(legend.background = element_blank()) + theme(axis.text.y = element_text(size = 8, 
                                                                                  colour = "black", hjust = 0.5, angle = 90))+
    coord_cartesian(ylim = c(-3.5, 3.5))+
    scale_y_continuous(limits = c(-3.5, 3.5),
                       breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    theme(legend.position = "none") 
}

# duration
data_es_rr_temp_native <- data_es_rr %>% filter(Status_parameter1 == "Native")
data_es_rr_temp_native <- data_es_rr_temp_native %>% drop_na(Duration_months._FS)
data_es_rr_temp_native <- data_es_rr_temp_native[order(data_es_rr_temp_native$CHT_ID), ]
vcv_RR_0.5_temp_native <- make_VCV_matrix(data = data_es_rr_temp_native, V = "lnRR.sv",
                              cluster = "CHT_ID", obs = "ES_ID",
                              type = "vcv", rho = 0.5)

lat_effect_uni_temporal <- rma.mv(yi = lnRR, V = vcv_RR_0.5_temp_native, mods = ~Duration_months._FS, test = "t", 
                              random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr_temp_native)
a1<-uni_mod_plot_ns(lat_effect_uni_temporal, data_es_rr_temp_native, log_ratio = data_es_rr_temp_native$lnRR, response = data_es_rr_temp_native$Duration_months._FS, variance = data_es_rr_temp_native$lnRR.sv)+ylab("Log (Response ratio) (lnRR) ")+xlab("Temporal duration (months)")


data_es_rr_temp_non_native <- data_es_rr %>% filter(Status_parameter1 == "Non-native")
data_es_rr_temp_non_native <- data_es_rr_temp_non_native %>% drop_na(Duration_months._FS)
data_es_rr_temp_non_native <- data_es_rr_temp_non_native[order(data_es_rr_temp_non_native$CHT_ID), ]
vcv_RR_0.5_temp_non_native <- make_VCV_matrix(data = data_es_rr_temp_non_native, V = "lnRR.sv",
                                          cluster = "CHT_ID", obs = "ES_ID",
                                          type = "vcv", rho = 0.5)

lat_effect_uni_temporal_non_native <- rma.mv(yi = lnRR, V = vcv_RR_0.5_temp_non_native, mods = ~Duration_months._FS, test = "t", 
                                  random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr_temp_non_native)
a2<-uni_mod_plot_ns(lat_effect_uni_temporal_non_native, data_es_rr_temp_non_native, log_ratio = data_es_rr_temp_non_native$lnRR, response = data_es_rr_temp_non_native$Duration_months._FS, variance = data_es_rr_temp_non_native$lnRR.sv)+ylab("Log (Response ratio) (lnRR) ")+xlab("Temporal duration (months)")

(a1 + a2)  + plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
ggsave("Duration.jpg", height=4, width=7, units='in', dpi=1200)
# getting marginal R2
round(r2_lat_effect_uni_temporal <- r2_ml(lat_effect_uni_temporal)*100, 2)
round(r2_lat_effect_uni_temporal_non_native <- r2_ml(lat_effect_uni_temporal_non_native)*100, 2)

# temporal
data_es_rr_spatial_native <- data_es_rr %>% filter(Status_parameter1 == "Native")
data_es_rr_spatial_native <- data_es_rr_spatial_native %>% drop_na(Plot_size_m2)
data_es_rr_spatial_native <- data_es_rr_spatial_native[order(data_es_rr_spatial_native$CHT_ID), ]
vcv_RR_0.5_spatial_native <- make_VCV_matrix(data = data_es_rr_spatial_native, V = "lnRR.sv",
                                          cluster = "CHT_ID", obs = "ES_ID",
                                          type = "vcv", rho = 0.5)

lat_effect_uni_spatial <- rma.mv(yi = lnRR, V = vcv_RR_0.5_spatial_native, mods = ~log10(Plot_size_m2), test = "t", 
                                  random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr_spatial_native)
a1<-uni_mod_plot_ns(lat_effect_uni_spatial, data_es_rr_spatial_native, log_ratio = data_es_rr_spatial_native$lnRR, response = log10(data_es_rr_spatial_native$Plot_size_m2), variance = data_es_rr_spatial_native$lnRR.sv)+ylab("Log (Response ratio) (lnRR) ")+xlab(expression("Spatial scale" (log(m^2))))


data_es_rr_spatial_non_native <- data_es_rr %>% filter(Status_parameter1 == "Non-native")
data_es_rr_spatial_non_native <- data_es_rr_spatial_non_native %>% drop_na(Plot_size_m2)
data_es_rr_spatial_non_native <- data_es_rr_spatial_non_native[order(data_es_rr_spatial_non_native$CHT_ID), ]
vcv_RR_0.5_spatial_non_native <- make_VCV_matrix(data = data_es_rr_spatial_non_native, V = "lnRR.sv",
                                              cluster = "CHT_ID", obs = "ES_ID",
                                              type = "vcv", rho = 0.5)

lat_effect_uni_spatial_non_native <- rma.mv(yi = lnRR, V = vcv_RR_0.5_spatial_non_native, mods = ~log10(Plot_size_m2), test = "t", 
                                             random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr_spatial_non_native)
a2<-uni_mod_plot_ns(lat_effect_uni_spatial_non_native, data_es_rr_spatial_non_native, log_ratio = data_es_rr_spatial_non_native$lnRR, response = log10(data_es_rr_spatial_non_native$Plot_size_m2), variance = data_es_rr_spatial_non_native$lnRR.sv)+ylab("Log (Response ratio) (lnRR) ")+xlab(expression("Spatial scale" (log(m^2))))

(a1 + a2)  + plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
ggsave("Plot_size.jpg", height=4, width=7, units='in', dpi=1200)
# getting marginal R2
round(r2_lat_effect_uni_spatial <- r2_ml(lat_effect_uni_spatial)*100, 2)
round(r2_lat_effect_uni_spatial_non_native <- r2_ml(lat_effect_uni_spatial_non_native)*100, 3)

# univariate egger regressoin

egger_uni_rr <- rma.mv(yi = lnRR, V = lnRR.sv, mods = ~sqrt(lnRR.sv), test = "t", random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)



### plotting funnel plots


funnel(egger_uni_rr, yaxis = "seinv", xlim = c(-7, 7), ylim = c(0.05, 1.2),level = c(90, 95, 99), shade = c("white", "gray55", "gray75"), refline = 0, legend = TRUE)
# getting marginal R2 for lnRR
round(r2_egger_uni_rr <- r2_ml(egger_uni_rr)*100, 2)

# temporal effect
time_lag_effect_uni_lnRR <- rma.mv(yi = lnRR, V = vcv_RR_0.5, mods = ~Year, test = "t", 
                                   random = list(~1 | PR_ID, ~1 | ES_ID), method = "REML", data = data_es_rr)
# getting marginal R2
round(r2_time_lag_effect_uni_lnRR <- r2_ml(time_lag_effect_uni_lnRR)*100, 2)



#plotting figures

pred_pub_lnRR <- predict.rma(egger_uni_rr)
pub_all_lnRR <- data_es_rr %>% mutate(ymin = pred_pub_lnRR$ci.lb, ymax = pred_pub_lnRR$ci.ub,
                                      ymin2 = pred_pub_lnRR$cr.lb, ymax2 = pred_pub_lnRR$cr.ub, pred = pred_pub_lnRR$pred) %>%
  ggplot(aes(x = sqrt(lnRR.sv), y = lnRR, size = sqrt(1/lnRR.sv))) + geom_point(shape = 21,
                                                                                fill = "grey90") + geom_smooth(aes(y = ymin2), method = "loess",
                                                                                                               se = FALSE, lty = "dotted", lwd = 0.75, colour = "#0072B2") +
  geom_smooth(aes(y = ymax2), method = "loess", se = FALSE,
              lty = "dotted", lwd = 0.75, colour = "#0072B2") + geom_smooth(aes(y = ymin),
                                                                            method = "loess", se = FALSE, lty = "dotted", lwd = 0.75,
                                                                            colour = "#D55E00") + geom_smooth(aes(y = ymax), method = "loess",
                                                                                                              se = FALSE, lty = "dotted", lwd = 0.75, colour = "#D55E00") +
  geom_smooth(aes(y = pred), method = "loess", se = FALSE,
              lty = "dashed", lwd = 1, colour = "black") + ylim(-2.5,
                                                                2.5) + xlim(0, 2) + labs(x = "sqrt(sampling variance)",
                                                                                            y = "log(Response ratio) (lnRR)", size = "Precision (1/SE)") + guides(fill = "none",
                                                                                                                                                                  colour = "none") + # themes
  theme_bw() + theme(legend.position = c(0, 1), legend.justification = c(0,
                                                                         1)) + theme(legend.direction = "horizontal") + theme(legend.background = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black",
                                   hjust = 0.5, angle = 90))

pred_time <- predict.rma(time_lag_effect_uni_lnRR)
time_lag_lnRR <- data_es_rr %>% mutate(ymin = pred_time$ci.lb, ymax = pred_time$ci.ub,
                                       ymin2 = pred_time$cr.lb, ymax2 = pred_time$cr.ub, pred = pred_time$pred) %>%
  ggplot(aes(x = Year, y = lnRR, size = sqrt(1/lnRR.sv))) + geom_point(shape = 21,
                                                                       fill = "grey90") + 
  geom_smooth(aes(y = ymin2), method = "loess", se = FALSE, lty = "dotted",
              lwd = 0.75, colour = "#0072B2") + geom_smooth(aes(y = ymax2),
                                                            method = "loess", se = FALSE, lty = "dotted", lwd = 0.75,
                                                            colour = "#0072B2") + geom_smooth(aes(y = ymin), method = "loess",
                                                                                              se = FALSE, lty = "dotted", lwd = 0.75, colour = "#D55E00") +
  geom_smooth(aes(y = ymax), method = "loess", se = FALSE,
              lty = "dotted", lwd = 0.75, colour = "#D55E00") + geom_smooth(aes(y = pred),
                                                                            method = "loess", se = FALSE, lty = "dashed", lwd = 1, colour = "black") +
  ylim(-2.5, 2.5) + xlim(1998, 2023) + labs(x = "Year", y = "log(Response ratio) (lnRR)",
                                            size = "Precision (1/SE)") + guides(fill = "none", colour = "none") +
  # themes
  theme_bw() + theme(legend.position = c(0, 1), legend.justification = c(0,
                                                                         1)) + theme(legend.direction = "horizontal") + theme(legend.background = element_blank()) +
  theme(axis.text.y = element_text(size = 10, colour = "black",
                                   hjust = 0.5, angle = 90))
(pub_all_lnRR + time_lag_lnRR)  + plot_annotation(tag_prefix = "(", tag_levels = "a", tag_suffix = ")")
ggsave("Bias.jpg", height=4, width=7, units='in', dpi=1200)

# R session info
library(pander)
sessionInfo() %>% pander()
