
sample_heatmap_plot <- function(data_source, 
                                genes = c("Hspa8","Snap25","Gad2","Vip"), 
                                group_by = "cluster", 
                                groups = 1:10,
                                top_labels = "layer",
                                sample_n = 0,
                                scale_mode = "scale.log",
                                showall = F,
                                autorange = "auto", 
                                minrange = 0, 
                                maxrange = 10,
                                pfontsize = 14,
                                expand = F,
                                rotatelabels = F,
                                showlines = F,
                                showids = T) {
  
  library(feather)
  library(dplyr)
  library(ggplot2)
  
  data <- read_feather(file.path(data_source, "data.feather"), columns = c("sample_id",genes))
  anno <- read_feather(file.path(data_source, "anno.feather"))
  desc_table <- read_feather(file.path(data_source, "desc.feather"))
  
  primary <- list(base = group_by,
                  id = paste0(group_by,"_id"),
                  label = paste0(group_by,"_label"),
                  color = paste0(group_by,"_color"))
  secondary <- list(base = top_labels,
                    id = paste0(top_labels,"_id"),
                    label = paste0(top_labels,"_label"),
                    color = paste0(top_labels,"_color"))
  
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  genes.df <- data

  ###############
  ## Filtering ##
  ###############
  
  # Join the annotation and genes data frames
  sub.df <- anno %>% 
    # Filter for the selected clusters
    filter_(paste0(primary$id, " %in% c(", paste(groups, collapse = ","), ")")) %>%
    left_join(genes.df)
  
  # Subsample data per primary group if Sample N is set to something other than 0
  if(sample_n > 0) {
    sub.df <- sub.df %>%
      group_by_(primary$id) %>%
      sample_n(sample_n, replace = T) %>%
      unique() %>%
      ungroup() %>%
      as.data.frame()
  }
  
  #############
  ## Scaling ##
  #############

  ## If the y-axis is plotted on a log scale, add 1 to the data values to plot data + 1
  if(scale_mode == "scale.log") {
    for(gene in genes) {
      sub.df[,gene] <- log10(sub.df[,gene] + 1)
    }
  }
  if(scale_mode == "scale.rel") {
    for(gene in genes) {
      sub.df[,gene] <- sub.df[,gene]/max(sub.df[,gene])
    }
  }
  if(scale_mode == "scale.log.rel") {
    for(gene in genes) {
      sub.df[,gene] <- log10(sub.df[,gene]+1)/log10((max(sub.df[,gene]+1)))
    }
  }
  
  #############
  ## Sorting ##
  #############
  
  cluster_order <- data.frame(clust = groups,
                              plot_order = 1:length(groups))
  
  names(cluster_order)[1] <- primary$id
  
  sub.df <- sub.df %>% 
    left_join(cluster_order, by = primary$id)
  
  sort.df <- sub.df %>% 
    arrange_(.dots = c("plot_order", secondary$id)) %>% 
    mutate(xpos = 1:nrow(sub.df))
  
  
  
  # Start buildplot
  genes <- sub("-", ".", genes)
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  names(sort.df) <- sub("-",".",names(sort.df))
  
  colors <- colorRampPalette(c("darkblue", "white", "red"))(1001)
  
  if(autorange == "auto") {
    min.val <- 0
    max.val <- max(unlist(sort.df[ ,genes]))
  } else if (autorange == "manual") {
    min.val <- as.numeric(minrange)
    max.val <- as.numeric(maxrange)
  }
  
  ## Convert data to geom_rect() compatible table
  plot.df <- data.frame(xmin=numeric(),xmax=numeric(),ymin=numeric(),ymax=numeric(),fill=character())
  
  for(i in 1:length(genes)) {
    
    fill_ids <- unlist(round( (sort.df[,genes[i]] - min.val) / (max.val - min.val) * 1000 ) + 1)
    fill_ids[fill_ids < 1] <- 1
    fill_ids[fill_ids > 1001] <- 1001
    
    gene.plot <- data.frame(xmin = sort.df$xpos - 1,
                            xmax = sort.df$xpos,
                            ymin = length(genes) - i,
                            ymax = length(genes) - i + 1,
                            fill = colors[fill_ids])
    
    plot.df <- rbind(plot.df, gene.plot)
    
  }
  
  primary.plot <- data.frame(xmin = sort.df$xpos - 1,
                             xmax = sort.df$xpos,
                             ymin = -0.5,
                             ymax = 0,
                             fill = unlist(sort.df[ ,primary$color]))
  
  
  
  plot.df <- rbind(plot.df, primary.plot)
  
  ## add additional secondary color bars
  all.desc <- desc_table
  primary.name   <- all.desc$name[all.desc$base == primary$base]
  secondary.name <- all.desc$name[all.desc$base %in% secondary$base]
  
  primary.desc <- all.desc[all.desc$base == primary$base,]
  secondary.desc <- all.desc[all.desc$base %in% secondary$base,]
  other.desc <- all.desc[!all.desc$base %in% c(primary$base,secondary$base),]
  
  sec.color <- paste0(secondary$base, "_color")
  anno.color <- paste0(other.desc$base, "_color")
  anno_y_labels <- data.frame(breaks = numeric(), labels = character())
  
  if(showall) {
    
    # scale the plot so it's evenly divided between annotations and genes
    #s <- length(genes)/nrow(other.desc)*0.75
    s <- 1
    
    for(j in 1:length(secondary$base)) {
      anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                              xmax = sort.df$xpos,
                              ymin = length(genes) + (j - 1) * s,
                              ymax = length(genes) + (j - 1) * s + s,
                              fill = unlist(sort.df[ ,sec.color[j]]))
      plot.df <- rbind(plot.df, anno.plot)
      
      anno_y <- data.frame(breaks = length(genes) + (j - 1)*s + 0.5*s,
                           labels = secondary.desc$name[secondary.desc$base == secondary$base[j]])
      
      anno_y_labels <- rbind(anno_y_labels,anno_y)
    }
    
    sec_top <- length(secondary$base)
    
    for(j in 1:nrow(other.desc)) {
      anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                              xmax = sort.df$xpos,
                              ymin = length(genes) + (j + sec_top - 1) * s,
                              ymax = length(genes) + (j + sec_top - 1) * s + s,
                              fill = unlist(sort.df[ ,anno.color[j]]))
      plot.df <- rbind(plot.df, anno.plot)
      
      anno_y <- data.frame(breaks = length(genes) + (j + sec_top - 1)*s + 0.5*s,
                           labels = other.desc$name[j])
      
      anno_y_labels <- rbind(anno_y_labels,anno_y)
    }
  } else {
    
    if(length(secondary) > 0) {
      
      for(j in 1:length(secondary$base)) {
        anno.plot <- data.frame(xmin = sort.df$xpos - 1,
                                xmax = sort.df$xpos,
                                ymin = length(genes) + (j - 1) * 0.5,
                                ymax = length(genes) + (j - 1) * 0.5 + 0.5,
                                fill = unlist(sort.df[,sec.color[j]]))
        plot.df <- rbind(plot.df,anno.plot)
        
        anno_y <- data.frame(breaks = length(genes) + (j - 1) * 0.5 + 0.25,
                             labels = secondary.desc$name[secondary.desc$base == secondary$base[j]])
        anno_y_labels <- rbind(anno_y_labels,anno_y)
      }
    }
  }
  
  ## build new, more complex y-axis labels
  y_labels <- data.frame(breaks = (1:length(genes) - 0.5),
                         labels = rev(genes))
  primary_y_label <- data.frame(breaks = -0.25,
                                labels = primary.name)
  # secondary_y_label <- data.frame(breaks = mean(c(anno.plot$ymin, secondary.plot$ymax)),
  #                                 labels = secondary.name)
  y_labels <- rbind(y_labels, 
                    primary_y_label, 
                    # secondary_y_label,
                    anno_y_labels)
  
  hlines <- data.frame(yintercept = c(-0.5, max(plot.df$ymax)))
  
  sort.lab <- sort.df %>%
    group_by_(primary$id, primary$label) %>%
    summarise(xmean = mean(c(min(xpos) - 1, xpos)),
              y = length(genes) + 1,
              angle = 90) %>%
    ungroup() %>%
    select_(primary$id, primary$label,"xmean","y","angle")
  names(sort.lab)[1:2] <- c("primary_id","primary_label")
  
  if(showids) {
    sort.lab$primary_label <- paste(sort.lab$primary_id, sort.lab$primary_label)
  }
  
  if(rotatelabels) {
    sort.lab$primary_label <- gsub("[_|;| ]","\n",sort.lab$primary_label)
    sort.lab$angle <- 0
  }
  
  # Segments that divide groups
  segment_lines <- sort.df %>%
    group_by_(primary$id, primary$label) %>%
    summarise(x = max(xpos)) %>%
    mutate(xend = x,
           y = -0.5,
           yend = max(plot.df$ymax)) %>%
    select(x, xend, y , yend) %>%
    as.data.frame()
  
  ##############
  ## Plotting ##
  ##############
  
  p <- ggplot() + 
    # Main heatmap
    geom_rect(data = plot.df,
              aes(xmin = xmin, 
                  xmax = xmax, 
                  ymin = ymin, 
                  ymax = ymax, 
                  fill = fill)) +
    # Axis labels
    scale_x_continuous(breaks = sort.lab$xmean,
                       labels = sort.lab$primary_label,
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = y_labels$breaks,
                       labels = y_labels$labels,
                       expand = c(0, 0)) +
    # fill and theme
    scale_fill_identity(guide=F) +
    theme_classic(base_size=pfontsize) +
    theme(axis.ticks = element_line(size=0.2),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank())
  
  
  # cluster guide lines
  if(showlines) {
    
    p <- p + 
      geom_segment(data = segment_lines,
                   aes(x = x, 
                       xend = xend, 
                       y = y, 
                       yend = yend),
                   size = 0.2) +
      geom_vline(aes(xintercept = 0), 
                 size = 0.2) + 
      geom_hline(data = hlines,
                 aes(yintercept = yintercept), 
                 size = 0.2)
    
  }
  
  
  # Rotate labels checkbox
  if(rotatelabels) {
    p <- p + 
      theme(axis.text.x = element_text(angle = 0, 
                                       hjust = 0.5, 
                                       vjust=1))
  } else {
    p <- p + 
      theme(axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5))
  }
  
  p
  
}


build_legend_plot <- function(data_source, 
                              genes = c("Hspa8","Snap25","Gad2","Vip"), 
                              autorange = "auto", 
                              minrange = 0, 
                              maxrange = 10,
                              scale_type = "scale.log",
                              pfontsize = 14) {
  
  library(dplyr)
  library(ggplot2)
  library(feather)
  
  data <- read_feather(file.path(data_source,"data.feather"), columns = genes)
  
  genes <- sub("-", ".", genes)
  genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
  
  names(data) <- sub("-",".",names(data))
  colors <- colorRampPalette(c("darkblue","white","red"))(1001)
  
  if(autorange == "auto") {
    min.val <- 0
    max.val <- max(unlist(data[, genes]))
  } else if (autorange == "manual") {
    min.val <- as.numeric(minrange)
    max.val <- as.numeric(maxrange)
  }
  
  ## Build geom_rect() compatible table
  legend_data <- data.frame(xmin = 1:1001,
                            xmax = 1:1001+1,
                            ymin = 0,
                            ymax = 1,
                            fill = colors)
  
  if(scale_type == "scale.abs") {
    scale_name <- "RPKM"
  } else if(scale_type == "scale.log") {
    scale_name <- "log10(RPKM + 1)"
    min.val <- log10(min.val + 1)
    max.val <- log10(max.val + 1)
  } else if(scale_type == "scale.rel") {
    scale_name <- "RPKM/max(RPKM)"
    min.val <- min.val/max.val
    max.val <- 1
  } else if(scale_type == "scale.log.rel") {
    scale_name <- "log10(RPKM + 1)/max(log10(RPKM + 1))"
    min.val <- log10(min.val + 1)
    max.val <- log10(max.val + 1)
    min.val <- min.val/max.val
    max.val <- 1
  }
  
  segment_data <- data.frame()
  
  legend_plot <- ggplot(legend_data) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
    geom_segment(aes(x = min(xmin), xend = max(xmax), y = 0, yend = 0)) +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(scale_name, breaks=c(0,250,500,750,1000),
                       labels=round(seq(min.val, max.val, by = (max.val-min.val)/4),2)) +
    theme_classic(base_size = pfontsize) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank())
  
  return(legend_plot)
}  
