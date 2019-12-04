
Make_matrix_plot <- function(Mat_data,Set_size_data, Main_bar_data, point_size, line_size, text_scale, labels,
                             shading_data, shade_alpha){
  requireNamespace("UpSetR")
  if(length(text_scale) == 1){
    name_size_scale <- text_scale
  }
  if(length(text_scale) > 1 && length(text_scale) <= 6){
    name_size_scale <- text_scale[5]
  }
  
  Mat_data$line_col <- 'black'
  
  Matrix_plot <- (ggplot()
                  + theme(panel.background = element_rect(fill = "white"),
                          plot.margin=unit(c(-0.2,0.5,0.5,0.5), "lines"),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_text(colour = "gray0",
                                                     size = 7*name_size_scale, hjust = 0.4),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                  + xlab(NULL) + ylab("   ")
                  + scale_y_continuous(breaks = c(1:nrow(Set_size_data)),
                                       limits = c(0.5,(nrow(Set_size_data) +0.5)),
                                       labels = labels, expand = c(0,0)) # changed expand from c(0,0)
                  + scale_x_continuous(limits = c(0,(nrow(Main_bar_data)+1 )), expand = c(0,0))
                  + geom_rect(data = shading_data, aes_string(xmin = "min", xmax = "max",
                                                              ymin = "y_min", ymax = "y_max"),
                              fill = shading_data$shade_color, alpha = shade_alpha)
                  + geom_line(data= Mat_data, aes_string(group = "Intersection", x="x", y="y",
                                                         colour = "line_col"), size = line_size)
                  + geom_point(data= Mat_data, aes_string(x= "x", y= "y"), colour = Mat_data$color,
                               size= point_size, alpha = Mat_data$alpha, shape=16)
                  + scale_color_identity())
  Matrix_plot <- ggplot_gtable(ggplot_build(Matrix_plot))
  return(Matrix_plot)
}

Make_main_bar<-function (Main_bar_data, Q, show_num, ratios, customQ, number_angles, 
          ebar, ylabel, ymax, scale_intersections, text_scale, attribute_plots,margin1scale) {
  bottom_margin <- (-1) * margin1scale
  if (is.null(attribute_plots) == FALSE) {
    bottom_margin <- (-1) * 0.45
  }
  if (length(text_scale) > 1 && length(text_scale) <= 6) {
    y_axis_title_scale <- text_scale[1]
    y_axis_tick_label_scale <- text_scale[2]
    intersection_size_number_scale <- text_scale[6]
  }
  else {
    y_axis_title_scale <- text_scale
    y_axis_tick_label_scale <- text_scale
    intersection_size_number_scale <- text_scale
  }
  if (is.null(Q) == F) {
    inter_data <- Q
    if (nrow(inter_data) != 0) {
      inter_data <- inter_data[order(inter_data$x), ]
    }
    else {
      inter_data <- NULL
    }
  }
  else {
    inter_data <- NULL
  }
  if (is.null(ebar) == F) {
    elem_data <- ebar
    if (nrow(elem_data) != 0) {
      elem_data <- elem_data[order(elem_data$x), ]
    }
    else {
      elem_data <- NULL
    }
  }
  else {
    elem_data <- NULL
  }
  if (is.null(ymax) == T) {
    ten_perc <- ((max(Main_bar_data$freq)) * 0.1)
    ymax <- max(Main_bar_data$freq) + ten_perc
  }
  if (ylabel == "Intersection Size" && scale_intersections != 
      "identity") {
    ylabel <- paste("Intersection Size", paste0("( ", scale_intersections, 
                                                " )"))
  }
  if (scale_intersections == "log2") {
    Main_bar_data$freq <- round(log2(Main_bar_data$freq), 
                                2)
    ymax <- log2(ymax)
  }
  if (scale_intersections == "log10") {
    Main_bar_data$freq <- round(log10(Main_bar_data$freq), 
                                2)
    ymax <- log10(ymax)
  }
  Main_bar_plot <- (ggplot(data = Main_bar_data, aes_string(x = "x", 
                                                            y = "freq")) + 
                      scale_y_continuous(trans = scale_intersections,expand=expand_scale(add=0.1)) + 
                      ylim(0, ymax) + geom_bar(stat = "identity", width = 0.6, 
                                               fill = Main_bar_data$color) + 
                      scale_x_continuous(limits = c(0,(nrow(Main_bar_data) + 1)), expand = c(0, 0), breaks = NULL) + 
                      xlab(NULL) + ylab(ylabel) + labs(title = NULL) + 
                      theme(panel.background = element_rect(fill = "transparent"), 
                            plot.margin = unit(c(0.5, 0.5, bottom_margin, 0.5), "lines"), 
                            panel.border = element_blank(), 
                            axis.title.y = element_text(vjust = -0.8, size = 8.3 * y_axis_title_scale), 
                            axis.text.y = element_text(vjust = 0.3, size = 7 * y_axis_tick_label_scale)))
  if ((show_num == "yes") || (show_num == "Yes")) {
    Main_bar_plot <- (Main_bar_plot + geom_text(aes_string(label = "freq"), 
                                                size = 2.2 * intersection_size_number_scale, vjust = -1, 
                                                angle = number_angles, colour = Main_bar_data$color))
  }
  bInterDat <- NULL
  pInterDat <- NULL
  bCustomDat <- NULL
  pCustomDat <- NULL
  bElemDat <- NULL
  pElemDat <- NULL
  if (is.null(elem_data) == F) {
    bElemDat <- elem_data[which(elem_data$act == T), ]
    bElemDat <- bElemDat[order(bElemDat$x), ]
    pElemDat <- elem_data[which(elem_data$act == F), ]
  }
  if (is.null(inter_data) == F) {
    bInterDat <- inter_data[which(inter_data$act == T), ]
    bInterDat <- bInterDat[order(bInterDat$x), ]
    pInterDat <- inter_data[which(inter_data$act == F), ]
  }
  if (length(customQ) != 0) {
    pCustomDat <- customQ[which(customQ$act == F), ]
    bCustomDat <- customQ[which(customQ$act == T), ]
    bCustomDat <- bCustomDat[order(bCustomDat$x), ]
  }
  if (length(bInterDat) != 0) {
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bInterDat, 
                                              aes_string(x = "x", y = "freq"), fill = bInterDat$color, 
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if (length(bElemDat) != 0) {
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bElemDat, 
                                              aes_string(x = "x", y = "freq"), fill = bElemDat$color, 
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if (length(bCustomDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_bar(data = bCustomDat, 
                                               aes_string(x = "x", y = "freq2"), fill = bCustomDat$color2, 
                                               stat = "identity", position = "identity", width = 0.6))
  }
  if (length(pCustomDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pCustomDat, 
                                                 aes_string(x = "x", y = "freq2"), colour = pCustomDat$color2, 
                                                 size = 2, shape = 17, position = position_jitter(width = 0.2, 
                                                                                                  height = 0.2)))
  }
  if (length(pInterDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pInterDat, 
                                                 aes_string(x = "x", y = "freq"), position = position_jitter(width = 0.2, 
                                                                                                             height = 0.2), colour = pInterDat$color, size = 2, 
                                                 shape = 17))
  }
  if (length(pElemDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pElemDat, 
                                                 aes_string(x = "x", y = "freq"), position = position_jitter(width = 0.2, 
                                                                                                             height = 0.2), colour = pElemDat$color, size = 2, 
                                                 shape = 17))
  }
  Main_bar_plot <- (Main_bar_plot + geom_vline(xintercept = 0, 
                                               color = "gray0") + geom_hline(yintercept = 0, color = "gray0"))
  Main_bar_plot <- ggplotGrob(Main_bar_plot)
  return(Main_bar_plot)
}

Make_size_plot<-function (Set_size_data, sbar_color, ratios, ylabel, scale_sets, 
          text_scale, set_size_angle, set_size.show, set_size.scale_max, 
          set_size.number_size) 
{
  #browser()
  if (length(text_scale) > 1 && length(text_scale) <= 6) {
    x_axis_title_scale <- text_scale[3]
    x_axis_tick_label_scale <- text_scale[4]
  }
  else {
    x_axis_title_scale <- text_scale
    x_axis_tick_label_scale <- text_scale
  }
  if (ylabel == "Set Size" && scale_sets != "identity") {
    ylabel <- paste("Set Size", paste0("( ", scale_sets, 
                                       " )"))
    if (scale_sets == "log2") {
      Set_size_data$y <- log2(Set_size_data$y)
    }
    if (scale_sets == "log10") {
      Set_size_data$y <- log10(Set_size_data$y)
    }
  }
  if (!is.null(set_size.number_size)) {
    num.size <- (set_size.number_size/2.845276) * x_axis_tick_label_scale
  }
  else {
    num.size <- (7/2.845276) * x_axis_tick_label_scale
  }
  Size_plot <- (ggplot(data = Set_size_data, aes_string(x = "x", 
                                                        y = "y")) + 
                  geom_bar(stat = "identity", colour = sbar_color, 
                           width = 0.4, fill = sbar_color, position = "identity") + 
                  scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) + 
                                                        0.5)), breaks = c(0, max(Set_size_data)), expand = expand_scale(add=0.1)) + 
                  theme(panel.background = element_rect(fill = "transparent"),
                        plot.margin = unit(c(-0.11, -1.3, 0.5, 0.5), "lines"), 
                        axis.title.x = element_text(size = 8.3 * x_axis_title_scale), 
                        axis.text.x = element_text(size = 7 * x_axis_tick_label_scale, 
                                                   vjust = 1, hjust = 0.5), axis.line = element_line(colour = "gray0"), 
                        axis.line.y = element_blank(), axis.line.x = element_line(colour = "gray0", size = 0.3),
                        axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
                  xlab(NULL) + ylab(ylabel) + coord_flip())
  if (set_size.show == TRUE) {
    Size_plot <- (Size_plot + geom_text(aes(label = y, vjust = 0.5, 
                                            hjust = 1.2, angle = set_size_angle), size = num.size))
  }
  if (scale_sets == "log10") {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
                                                              0), trans = log10_reverse_trans()))
    }
    else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log10_reverse_trans()))
    }
  }
  else if (scale_sets == "log2") {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
                                                              0), trans = log2_reverse_trans()))
    }
    else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log2_reverse_trans()))
    }
  }
  else {
    if (!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 
                                                              0), trans = "reverse"))
    }
    else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = "reverse"))
    }
  }
  Size_plot <- ggplot_gtable(ggplot_build(Size_plot))
  return(Size_plot)
}

FindSetFreqs<-function (data, start_col, num_sets, set_names, keep_order){
  end_col <- as.numeric(((start_col + num_sets) - 1))
  temp_data <- data[, start_col:end_col]
  temp_data <- temp_data[set_names]
  temp_data <- as.data.frame(colSums(temp_data))
  colnames(temp_data) <- c("y")
  if (keep_order == FALSE) {
    set_order<-order(temp_data$y, decreasing = T)
    temp_data <- temp_data[set_order,]
    set_names<-set_names[set_order]
  }
  else {
    temp_data <- temp_data$y
  }
  x <- seq(1:num_sets)
  temp_data <- cbind(temp_data, x)
  colnames(temp_data) <- c("y", "x")
  return(list(as.data.frame(temp_data),set_names))
}

upset<-function(data, nsets = 5, nintersects = 40, sets = NULL, keep.order = F, set.metadata = NULL, intersections = NULL,
                matrix.color = "gray23", main.bar.color = "gray23", mainbar.y.label = "Intersection Size", mainbar.y.max = NULL,
                sets.bar.color = "gray23",sets.pt.color="gray23", sets.x.label = "Set Size", point.size = 2.2, line.size = 0.7,
                mb.ratio = c(0.70,0.30), expression = NULL, att.pos = NULL, att.color = main.bar.color, order.by = c("freq", "degree"),
                decreasing = c(T, F), show.numbers = "yes", number.angles = 0, group.by = "degree",cutoff = NULL,
                queries = NULL, query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha =0.5,
                empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity",
                scale.sets = "identity", text.scale = 1, set_size.angles = 0 , set_size.show = FALSE, set_size.numbers_size = NULL, set_size.scale_max = NULL,
                margin1scale=0.65){
  requireNamespace("UpSetR")
  #browser()
  startend <- UpSetR:::FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]
  
  if(color.pal == 1){
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF")
  } else{
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                 "#CC79A7")
  }
  
  if(is.null(intersections) == F){
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::specific_intersections(data, first.col, last.col, intersections, order.by, group.by, decreasing,
                                                 cutoff, main.bar.color, Set_names)
  } else if(is.null(intersections) == T){
    Set_names <- sets
    if(is.null(Set_names) == T || length(Set_names) == 0 ){
      Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, nsets)
    }
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, Set_names, nintersects, main.bar.color,
                                  order.by, group.by, cutoff, empty.intersections, decreasing)
  }
  
  Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
  labels <- UpSetR:::Make_labels(Matrix_setup)
  #Chose NA to represent NULL case as result of NA being inserted when at least one contained both x and y
  #i.e. if one custom plot had both x and y, and others had only x, the y's for the other plots were NA
  #if I decided to make the NULL case (all x and no y, or vice versa), there would have been alot more if/else statements
  #NA can be indexed so that we still get the non NA y aesthetics on correct plot. NULL cant be indexed.
  att.x <- c(); att.y <- c();
  if(is.null(attribute.plots) == F){
    for(i in seq_along(attribute.plots$plots)){
      if(length(attribute.plots$plots[[i]]$x) != 0){
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if(length(attribute.plots$plots[[i]]$x) == 0){
        att.x[i] <- NA
      }
      if(length(attribute.plots$plots[[i]]$y) != 0){
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if(length(attribute.plots$plots[[i]]$y) == 0){
        att.y[i] <- NA
      }
    }
  }
  
  BoxPlots <- NULL
  if(is.null(boxplot.summary) == F){
    BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, Set_names)
    BoxPlots <- list()
    for(i in seq_along(boxplot.summary)){
      BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], att.color)
    }
  }
  
  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if(is.null(queries) == F){
    custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
    customDat <- UpSetR:::customQueries(New_data, custom.queries, Set_names)
    legend <- UpSetR:::GuideGenerator(queries, palette)
    legend <- UpSetR:::Make_legend(legend)
    if(is.null(att.x) == F && is.null(customDat) == F){
      customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
    }
    customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, All_Freqs, custom.queries)
  }
  if(is.null(queries) == F){
    Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
    Matrix_col <- UpSetR:::intersects(QuerieInterData, Intersection, New_data, first.col, Num_of_set,
                                      All_Freqs, expression, Set_names, palette)
    Element <- UpSetR:::SeperateQueries(queries, 1, palette)
    EBar_data <-UpSetR:::ElemBarDat(Element, New_data, first.col, expression, Set_names,palette, All_Freqs)
  } else{
    Matrix_col <- NULL
  }
  Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, Matrix_col, matrix.dot.alpha)
  
 
  Set_sizes <- FindSetFreqs(New_data, first.col, Num_of_set, Set_names, keep.order)
  Set_names<-Set_sizes[[2]]
  Set_sizes<-Set_sizes[[1]]
  
  ###################### this was added ########################
  # originally from https://www.r-bloggers.com/hacking-our-way-through-upsetr/
  # and then SPF modified further
  sets.pt.color<-sets.pt.color[Set_names]
  sets.bar.color<-sets.bar.color[Set_names]
  labels<-labels[labels %in% Set_names]
  if(length(sets.pt.color) > 1){
    for(i in 1:length(sets.pt.color)) {
      j <- which(Matrix_layout$y == i & Matrix_layout$value == 1)
      if(length(j) > 0) Matrix_layout$color[j] <- sets.pt.color[i]
    }
  }
  ### end modified
  
  
  Bar_Q <- NULL
  if(is.null(queries) == F){
    Bar_Q <- UpSetR:::intersects(QuerieInterBar, Intersection, New_data, first.col, Num_of_set, All_Freqs, expression, Set_names, palette)
  }
  QInter_att_data <- NULL
  QElem_att_data <- NULL
  if((is.null(queries) == F) & (is.null(att.x) == F)){
    QInter_att_data <- UpSetR:::intersects(QuerieInterAtt, Intersection, New_data, first.col, Num_of_set, att.x, att.y,
                                           expression, Set_names, palette)
    QElem_att_data <- UpSetR:::elements(QuerieElemAtt, Element, New_data, first.col, expression, Set_names, att.x, att.y,
                                        palette)
  }
  AllQueryData <- UpSetR:::combineQueriesData(QInter_att_data, QElem_att_data, customAttDat, att.x, att.y)
  
  ShadingData <- NULL
  
  if(is.null(set.metadata) == F){
    ShadingData <- get_shade_groups(set.metadata, Set_names, Matrix_layout, shade.alpha)
    output <- Make_set_metadata_plot(set.metadata, Set_names)
    set.metadata.plots <- output[[1]]
    set.metadata <- output[[2]]
    
    if(is.null(ShadingData) == FALSE){
      shade.alpha <- unique(ShadingData$alpha)
    }
  } else {
    set.metadata.plots <- NULL
  }
  if(is.null(ShadingData) == TRUE){
    ShadingData <- UpSetR:::MakeShading(Matrix_layout, shade.color)
  }
  
  
  Main_bar <- suppressMessages(Make_main_bar(All_Freqs, Bar_Q, show.numbers, mb.ratio, customQBar, number.angles, EBar_data, mainbar.y.label,
                                                      mainbar.y.max, scale.intersections, text.scale, attribute.plots,margin1scale=margin1scale))
  Matrix <- Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, point.size, line.size,
                             text.scale, labels, ShadingData, shade.alpha)
  Sizes <- Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, sets.x.label, scale.sets, text.scale, set_size.angles,set_size.show,
                                   set_size.scale_max, set_size.numbers_size)
  
  # Make_base_plot(Main_bar, Matrix, Sizes, labels, mb.ratio, att.x, att.y, New_data,
  #                expression, att.pos, first.col, att.color, AllQueryData, attribute.plots,
  #                legend, query.legend, BoxPlots, Set_names, set.metadata, set.metadata.plots)
  
  structure(class = "upset",
            .Data=list(
              Main_bar = Main_bar,
              Matrix = Matrix,
              Sizes = Sizes,
              labels = labels,
              mb.ratio = mb.ratio,
              att.x = att.x,
              att.y = att.y,
              New_data = New_data,
              expression = expression,
              att.pos = att.pos,
              first.col = first.col,
              att.color = att.color,
              AllQueryData = AllQueryData,
              attribute.plots = attribute.plots,
              legend = legend,
              query.legend = query.legend,
              BoxPlots = BoxPlots,
              Set_names = Set_names,
              set.metadata = set.metadata,
              set.metadata.plots = set.metadata.plots)
  )
}
