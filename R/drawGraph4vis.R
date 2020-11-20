
#' Visualize GO DAG
#'
#' @param Elegans The igraph object
#' @param nodeinfo  The graph structure
#' @param title     The title of the plot
#' @param species   The organism of interest
#' @param ont       The ontology
#' @param val       The important GO categories
#' @param org       The organism
#'
#' @return A graph of the simplified GO tree
#' @keywords internal
#' @importFrom network network network.size %v%<-
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes geom_nodetext
#' @importFrom ggplot2 ggplot scale_fill_manual theme scale_alpha geom_hline scale_color_discrete scale_x_discrete scale_y_continuous guides annotate aes arrow unit  element_blank element_text rel labs element_line
#' @importFrom dplyr distinct %>% rename mutate arrange id
#' @importFrom ggplot2 scale_x_continuous

drawGraph4Vis <- function(Elegans, nodeinfo, title, species, ont, val, org){
  xend <- yend <- level <- label <- NULL
  if(toupper(ont) == "BP"){
    y <- xx.ch
    jump <- 0.03
  }else if(toupper(ont) == "MF"){
    y <- xx.ch1
    jump <- 0.02
  }else if(toupper(ont) == "CC"){
    y <- xx.ch2
    jump <- 0.02
  }


  # Preprocess data
  title <- title
  edges <- unique(as.data.frame(Elegans))

  # Get levels
  levels <- as.data.frame(nodeinfo)
  num <- max(unlist(levels))

  # Count total number of terms

  GOTotal <- lapply(1:length(species), function(x){
    length(species[[as.character(x)]])})
  GOTotal <- unlist(GOTotal)

  # Edges
  colnames(edges) <- c("source", "target")

  source <- edges %>%
    distinct(source) %>%
    rename(id = source)

  target <- edges %>%
    distinct(target) %>%
    rename(id = target)

  # Nodes
  node.data <- data.frame(id=1:num)
  nodes <- node.data %>%
    mutate(label = GOTotal) %>%
    arrange(id)

  #node.data <- data.frame(id=1:num)
  #nodes <- node.data %>%
   # mutate(label = paste("N", id, sep="")) %>%
  #  arrange(id)

  #nodes

  # TAKEN
  valEdges <- c()
  for (i in 1:length(val)) {
    v <- which(val[i] == Elegans[,2])
    valEdges <- c(valEdges, v)
  }
  edges <- edges[valEdges,]


  # Network
  n <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE)

  tags <- c()
  for (i in 1:network.size(n)) {
    tags[i] <- which(levels == i, arr.ind = T)[1]
  }

  n %v% "level" <- tags

  # Generate layout
  start <- 1
  increment <- 1/max(tags)
  traces <- list()
  for (i in 1:network.size(n)) {
    position <- which(levels == i, arr.ind = T)
    if(position[2] == 1){
      traces[[i]] <- c(0, start)
    }
    else if(position[2] == 2){
      traces[[i]] <- c(0.5, start)
    }
    else{
      traces[[i]] <- c(1, start)
    }

    next.pos <- which(levels == i+1, arr.ind = T)

    if(!is.na(next.pos[1])){
      if(next.pos[1] != position[1]){
        start <- start - increment
      }
    }
  }

  layout <- do.call(rbind, traces)

  net <- ggnetwork(n, layout = layout)
  net$level <- factor(net$level, levels=c(1:max(tags)))

  ### category
  #p1 <- rep(1, length(unique(net$vertex.names)))
  #p1[val] <- 2
  num1 <- network.size(n)
  p1 <- sample(0.1, num1, replace = TRUE)
  p1[val] <- 1

  node.data1 <- data.frame(id=1:num1)
  nodes1 <- node.data1 %>%
    mutate(label1 = p1) %>%
    arrange(id)

  # Connection between nodes

  a <- nodeinfo
  connect <- a[,2]
  storeCount <- matrix(0, nrow = dim(a)[1],ncol = dim(a)[2])

  for(x in 1:length(connect)){
    if(connect[x] != 0 && x != dim(a)[1]){
      childTerms <- lapply(species[[as.character(connect[x])]], function(x){ # Finding children terms
        y[[x]]
      })
      childTerms <- unique(unlist(childTerms))
      currentPos <- which(connect[x] == connect)
      nextPos <- a[currentPos + 1,]
      for(i in 1:length(nextPos)){
        if (nextPos[i] != 0){
          storeCount[currentPos,i] <- length(which(species[[as.character(nextPos[i])]] %in% childTerms))
        }
      }
    }
  }

  storeCount <- as.data.frame(storeCount)
  colnames(storeCount) <- c("a","b","c")

  #fadeEdge <- sample(0.1, length(Elegans[,1]), replace = TRUE)
  #valEdges <- c()

  #for (i in 1:length(val)) {
   # v <- which(val[i] == Elegans[,2])
  #  valEdges <- c(valEdges, v)
  #}

  #fadeEdge[n$iel[[53]]] <- 1
  ## the coordinates for annotating text
  ind <- sort(unique(as.integer(net$level)),index.return = TRUE)$ix
  y.ind <- unique(net$y)

  #rearrange "net" data.frame object
  net <- net[order(net[,6]),]
  # Plot network

  p <- ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))  + geom_hline(yintercept = unique(net$y), color = "white", size=0.5) +
    geom_edges(aes(color = level), arrow = arrow(length = unit(3, "pt"), type = "closed"), curvature = 0.2) +
    #geom_nodes(aes(color = as.factor(net$vertex.names), size = 10, stroke = 4)) +
    geom_nodes(aes(color = level), alpha = nodes1$label1, size = 7) +
    #geom_nodes(size = 3, shape = 21, color = "steelblue") + alpha = nodes1$label1
    geom_nodetext(aes(label = label), alpha = nodes1$label1, color= "#3f2a1d", fontface = "bold", size = 2.7)  +
    theme(plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0.4), plot.caption = element_blank(),
          axis.ticks.y =element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_blank(),legend.position = "none", panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank()) + labs(caption = paste("Visualization of the importance", title, "GO categories for reduced", tolower(org) , "DAG", sep = " ")) +
    scale_x_continuous(limits = c(-0.2,2.4),breaks=c(0,0.5,1),labels = c("JN","RN","LN")) +
    guides(size = F) + annotate("text", x = 1.2, y = y.ind[ind] + jump , label = paste("L", 0:(max(tags)-1), ": ", "   "), fontface = 2, hjust = 0) +
    annotate("text", x = 1.4, y = y.ind[ind] + jump, label = paste("J = ", storeCount$a[1:max(tags)], sep = ""), hjust = 0) +
    annotate("text", x = 1.7, y = y.ind[ind] + jump, label = paste("R = ", storeCount$b[1:max(tags)], sep = ""), hjust = 0) +
    annotate("text", x = 2.0, y = y.ind[ind] + jump, label = paste("L = ", storeCount$c[1:max(tags)], sep = ""), hjust = 0)

  listVal <- list("Importance categories" = sort(val))
  #print(listVal)
  return(p)

}






