library(readxl)
library(plotly)
library(dplyr)


execel_1 <- readxl::read_xlsx("/home/kevhu/data/imaging/2020_03_27_KevinHu_Grid.xlsx", sheet = 2, col_names = FALSE)

extract_vars <- function(df){
  tc_df <- matrix(rep(0, nrow(execel_1) * ncol(execel_1)),nrow = nrow(execel_1), ncol = ncol(execel_1))
  cell_df <- matrix(rep(0, nrow(execel_1) * ncol(execel_1)),nrow = nrow(execel_1), ncol = ncol(execel_1))
  edge_df <- matrix(rep(0, nrow(execel_1) * ncol(execel_1)),nrow = nrow(execel_1), ncol = ncol(execel_1))
  for (i in 1:nrow(execel_1)) {
    for (j in 1:ncol(execel_1)) {
      tmpVar <- unlist(execel_1[i,j])
      res <- strsplit(tmpVar, ",")
      #tmpTc <- as.numeric(unlist(res)[1])
      #tmpCell <- as.numeric(unlist(res)[2])
      #tmpEdge <- as.numeric(unlist(res)[3])
      #ifelse(tmpTc == 0, 
      #       tmpRow_tc <- list(i, j, c(0,0,0,0,0)),
      #       tmpRow_tc <- list(i,j, c(tmpTc - 0.02, tmpTc - 0.01, tmpTc, tmpTc + 0.01, tmpTc + 0.02)))
      
      #ifelse(tmpCell == 0, 
      #       tmpRow_cel <- list(i, j, c(0,0,0,0,0)),
      #       tmpRow_cel <- list(i,j, c(tmpCell - 0.02, tmpCell - 0.01, tmpCell, tmpCell + 0.01, tmpCell + 0.02)))
      
      #ifelse(tmpEdge == 0, 
      #       tmpRow_edge <- list(i, j, c(0,0,0,0,0)),
      #       tmpRow_edge <- list(i,j, c(tmpEdge - 0.02, tmpEdge - 0.01, tmpEdge, tmpEdge + 0.01, tmpEdge + 0.02)))
      
      #tmpRow_tc <- c(i, j, as.numeric(unlist(res)[1]))
      #tmpRow_cel <- c(i, j , as.numeric(unlist(res)[2]))
      #tmpRow_edge <- c(i, j , as.numeric(unlist(res)[3]))
      
      
      tc_df[i,j] <- as.numeric(unlist(res)[1])
      cell_df[i,j] <- as.numeric(unlist(res)[2])
      edge_df[i,j] <- as.numeric(unlist(res)[3])
      
      #tc_df <- rbind(tc_df, tmpRow_tc)
      #cell_df <- rbind(cell_df, tmpRow_cel)
      #edge_df <- rbind(edge_df, tmpRow_edge)
    }
  }
  #colnames(tc_df) <- c("row", "column", "tumor_content")
  #rownames(tc_df) <- NULL
  #colnames(cell_df) <- c("row", "column", "tumor_cellularity")
  #rownames(cell_df) <- NULL
  #colnames(edge_df) <- c("row", "column", "edge")
  #rownames(edge_df) <- NULL
  res <- list(tc_df, cell_df, edge_df)
}

vars_1 <- extract_vars(execel_1)

x <- 1:nrow(execel_1)
y <- 1:ncol(execel_1)
z <- vars_1[[1]]/100
z2 <- vars_1[[2]]/100

plot_ly() %>% 
  add_surface(x = ~x, y = ~y, z = ~z, opacity = 0.85,
                          colorscale = list(c(0,"rgb(107,184,214)"),c(1,"rgb(0,90,124)"))) %>%
  add_surface(x = ~x, y = ~y, z = ~z2, opacity = 0.85,
              colorscale = list(c(0,"rgb(169, 38, 38)"),c(1,"rgb(255, 146, 146)")))





### trying to create data necessary to use a linear optimiztion
### variables are the base ones + others: (base) cellularity, tumor content, edge (window-based) tumor-content avg, cellularity avg, avg edge(?)
### window based is used to look at surrounding patches



### need to first pad ends with so I windowing can work properly




