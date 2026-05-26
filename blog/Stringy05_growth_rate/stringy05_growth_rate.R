library(cassowaryr)
library(igraph)

# Compute diameter, length, and ratio

compute_stringy_metrics <- function(x, y) {
  
  sc <- scree(x, y, out.rm = FALSE, binner = NULL)
  mst <- cassowaryr:::gen_mst(sc$del, sc$weights)
  
  diameter_length <- igraph::diameter(mst)
  
  total_length <- sum(igraph::E(mst)$weight)
  
  ratio <- diameter_length / total_length
  
  data.frame(
    diameter = diameter_length,
    length = total_length,
    ratio = ratio
  )
}




run_stringy_simulation <- function(n_start,
                                   n_end,
                                   step,
                                   B = 1000,
                                   seed = 403,
                                   file = "stringy05_results.csv") {
  
  
  n_values <- seq(n_start, n_end, by = step)
  
  all_results <- list()
  counter <- 1
  
  for (n in n_values) {
    
    message("Running n = ", n)
    
    for (sim in seq_len(B)) {
      
      current_seed <- seed + counter
      set.seed(current_seed)
      
      x <- rnorm(n, mean = 0, sd = 1)
      y <- rnorm(n, mean = 0, sd = 1)
      
      metrics <- compute_stringy_metrics(x, y)
      
      all_results[[counter]] <- data.frame(
        n = n,
        simulation_run = sim,
        diameter = metrics$diameter,
        length = metrics$length,
        ratio = metrics$ratio
      )
      
      counter <- counter + 1
    }
  }
  
  results <- do.call(rbind, all_results)
  
  write.csv(results, file = file, row.names = FALSE)
  
  return(results)
}


results <- run_stringy_simulation(
  n_start = 8000,
  n_end = 10000,
  step = 1000,
  B = 1000,
  seed = 403,
  file = "blog/stringy05_growth_rate/stringy05_7000_10000.csv"
)
