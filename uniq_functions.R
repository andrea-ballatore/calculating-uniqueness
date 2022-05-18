# 
# Functions to calculate and visualise uniqueness.
# 

#' @return x rescaled between 0 and 1
range01 <- function(x){
  den = (max(x,na.rm = T)-min(x,na.rm = T))
  stopifnot(den != 0)
  (x-min(x,na.rm = T))/den
}

#' @return vector of categorical data from percentages and labels
generate_cat_data_from_pc = function(pcs, labs, size){
  stopifnot(length(pcs)==length(labs))
  x = c()
  for (i in 1:length(pcs)){
    pc = pcs[i]
    lab = labs[i]
    n = size * pc/100
    x = c(x, rep(lab,n))
  }
  sample(x)
}

#' rename scale as z score for clarity
zscore <- function(x){scale(x)}

#' @return a tibble with the given value in each cell
gen_const_matrix = function(n_row, n_col, val){
  tib = as_tibble(matrix(nrow=n_row,ncol=n_col))
  tib[,] = val
  return(tib)
}

#' Generate random data for testing.
#' @return a tibble with random values following a given distribution
gen_random_matrix = function(distrib, n_row, n_col, scale=1, group_split=.2){
  #print(paste('gen_random_matrix',distrib, n_row, n_col, scale))
  stopifnot(distrib %in% c('normal','uniform','2groups','3groups'))
  stopifnot(n_row > 0 && n_col > 0)
  stopifnot(group_split > 0 && group_split < 1)
  if (distrib=='normal'){
    rand = as_tibble(matrix(rnorm(n_row*n_col), nrow=n_row))
    return(rand*scale)
  }
  if (distrib=='uniform'){
    return(as_tibble(matrix(runif(n_row*n_col), nrow=n_row))*scale)
  }
  if (distrib=='2groups'){
    # create data with some common and rare groups
    common_tib = gen_random_matrix('normal', n_row * (1-group_split), n_col, scale)
    rare_tib =   gen_random_matrix('normal', n_row * group_split, n_col, scale*1000)
    tib = bind_rows(common_tib, rare_tib)
    return(tib %>% sample_frac(1))
  }
  if (distrib=='3groups'){
    # create data with some common and rare groups
    common_tib = gen_random_matrix('normal', n_row * .7, n_col, scale)
    med_tib =   gen_random_matrix('normal', n_row * .2, n_col, scale*100)
    rare_tib =   gen_random_matrix('normal', n_row * .1, n_col, scale*1000)
    tib = bind_rows(common_tib, med_tib, rare_tib)
    return(tib %>% sample_frac(1))
  }
}

#' Generate a data frame with a mix of points with gaussian distribution to simulate
#' clustered data
#' @returns a data frame with x and y coordinates.
#' Based on https://www.r-bloggers.com/2018/11/generate-datasets-to-understand-some-clustering-algorithms-behavior/
generate_gaussian_mix <- function(n_clusters, centres_list, sizes_list){
  generate_gaussian_data <- function(n, center, sigma, label) {
    data = rmvnorm(n, mean = center, sigma = sigma)
    data = data.frame(data)
    names(data) = c("x", "y")
    data
  }
  mix = data.frame()
  for (c in seq(n_clusters)){
    n = sizes_list[[c]]
    center = centres_list[[c]]
    sigma = matrix(c(1, 0, 0, 1), nrow = 2)
    gaus = generate_gaussian_data(n, center, sigma, 1)  
    gaus$cluster = c
    mix = bind_rows(mix, gaus)
  }
  mix
}

#' Impute NAs in all numeric columns in tib
impute_numeric_cols_median = function(tib){
  tib %>% mutate_if(is.numeric, 
                function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
}

#' @return TRUE if vector x is unique
is_unique <- function(x){ return(!any(duplicated(x))) }



#' 
OLD_plot_scaled_variable_matrix = function(tib){
  tib_num = tib %>% arrange(-sum_dists) %>% st_drop_geometry() %>% select(-sum_dists, -sum_dists_p_thresh)
  ggplot(dt2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile()
}

#' Plots uniqueness maps using `tmap` package.
plot_uniq_map = function(uniq_geo){
  # table with most unique rows
  #uniq_geo %>% arrange(-sum_dists) %>% head(20) %>% flextable() %>% print()
  scale_bar_sz = .4
  scale_bar_width = .2
  # maps
  tmap_mode('plot')
  # distances
  p = tm_shape(uniq_geo) +
  tm_polygons(col = 'sum_dists',
    legend.hist = TRUE,
    border.alpha = 0,
    frame.lwd = NA, panel.label.bg.color = NA) +
  tm_layout(legend.outside = TRUE, frame = F) +   
  tm_scale_bar(text.size = scale_bar_sz, width=scale_bar_width, position = c("RIGHT", "BOTTOM"))
  print(p)
  
  # distances z values
  p = tm_shape(uniq_geo) +
  tm_polygons(col = 'sum_dists_z',
    legend.hist = TRUE,
    border.alpha = 0,
    frame.lwd = NA, panel.label.bg.color = NA) +
  tm_layout(legend.outside = TRUE, frame = F) + 
  tm_scale_bar(text.size = scale_bar_sz, width=scale_bar_width, position = c("RIGHT", "BOTTOM"))
  print(p)
  
  # classes
  p = tm_shape(uniq_geo) +
  tm_polygons(col = 'sum_dists_p_thresh',
    legend.hist = TRUE,
    #style="cat",
    border.alpha = 0,
    frame.lwd = NA, panel.label.bg.color = NA) +
  tm_layout(legend.outside = TRUE, frame = F) + 
  tm_scale_bar(text.size = scale_bar_sz, width=scale_bar_width, position = c("RIGHT", "BOTTOM"))
  print(p)
}

#' Generate a variable heatmap to show the variability in the variables for the
#' most unique observations.
#' @ show_head_rows: TRUE for top results, FALSE for bottom results
#' @ max_rows: maximum number of rows in the heatmap.
gen_variable_heatmap = function(uniq_tib, show_head_rows=T, max_rows = 30){
  stopifnot(max_rows > 0)
  stopifnot(nrow(uniq_tib) > 0)
  stopifnot("name" %in% names(uniq_tib))
  # filter irrelevant cols out
  uniq_tib = uniq_tib %>% arrange(sum_dists) %>%
    select(-sum_dists_dec,-sum_dists_class,-sum_dists_p_thresh)
  # head or tail of the dataset
  if(show_head_rows)
    uniq_tib = uniq_tib %>% arrange(-sum_dists)
  
  # get order of rows and columns
  row_order = uniq_tib %>% head(max_rows) %>% .[['name']] %>% rev()
  col_order = uniq_tib %>% st_drop_geometry() %>% select_if(is.numeric) %>% names()
  
  # generate long data for heatmap
  tib = uniq_tib %>% st_drop_geometry() %>% head(max_rows) %>%
    select(-sum_dists) %>% column_to_rownames(var='name') %>%
    select_if(is.numeric) %>% rownames_to_column('row_id') %>% 
    pivot_longer(-row_id, names_to = "col_name", values_to = 'value')
  print(tib)
  # fix column names
  tib$col_name = gsub("_"," ",tib$col_name)
  col_order = gsub("_"," ",col_order)
  
  plot_title = paste0("Uniqueness [", uniq_tib$dist_method[[1]], '] ',
                  ifelse(show_head_rows, "(top", "(bottom"), max_rows, "results)")
  
  # plot heatmap
  p = ggplot(tib, aes(x = factor(col_name, level = col_order), 
                  y = factor(row_id, level = row_order),
                  fill = value)) + ggtitle(plot_title) + 
    geom_tile() + theme_bw() + xlab('Scaled variables') + ylab('Observations') + 
    scale_x_discrete(position = "top", labels = wrap_format(10)) +
    geom_text(aes(label=round(value,1)), size=2) +
    colorspace::scale_fill_continuous_divergingx(
      name = "z score",
      mid = 0,
      palette = "RdBu"
    ) + theme(panel.border = element_blank(), panel.grid.major = element_blank())
  
  print(p)
}
