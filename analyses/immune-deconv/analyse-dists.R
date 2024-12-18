# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  

  
})





##interested in immune cell composition only, so quantiseq is a more appropriate dataset

option_list <- list(
  make_option(c("--fracs_file"), type = "character",
              help = "fractions file, usually quantiseq output (.rds)"),
  make_option(c("--output_dir"), type = "character", 
              help = "output directory")
)


opt <- parse_args(OptionParser(option_list = option_list))
fracs_file <- opt$fracs_file
output_dir <- opt$output_dir

##for development
#fracs_file <- "OpenPedCan-Project/analyses/immune-deconv/results/quantiseq_output.rds"
#output_dir <- "OpenPedCan-Project/analyses/immune-deconv/results"


## column molecular_subtype
## options 
#MB, WNT
#MB, SHH
#MB, Group3
#MB, Group4


fracs <- readRDS(fracs_file)

fracs <- fracs[fracs$molecular_subtype %in% c("MB, WNT", "MB, SHH", "MB, Group3", "MB, Group4"),]
fracs$molecular_subtype <- gsub("MB, ", "", fracs$molecular_subtype)
#fracs[] <- lapply(fracs, function(x) gsub("MB, ", "", x)) #MB, is redundant 
newFracs <- fracs[, c("Kids_First_Biospecimen_ID","molecular_subtype", "cell_type", "fraction")]

data_wide <- newFracs %>%
  pivot_wider(
    names_from = cell_type,        # Spread immune cell types into columns
    values_from = fraction,        # Use proportion values for the cells
    values_fill = list(fraction = NA)  # Fill any missing values with NA
  )

##above gets the data into the desired format. 

library(ggplot2)
##creates a plot for exploratory analysis 
boxplot <- ggplot(newFracs, aes(x = molecular_subtype, y = fraction, fill = cell_type)) +
  geom_boxplot() +
  labs(title = "Immune Cell Proportions Across Different Medulloblastoma Subtypes", y = "Proportion", x = "Subtype", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

output_bp <- file.path(output_dir, "dists_boxplot.png")
ggsave(output_bp, plot = boxplot, width = 10, height = 8, dpi = 300)


#from the exploratory analysis, uncharacterised cells seem to make up a high proportion of the dataset
#should these be excluded from other visualisations? -- using plotly so this becomes a non issue as you can just zoom in to the desired data
#appears from exploratory analysis that there is a bit of a skew - data is non parametric.
#KW test on individual cell types?? 

#
#kruskal.test(`B cell` ~ molecular_subtype, data = data_wide)

#do a loop on column names that are unique cell types, outputting KW results to a dataframe

cell_types_list <- unique(newFracs$cell_type)
#define dataframe
kruskal_results_df <- data.frame(
  Immune_Cell = character(),
  Kruskal_Wallis_Chi_Squared = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each immune cell column and perform Kruskal-Wallis test
for (cell in cell_types_list) {
  # Perform Kruskal-Wallis test for each immune cell type
  result <- kruskal.test(data_wide[[cell]] ~ data_wide$molecular_subtype)
  
  # Store the result in the dataframe
  kruskal_results_df <- rbind(kruskal_results_df, data.frame(
    Immune_Cell = cell,
    Kruskal_Wallis_Chi_Squared = result$statistic,
    P_Value = result$p.value
  ))
}


##need to write dataframe to an output file without the index
KW_file <- file.path(output_dir, "kruskal_results.csv")
write.table(kruskal_results_df, file = KW_file, sep = ",", row.names = FALSE, quote = FALSE)

##you want to make an interactive stacked bar plot to visualise the outputs. 
##need to calculate the mean across different cancer types - calculate on data wide

data_mean_proportion <- newFracs %>%
  group_by(molecular_subtype, cell_type) %>%
  summarize(Mean_Proportion = mean(fraction), .groups = 'drop')

library(plotly)
plot <- plot_ly(data_mean_proportion, 
                x = ~molecular_subtype, 
                y = ~Mean_Proportion, 
                color = ~cell_type, 
                type = 'bar', 
                text = ~paste(cell_type, ": ", round(Mean_Proportion * 100, 1), "%"),
                hoverinfo = 'text') %>%
  layout(
    title = "Immune Cell Proportions Across Different Medulloblastoma Subtypes",
    barmode = 'stack',
    xaxis = list(title = 'Subtype'),
    yaxis = list(title = 'Proportion'),
    legend = list(title = list(text = 'Immune Cell Type'))
  )

plotly_file <- file.path(output_dir, "immune_cell_dist_bar.html")

htmlwidgets::saveWidget(plot, plotly_file)

