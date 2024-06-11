# GenovisR

`GenovisR` is an R package designed for visualizing and analyzing genotypic data. It provides comprehensive tools to import, visualize, and summarize genotype, haplotype, and dosage data. The package also includes a Shiny application for interactive data exploration.

## Installation

Before installing `GenovisR`, you need to install some prerequisite packages. You can install them using the following commands:

```r
install.packages(c("shiny", "ggplot2", "plotly", "shinyjs"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("Biostrings")
```

To install `GBScleanR`, please visit the [GitHub repository](https://github.com/tomoyukif/GBScleanR)

To install the development version of `GenovisR` from GitHub, you can use the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("ipsr-igb/genovisr", build_vignettes = TRUE)
```

## Vignette

For more information, run the following code on a R console to see a vignette.
```R
browseVignettes(package = "GenovisR")
```


## Features

With `GenovisR`, you can:

- **Create and manage genovis objects**: Generate dummy data for testing or import data directly from GDS files.
- **Visualize genotypic data**: Generate classical graphical genotypes, interactive plots, and histograms of segment lengths.
- **Summarize data**: Compute and visualize statistical summaries of genotypic, haplotypic, and dosage data.
- **Interactive data exploration**: Use the built-in Shiny application for dynamic and interactive exploration of your genotypic data.

## Usage

Below is a quick start guide to get you up and running with `GenovisR`.

### Load the package

```r
library(GenovisR)
```

### Create a genovis object with dummy data

```r
object <- simGenovis()
```

### Import data from a GDS file

```r
object <- gbscleanr2genovis(gds_fn = "/path/to/your/gds_file.gds")
```

### Visualize data

Evaluate segments in the data:

```r
object <- evalSegments(object)
object <- statsGeno(object = object)
```

Generate graphical genotype plots:

```r
plotGraphGeno(object, direction = "h", data = "haplotype", sample = 1:5, width = 0.1)
```

Create histograms of segment lengths:

```r
plotHist(object, data = "haplotype", binwidth = 1e6, fill = "darkgreen")
```

### Interactive Shiny Application

Launch the Shiny application for interactive data exploration:

```r
shinyGenovir()
```

## Contributing

We welcome contributions! If you have any suggestions or find any issues, please feel free to open an issue or submit a pull request on GitHub.

## License

This project is licensed under the GPL-3 License.
```
