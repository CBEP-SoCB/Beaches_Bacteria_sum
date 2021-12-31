Graphics for Bacteria Levels at Casco Bay Beaches
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
01/23/2021

-   [Introduction](#introduction)
    -   [Handling non-detects](#handling-non-detects)
    -   [Standards](#standards)
        -   [Beaches Program](#beaches-program)
        -   [Maine State Class SB Waters
            Standards](#maine-state-class-sb-waters-standards)
-   [Import Libraries](#import-libraries)
-   [Data Preparation](#data-preparation)
    -   [Initial Folder References](#initial-folder-references)
    -   [Load Data](#load-data)
    -   [Add a “Beach” Identifier](#add-a-beach-identifier)
    -   [Simplify and Correct Beach
        Names](#simplify-and-correct-beach-names)
    -   [Restrict to Recent Conditions](#restrict-to-recent-conditions)
-   [Create Geometric Mean Function](#create-geometric-mean-function)
-   [Jitter Plot of Recent
    Conditions](#jitter-plot-of-recent-conditions)
    -   [Add Geometric Means](#add-geometric-means)
    -   [Add Annotations](#add-annotations)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

We present code for generating draft graphics for the “State of Casco
Bay” report, relating specifically to levels of bacteria observed at
swimming beaches.

The Beaches program monitors bacteria levels (currently) at six Casco
Bay beaches. Data is collected periodically (usually weekly) at each
beach, to inform beach managers and the public about possible risk of
swimming in water that may be polluted by certain pathogens.

The Beaches program measures “enterococci” bacteria, while DMR’s
shellfish program monitor’s “fedal coliform” bacteria. the two measures
are generally correlated, but are not directly comparable because of
different methods.

## Handling non-detects

**All graphics here treat non-detects as equal to the Reporting Limit!**

## Standards

### Beaches Program

104 CFU / 100 ml, for individual observations.

### Maine State Class SB Waters Standards

> the number of enterococcus bacteria in these waters may not exceed a
> geometric mean of 8 CFU per 100 milliliters in any 90-day interval or
> 54 CFU per 100 milliliters in more than 10% of the samples in any
> 90-day interval.

38 M.R.S. §465-B(2)(B)

# Import Libraries

``` r
#library(fitdistrplus)
library(tidyverse)
#> Warning: package 'tidyverse' was built under R version 4.0.5
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.4     v dplyr   1.0.7
#> v tidyr   1.1.3     v stringr 1.4.0
#> v readr   2.0.1     v forcats 0.5.1
#> Warning: package 'ggplot2' was built under R version 4.0.5
#> Warning: package 'tibble' was built under R version 4.0.5
#> Warning: package 'tidyr' was built under R version 4.0.5
#> Warning: package 'readr' was built under R version 4.0.5
#> Warning: package 'dplyr' was built under R version 4.0.5
#> Warning: package 'forcats' was built under R version 4.0.5
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
# library(GGally)

#library(mgcv)      # For GAMs and GAMMs; used here for seasonal smoothers
#library(emmeans)   # For marginal means

#library(mblm)      # for the Thiel-Sen estimators

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())

#library(LCensMeans)
```

# Data Preparation

## Initial Folder References

``` r
sibfldnm    <- 'Derived_Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

dir.create(file.path(getwd(), 'figures'), showWarnings = FALSE)
#dir.create(file.path(getwd(), 'models'),  showWarnings = FALSE)
```

## Load Data

``` r
fn <- "beaches_data.csv"
beach_data <- read_csv(file.path(sibling, fn))
#> Rows: 1895 Columns: 24
#> -- Column specification --------------------------------------------------------
#> Delimiter: ","
#> chr   (8): SiteCode, Sample_ID, Sample_Qualifier, Lab_Qualifier, Weather, Pa...
#> dbl  (11): Year, Month, DOY, Enterococci, Reporting_Limit, Bacteria, Rain24,...
#> lgl   (3): Censored_Flag, Tide_Stage, Current
#> dttm  (1): sdatetime
#> date  (1): sdate
#> 
#> i Use `spec()` to retrieve the full column specification for this data.
#> i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

## Add a “Beach” Identifier

``` r
fn = "beach_locations.csv"
beach_lookup = read_csv(file.path(sibling, fn),
                        col_types = cols(
                          Town = col_character(),
                          Beach_Name = col_character(),
                          SamplePoint = col_character(),
                          Latitude = col_double(),
                          Longitude = col_double()
                        )) %>%
  select(-Latitude, -Longitude)

beach_data <- beach_data %>%
  mutate(Beach = beach_lookup$Beach_Name[match(SiteCode, 
                                               beach_lookup$SamplePoint)])
```

## Simplify and Correct Beach Names

``` r
beach_data <- beach_data %>%
  mutate(Beach = if_else(Beach == 'Stovers Point Preserve', "Stover's Point", Beach)) %>%
  mutate(Beach = if_else(Beach == 'Broad Cove Reserve', "Broad Cove", Beach)) %>%
  mutate(Beach = if_else(Beach == 'Mitchell Field Beach', "Mitchell Field", Beach))
```

## Restrict to Recent Conditions

``` r
recent_data <- beach_data %>%
  filter(Year > 2015)
```

# Create Geometric Mean Function

We pass this to `stat_summary()`.

``` r
gm_mean <- function(x) {
  exp(mean(log(x)))
}
```

# Jitter Plot of Recent Conditions

``` r
jitter_plt <- recent_data %>%
  ggplot(aes(x = Beach, y = Bacteria)) +
  
  geom_jitter(aes(color = Censored_Flag),
              width = 0.3, 
              height = .05,
              alpha = 0.5) +
  
  scale_y_log10() +
  scale_color_manual(values = cbep_colors(), 
                     name = '', labels = c('Observed', 'Below Detection')) +
  
  theme_cbep(base_size = 12) +

  theme(axis.text.x = element_text(angle = 45, size = 9, hjust = 1)) +
  theme(legend.position = c(.6, .9)) +
  
  guides(color = guide_legend(override.aes = list(alpha = c(0.5,0.751) ) )) +

  ylab('Enterococci (MPN / 100ml)') +
  xlab('')
```

### Add Geometric Means

This looks like the geometric means annotation are not properly lined
up, but it is lined up better in the PDF version, which is what counts.

``` r
xanchor <- 3.75
yanchor <- 2200

jitter_plt <- jitter_plt + 
  stat_summary(fun = gm_mean, fill = 'red',shape = 22) 
  
  # annotate('point', x= xanchor, y = yanchor,
  #           size = 3, pch = 22, fill = 'red') +
  # annotate('text', x= xanchor + 0.25, y = yanchor,
  #          hjust = 0, size = 3.5, label = 'Geometric Mean')
```

### Add Annotations

``` r
jitter_plt +
  geom_hline(yintercept = 104, color = 'gray25', lty = 2) +
  #geom_hline(yintercept = 8, color = 'gray25', lty = 2) +
  
  annotate('text', x = 0, y  = 130, label = '104 MPN', 
           size = 2.5, hjust = 0) +
  annotate('text', x = 3.33, y = 2355, label = 'Geometric Mean', 
           size = 3.5, hjust = 0) +
  annotate('point', x = 3.01, y = 2355, shape = 22, fill = 'red')
#> Warning: Removed 6 rows containing missing values (geom_segment).
```

<img src="beaches_graphics_revised_files/figure-gfm/jitter_all-1.png" style="display: block; margin: auto;" />

``` r
ggsave('figures/recent_conditons_jitter_revised.pdf', device = cairo_pdf, 
       width = 5, height = 4)
#> Warning: Removed 6 rows containing missing values (geom_segment).
```
