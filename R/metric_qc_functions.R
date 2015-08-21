cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")

# Function to return values of pass or fail depending on whether a vector of
# values meets the specified condition
pass_fail_test <- function(value, range) {
    result <- ifelse((value < range[1] | value > range[2]), 
                     "fail", "pass")
    return(result)
}

# Function to calculate upper and lower limits for detecting outliers in a
# vector of values
outlier_limits <- function(values) {
    medVal <- median(values)
    outRange <- 1.5 * IQR(values) / 2
    outLims <- c(medVal - outRange, medVal + outRange)
    return(outLims)
}

# Function to simplify the data frame for downstream functions
format_plot_data <- function(df, metric = c("percentAligned", 
                                            "medianCVcoverage")) {
    df %>% 
        select(libID, x = fastqTotalReads, y = one_of(metric))
}

# Label data points that are outliers or that fail QC
label_libs <- function(libID, xPassFail, yPassFail, 
                       xValue, xOutLims, 
                       yValue, yOutLims) {
    label <- ifelse((xPassFail == "fail" |
                         yPassFail == "fail" |
                         xValue < xOutLims[1] | 
                         xValue > xOutLims[2] |
                         yValue < yOutLims[1] |
                         yValue > yOutLims[2]),
                    libID, "")
}

# Build the final metric plot
plot_metric <- function(df, metric, yRange, xRange) {
    plotDat <- format_plot_data(df, metric)
    xOutLims <- outlier_limits(plotDat$x)
    yOutLims <- outlier_limits(plotDat$y)
    
    plotDat %>% 
        mutate(yPassFail = as.factor(pass_fail_test(y, yRange)) %>% 
                   relevel("pass"),
               xPassFail = as.factor(pass_fail_test(x, xRange)) %>% 
                   relevel("pass"),
               label = label_libs(libID, xPassFail, yPassFail,
                                  x, xOutLims, y, yOutLims),
               label = label_libs(libID, xPassFail, yPassFail,
                                  x, xOutLims, y, yOutLims)) %>% 
        ggplot(aes(x, y)) +
        geom_hline(yintercept = yOutLims, linetype = 3) +
        geom_hline(yintercept = yRange, linetype = 1, colour = "red3", 
                   size = 1) +
        geom_vline(xintercept = xOutLims, linetype = 3) +
        geom_vline(xintercept = xRange, linetype = 1, colour = "red3", 
                   size = 1) +
        geom_point(aes(fill = yPassFail, alpha = xPassFail), 
                   shape = 21, size = 3, colour = "white") +
        geom_text(aes(label = label),
                  hjust = -0.2, vjust = 1.2, size = 3) +
        scale_alpha_manual(values = c(0.8, 0.4), guide = FALSE) +
        scale_fill_colorblind(guide = FALSE) + 
        scale_linetype_discrete(breaks = c(1, 3),
                                labels = c("outlier", "qcCutoff")) +
        xlab("fastqTotalReads") +
        ylab(metric) +
        theme_classic()
}