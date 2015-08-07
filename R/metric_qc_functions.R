cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")

# Function to return values of pass or fail depending on whether a vector of
# values meets the specified condition
pass_fail_test <- Vectorize(function(condition) {
    result <- factor(ifelse(condition, "pass", "fail"))
    return(result)
})

outlier_limit <- function(values, direction = c("low", "high")) {
    medVal <- median(values)
    outRange <- 1.5 * IQR(values) / 2
    if (direction == "low") {
        outLim <- medVal - outRange
    } else {
        outLim <- medVal + outRange
    }
    return(outLim)
}

qc_metric <- function(metricsSummary, metric, cutoff) {
    qcTest <- switch(metric, 
                     medianCVcoverage = function(x) { x < cutoff },
                     percentAligned = function(x) { x > cutoff })
    outlierTest <- switch(metric,
                          medianCVcoverage = function(x) { 
                              x > outlier_limit(x, "high") },
                          percentAligned = function(x) {
                              x < outlier_limit(x, "low") })
    direction <- switch(metric, 
                        medianCVcoverage = "high",
                        percentAligned = "low")
    
    metricQC <- metricSummary %>% 
        select(libID, fastqTotalReads,
               metricValue = one_of(metric)) %>% 
        mutate(qcPassFail = pass_fail_test(qcTest(metricValue)),
               metricOutlier = as.factor(ifelse(outlierTest(metricValue),
                                                1, 0)),
               readsOutlier = as.factor(ifelse(fastqTotalReads <
                                                   outlier_limit(fastqTotalReads, "low"),
                                               1, 0)),
               outlierLabel = ifelse(metricOutlier == 1 | 
                                         readsOutlier == 1,
                                     libID, ""))
}

plot_metric_qc <- function(metricsSummary, metric, cutoff) {
    metricQC <- qc_metric(metricQC, metric, cutoff)
    
    if ("fail" %in% levels(metricQC$qcPassFail)) {
        qcCol <- cbPalette[c(2, 6)]
    } else {
        qcCol <- cbPalette[6]
    }
    
    metricQC %>% 
        ggplot(aes(x = fastqTotalReads, y = metricValue)) +
        geom_point(aes(alpha = metricOutlier, shape = readsOutlier), 
                   size = 5, colour = "black") + 
        geom_point(aes(colour = qcPassFail, shape = readsOutlier), 
                   size = 3) +
        scale_colour_manual(values = qcCol) + 
        scale_shape_manual(values = c(16, 15)) +
        scale_alpha_manual(values = c(0, 1)) +
        geom_text(aes(label = outlierLabel),
                  hjust = -0.1, vjust = 1.1, size = 3) +
        theme_classic() +
        ggtitle(metric)
}