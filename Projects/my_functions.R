RMSE_func <- function(y, y_hat) { sqrt(mean((y-y_hat)^2)) / sqrt(mean((y-mean(y))^2)) }
MSE_func <- function(y, y_hat) { mean((y-y_hat)^2) }
R2_func <- function(y, y_hat){ 1 - sum((y - y_hat)^2) / sum((y - mean(y_hat))^2) }

R2adj_func <- function(y, y_hat, k){ 1 - (((1 - R2_func(y, y_hat))*(length(y_hat) - 1)) / (length(y_hat) - k - 1)) }



kable_format <- list(small.mark=",",
                     big.mark=',',
                     decimal.mark='.',
                     nsmall=3,
                     digits=3,
                     scientific=FALSE,
                     big.interval=3L)
theme_TS <- function(base_size=9, base_family="", face="plain"){
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          axis.text=element_text(size=base_size, face=face, family=base_family),
          axis.title=element_text(size=base_size, face=face, family=base_family),
          legend.text=element_text(size=base_size, face=face, family=base_family))
}
qqplot.data <- function(vec, conf = 0.95) {
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  d <- data.frame(resids = vec)
  # confidence intervals ----
  x <- d$resids
  ord <- order(x)
  ord.x <- x[ord]
  n <- dim(d)[1]
  #
  P <- ppoints(n)
  z <- qnorm(P)
  a <- coef(MASS::rlm(ord.x~z))[1]
  b <- coef(MASS::rlm(ord.x~z))[2]
  #
  d$z <- z
  zz<-qnorm(1-(1-conf)/2)
  SE <- (b/dnorm(d$z))*sqrt(P*(1-P)/n)     #[WHY?]
  fit.value <- a+b*d$z
  d$upper <- fit.value+zz*SE
  d$lower <- fit.value-zz*SE



  ggplot(d, aes(sample = resids)) +
    geom_abline(intercept = int, slope = slope, color = "grey50",linetype="dashed") +
    stat_qq() +
    geom_line(data=d, aes(z,upper),color = "grey50",linetype="dashed") +
    geom_line(data=d, aes(z,lower),color = "grey50",linetype="dashed") +
    labs(y="Std. residuals", x="Theoretical quantiles", color="") +
    theme_TS()
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots=length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol=cols, nrow=ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind=TRUE))

      print(plots[[i]], vp=viewport(layout.pos.row=matchidx$row,
                                    layout.pos.col=matchidx$col))
    }
  }
}

res_ana_plot <- function(res,fit, x_labs = "(log) Fitted values") {
  fit_df <- data.frame(t = 1:length(res),
                       res = res,
                       res_sd = (res-mean(res)) / sd(res),
                       res_sd_sqrt = sqrt(abs((res-mean(res)) / sd(res))),
                       fit = fit)



  no_bins <- as.integer(diff(range(fit_df$res)) / (2 * IQR(fit_df$res) / length(fit_df$res)^(1/3)))
  # res vs fit
  multiplot(ggplot(fit_df, aes(x=res)) +
              geom_histogram(bins = no_bins, aes(color = ""), alpha = 1/2) +
              labs(x="Residuals", y="Count",color = "") +
              theme_TS() +
              theme(legend.position="none"),
            ggplot(fit_df) +
              geom_abline(intercept = 0, slope = 0,color = "grey50",linetype="dashed") +
              geom_point(aes(fit,res)) +
              geom_smooth(aes(fit,res), fill=NA, method = "loess", span = 1) +
              labs(x=x_labs, y="Residuals") +
              theme_TS(),
            # scale location
            ggplot(fit_df) +
              geom_point(aes(fit,res_sd_sqrt)) +
              geom_smooth(aes(fit,res_sd_sqrt), fill=NA, method = "loess", span = 1) +
              labs(x=x_labs, y="sqrt(|Std. residuals|)") +
              theme_TS(),
            # qq
            qqplot.data(fit_df$res_sd), cols=2)

}
