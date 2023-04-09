scan1 <- function(y      = y,
                  geno_a = geno_a,
                  geno_d = geno_d,
                  tname  = tname,
                  ename  = ename,
                  qq     = qq,
                  LOD    = LOD)
{
  y_na <- which(is.na(y))
  if (length(y_na) != 0) {
    y <- y[-y_na]
    aa <- geno_a[, -y_na]
    dd <- geno_d[, -y_na]
  } else {
    aa <- geno_a
    dd <- geno_d
  }
  if (qq) y <- stats::qqnorm(y, plot.it = FALSE)$x
  y <- matrix(y, ncol = 1)
  qa <- unique(aa)
  all_zero_id <- as.numeric(which(apply(dd, 1, function(x){all(x == 0)})))
  if (length(all_zero_id) != 0) dd <- dd[-all_zero_id, ]
  qd <- unique(dd)
  ad <- rbind(qa, qd)
  ma <- nrow(qa)
  suppressWarnings(lasso <- glmnet::glmnet(x                = t(ad),
                                           y                = y,
                                           family           = "gaussian",
                                           alpha            = 1,
                                           standardize      = TRUE,
                                           intercept        = TRUE,
                                           pmax             = 100))
  best_id <- length(lasso$lambda)
  alpha <- lasso$a0[best_id]
  beta <- as.numeric(lasso$beta[, best_id])
  not_zero <- which(beta != 0)
  beta <- matrix(beta[not_zero], ncol = 1)
  deltat_y <- y - alpha - t(ad[not_zero, , drop = FALSE]) %*% beta
  lod <- sapply(1:length(not_zero), function(i){
    # i = 1
    id = not_zero[i]
    geno1 <- ad[id, ]
    deltat_y_1 <- deltat_y + geno1 * beta[i]
    # L0
    mu0 <- mean(deltat_y_1)
    s0  <- mean((deltat_y_1 - mu0)^2)
    l0 <- LogLikelihood(y = deltat_y_1, expec = mu0, sigma2 = s0)
    code <- sort(unique(geno1))
    if (length(intersect(code, c(-1, 0, 1))) == 3) {
      # L1
      id1 <- which(geno1 == -1)
      id2 <- which(geno1 ==  1)
      id3 <- which(geno1 ==  0)
      mu1 <- mean(deltat_y_1[id1])
      mu2 <- mean(deltat_y_1[id2])
      mu3 <- mean(deltat_y_1[id3])
      s11 <- sum((deltat_y_1[id1] - mu1)^2)
      s12 <- sum((deltat_y_1[id2] - mu2)^2)
      s13 <- sum((deltat_y_1[id3] - mu3)^2)
      s1 <- (s11 + s12 + s13) / length(y)
      l1 <- LogLikelihood(y = deltat_y_1[id1], expec = mu1, sigma2 = s1) +
        LogLikelihood(y = deltat_y_1[id2], expec = mu2, sigma2 = s1) +
        LogLikelihood(y = deltat_y_1[id3], expec = mu3, sigma2 = s1)
    } else {
      if (length(intersect(code, c(-1, 1))) == 2) {
        id1 <- which(geno1 == -1)
        id2 <- which(geno1 ==  1)
      } else if (length(intersect(code, c(-1, 0))) == 2) {
        id1 <- which(geno1 == -1)
        id2 <- which(geno1 ==  0)
      } else if (length(intersect(code, c(0, 1))) == 2) {
        id1 <- which(geno1 == 0)
        id2 <- which(geno1 == 1)
      }
      # L1
      mu1 <- mean(deltat_y_1[id1])
      mu2 <- mean(deltat_y_1[id2])
      s11 <- sum((deltat_y_1[id1] - mu1)^2)
      s12 <- sum((deltat_y_1[id2] - mu2)^2)
      s1 <- (s11 + s12) / length(y)
      l1 <- LogLikelihood(y = deltat_y_1[id1], expec = mu1, sigma2 = s1) +
        LogLikelihood(y = deltat_y_1[id2], expec = mu2, sigma2 = s1)
    }
    return(-2 * (l0 - l1) / (2 * log(10)))
  })
  if (any(lod > LOD)) {
    not_zero <- not_zero[lod > LOD]
    beta <- beta[lod > LOD]
    lod <- lod[lod > LOD]
    if (max(not_zero) <= ma) {
      # all a
      marker <- rownames(qa)[not_zero]
      type <- rep("additive", length(not_zero))
    } else if (min(not_zero) > ma) {
      # all d
      marker <- rownames(qd)[(not_zero - ma)]
      type <- rep("dominant", length(not_zero))
    } else {
      marker <- c(rownames(qa)[not_zero[not_zero <= ma]],
                  rownames(qd)[not_zero[not_zero >  ma] - ma])
      type <- c(rep("additive", length(which(not_zero <= ma))),
                rep("dominant", length(which(not_zero >  ma))))
    }
    result <- data.frame(Trait       = tname,
                         Environment = ename,
                         Marker      = marker,
                         Effect      = beta,
                         Type        = type,
                         LOD         = lod,
                         Same_marker = NA)
    result$Same_marker <- sapply(1:nrow(result), function(i){
      if (result$Type[i] == "additive") {
        marker_a <- result$Marker[i]
        id <- rownames(aa)[as.numeric(which(apply(aa, 1, function(x){
          all(x == aa[rownames(aa) == marker_a, ])
        })))]
        if (length(id) > 1) {
          sm <- paste0(setdiff(id, marker_a), collapse = ",")
        } else {
          sm <- NA
        }
      } else if (result$Type[i] == "dominant") {
        marker_d <- result$Marker[i]
        id <- rownames(dd)[as.numeric(which(apply(dd, 1, function(x){
          all(x == dd[rownames(dd) == marker_d, ])
        })))]
        if (length(id) > 1) {
          sm <- paste0(setdiff(id, marker_d), collapse = ",")
        } else {
          sm <- NA
        }
      }
      return(sm)
    })
    var_qtl <- sapply(1:length(not_zero), function(i){
      id = not_zero[i]
      yi <- ad[id, ] * beta[i]
      return(stats::var(yi))
    })
    ye <- y - alpha - t(ad[not_zero, , drop = FALSE]) %*% matrix(beta, ncol = 1)
    var_error <- stats::var(as.numeric(ye))
    result$`r^2(%)` <- round(var_qtl / (var_error + sum(var_qtl)) * 100, 4)
  } else {
    result <- NULL
  }
  return(result)
}
