#'@export
boots <- function(X, Y) {
    # generates B bootstrap samples
    B = 1
    bootstrap.samples <- matrix(sample(1:nrow(X), nrow(X)*B, replace=T), nrow=B)

    # to hold results (computes AUC
    aucs <- matrix(rep(0, B), nrow=B)

    # for each bootstrap sample, calculate AUC of naive Bayes
    for (bs in 1:B) {
      print(c(bs, B))
      fault <- F
      b.sample <- bootstrap.samples[bs,]

      # which data pts aren't in the bootstrap sample? this will be the test set
      oob <- which(!(1:nrow(X) %in% b.sample))

      X.b <- X[b.sample, ]
      Y.b <- Y[b.sample]

      nb <- naiveBayes(X.b, Y.b, laplace=1)

      X.oob <- X[oob, ]
      Y.oob <- Y[oob]

      # AUC is poorly defined when there are no "true" labels to predict
      if (sum(as.logical(Y.oob)) < 1) next

      pred <- predict(nb, X.oob, type="raw")[,"TRUE"]
      roc.pred <- roc(pred, Y.oob)
      aucs[bs] <- auc(roc.pred)#, min = 0, max = 1)
    }
    return(aucs)
}



#' @export
syntheticdata <- function(nsamples=2000, class=F, noisevars=0) {
  tmp <- cbind(matrix(runif(nsamples*6, 0, 1), ncol=6),
               matrix(runif(nsamples*4, 0.6, 1), ncol=4))
  colnames(tmp) <- c("x1", "x2", "x3", "x6", "x7", "x9", "x4", "x5", "x8", "x10")

  if (noisevars > 0) {
    noise <- matrix(rnorm(nsamples*noisevars), ncol=noisevars)
    mu <- runif(noisevars, -1, 1)
    noise <- t(t(noise)+mu)
  } else {
    noise <- NULL
  }

  #x6 is a noise variable in this thing
  response <- apply(tmp, 1,
                    function(x) { pi^(x["x1"]*x["x2"])*sqrt(2*x["x3"])      # 1,2,3,
                      - asin(x["x4"])                                       # 4
                      + log(x["x3"] + x["x5"])                              # 3, 5
                      - (x["x9"]/x["x10"])/sqrt(x["x7"]/x["x8"])            # 7, 8, 9, 10
                      - x["x2"]*x["x7"]                                     # 2, 7
                    }
  )

  m_resp <- -0.25
  if (class) {
    response <- as.factor(response > m_resp)
  }

  return(data.frame(tmp, noise, response))
}

#'@export
predictors <- function(dataset){
  return(as.matrix(dataset[, -ncol(dataset)]))
}

#'@export
responses <- function(dataset){
  return(dataset$response)
}
