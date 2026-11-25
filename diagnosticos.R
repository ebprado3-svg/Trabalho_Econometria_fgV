# Testes de heterocedasticidade das três equações de referência (seção 4.3).
testar_heterocedasticidade <- function(modelo) {
  X <- modelo$x
  e2 <- resid(modelo)^2

  # Breusch–Pagan clássico: regressão dos quadrados dos resíduos normalizados.
  bp <- lm.fit(X, e2 / mean(e2))
  bp_lm <- sum((bp$fitted.values - 1)^2) / 2
  bp_gl <- bp$rank - 1L

  # White: regressores originais, quadrados e produtos cruzados.
  regressores <- X[, -1, drop = FALSE]
  produtos <- do.call(cbind, lapply(seq_len(ncol(regressores)), function(k) {
    regressores[, k] * regressores[, k:ncol(regressores), drop = FALSE]
  }))
  white <- lm.fit(cbind(X, produtos), e2)
  white_lm <- length(e2) * (1 - sum(white$residuals^2) / sum((e2 - mean(e2))^2))
  white_gl <- white$rank - 1L

  data.frame(
    bp_lm, bp_gl, bp_p = pchisq(bp_lm, bp_gl, lower.tail = FALSE),
    white_lm, white_gl, white_p = pchisq(white_lm, white_gl, lower.tail = FALSE)
  )
}
