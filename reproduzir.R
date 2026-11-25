# Abra o projeto no RStudio e execute este arquivo.
options(encoding = "UTF-8", na.action = "na.fail")
dir.create("resultados", showWarnings = FALSE)
source("diagnosticos.R", encoding = "UTF-8")

# Base municipal, equações e inferência ----
ler <- function(nome) read.csv(paste0(nome, ".csv"), check.names = FALSE, fileEncoding = "UTF-8")
salvar <- function(nome, tabela) {
  if (!is.data.frame(tabela)) tabela <- do.call(rbind, tabela)
  write.csv(tabela, file.path("resultados", paste0(nome, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
}
base <- ler("base_municipal")
stopifnot(nrow(base) == 5570, !anyDuplicated(base$codigo_ibge), !anyNA(base),
          sum(base$nis) == 20065780,
          all(base$renda > 0 & base$pop > 0))
base <- transform(base, Y = 100*pt2022, H = 100*pt2018, F = folha/10,
                  primeiro = 100*haddad1, ln_renda = log(renda), ln_pop = log(pop),
                  dif = 100*(pt2022-pt2018), uf = factor(uf))
d <- droplevels(subset(base, uf != "DF"))
stopifnot(nrow(d) == 5569, nlevels(d$uf) == 26)
formulas <- list(Y ~ F, Y ~ F + uf, Y ~ F + uf + H,
                 Y ~ F + uf + H + ln_renda + edu,
                 Y ~ F + uf + H + ln_renda + edu + urb + jovens + idosos + ln_pop)
pesos <- c(igual = "um", população = "pop", votos = "validos2022")
base$um <- 1
d$um <- 1
# CR1 por grupo. Em MQP, o escore de cada linha inclui seu peso.
ajustar <- function(formula, dados, peso = "um", grupo = "uf") {
  dados$peso_ajuste <- dados[[peso]]
  m <- lm(formula, data = dados, weights = peso_ajuste, x = TRUE, y = TRUE, singular.ok = FALSE)
  X <- m$x
  w <- weights(m)
  N <- nrow(X)
  G <- length(unique(dados[[grupo]]))
  A <- summary(m)$cov.unscaled
  A <- A[colnames(X), colnames(X), drop = FALSE]
  S <- rowsum(X * (w * resid(m)), dados[[grupo]])
  correcao <- G/(G-1) * (N-1)/(N-ncol(X))
  V <- crossprod(S %*% A) * correcao
  list(m = m, se = sqrt(pmax(diag(V), 0)), G = G, n = N)
}
linha <- function(a, termo = "F", escala = 1, gl = a$G-1) {
  b <- unname(coef(a$m)[termo]) * escala
  ep <- unname(a$se[termo]) * escala
  margem <- qt(.975, gl) * ep
  data.frame(b, ep, lo = b-margem, hi = b+margem)
}
faixas <- function(x, n) cut(x, quantile(x, seq(0, 1, length.out = n+1)), include.lowest = TRUE, labels = FALSE)

# Resultados de referência ----
vars <- c("Y", "folha", "edu", "renda", "H", "dif")
desc <- lapply(vars, function(v) {
  x <- base[[v]]
  r <- data.frame(v, length(x), mean(x), sd(x), t(quantile(x, c(0,.25,.5,.75,1))), check.names = FALSE)
  names(r) <- c("variavel","count","mean","std","min","25%","50%","75%","max")
  r
})
salvar("descritivas", desc)
refs <- diags <- list()
equacoes <- list(Y ~ F, Y ~ F + edu, Y ~ F + edu + ln_renda)
for (j in seq_along(equacoes)) {
  a <- ajustar(equacoes[[j]], base)
  for (k in names(coef(a$m))) {
    est <- linha(a, k, if (k == "ln_renda") log(1.1) else 1)
    refs[[length(refs) + 1L]] <- data.frame(
      eq = j, termo = if (k == "(Intercept)") "intercepto" else k,
      erro = "UF", est, n = a$n, r2 = summary(a$m)$r.squared
    )
  }
  if (j == 1L) {
    hc1 <- ajustar(equacoes[[j]], base, grupo = "codigo_ibge")
    refs[[length(refs) + 1L]] <- data.frame(
      eq = j, termo = "F", erro = "HC1", linha(hc1, gl = df.residual(hc1$m)),
      n = hc1$n, r2 = summary(hc1$m)$r.squared
    )
  }
  diags[[j]] <- data.frame(eq = j, testar_heterocedasticidade(a$m),
                          acima_100 = sum(fitted(a$m) > 100))
}
salvar("referencia",refs)
salvar("diagnosticos_referencia",diags)
reta <- lm(Y ~ F, data = base)
quad <- ajustar(Y ~ F + I(F^2), base)
decil <- faixas(base$F, 10)
salvar("ajuste_por_decil", data.frame(
  decil = 1:10,
  voto_ajustado = as.vector(tapply(fitted(reta), decil, mean)),
  voto_observado = as.vector(tapply(base$Y, decil, mean))
))
termo_quadratico <- "I(F^2)"
salvar("ajuste_reta", data.frame(
  termo_quadratico = coef(quad$m)[[termo_quadratico]],
  ep_uf = quad$se[[termo_quadratico]],
  ponto_de_maximo_por_mil = -10 * coef(quad$m)[["F"]] / (2 * coef(quad$m)[[termo_quadratico]]),
  municipios_acima_de_350 = sum(base$folha > 350),
  reta_sem_excesso = coef(lm(Y ~ F, data = subset(base, familias <= dom)))[["F"]]
))
q <- quantile(base$F,c(.25,.75))
salvar("contraste_quartis",data.frame(quartil=c("Q1","Q3"),familias_por_mil=10*q,
  reta=predict(reta,data.frame(F=q)),quadratica=predict(quad$m,data.frame(F=q)),
  inclinacao_quadratica=coef(quad$m)[2]+2*coef(quad$m)[3]*q))

# Modelos complementares e sensibilidades ----
registro <- function(j,rotulo,form=formulas[[j]],dados=d,peso="um") {
  a <- ajustar(form,dados,peso)
  est <- linha(a)
  gamma <- if (j>=3) unname(coef(a$m)["H"]) else NA_real_
  data.frame(j,teste=rotulo,est[c("b","lo","hi")],r2=summary(a$m)$r.squared,gamma,n=a$n)
}
prog <- list(registro(2, "igual"))
for (j in 3:4) {
  for (p in names(pesos)) prog[[length(prog) + 1L]] <- registro(j, p, peso = pesos[[p]])
  if (j == 3L) prog[[length(prog) + 1L]] <- registro(j, "NIS", dados = transform(d, F = fav/10))
}
prog <- c(prog, list(
  registro(5, "igual"),
  registro(2, "antecedente", H ~ F + uf),
  registro(2, "diferença", dif ~ F + uf),
  registro(3, "1º turno", dados = transform(d, H = primeiro))
))
salvar("progressao", prog)

# Coeficientes e erros-padrão das colunas (4) a (7) da Tabela 2.
colunas_complementares <- lapply(2:5, function(j) {
  a <- ajustar(formulas[[j]], d)
  termos <- intersect(c("F", "H", "edu", "ln_renda"), names(coef(a$m)))
  do.call(rbind, lapply(termos, function(termo) {
    data.frame(coluna = j + 2L, termo, linha(a, termo), n = a$n,
               r2 = summary(a$m)$r.squared)
  }))
})
salvar("tabela_2_complementar", colunas_complementares)
sens <- function(j,dados,nome,peso="um") {
  a <- ajustar(formulas[[j]],dados,peso)
  data.frame(teste=nome,M=j,n=a$n,G=a$G,linha(a)[c("b","lo","hi")])
}
s <- list(sens(3,subset(d,familias<=dom),"sem excesso de famílias"))
maiores <- head(order(d$pop,d$codigo_ibge,decreasing=TRUE),ceiling(.01*nrow(d)))
for (p in c("igual", "população")) {
  s[[length(s) + 1L]] <- sens(3, d[-maiores, ], paste0("sem 1% maior, ", p), pesos[[p]])
}
# A observação exatamente no segundo corte pertence ao tercil superior no recorte arquivado.
cortes <- quantile(d$H,c(1/3,2/3))
tercil <- ifelse(d$H<=cortes[1],1L,ifelse(d$H>=cortes[2],3L,2L))
for (t in 1:3) s[[length(s)+1]] <- sens(3,d[tercil==t,],paste("tercil",t,"de Haddad"))
salvar("sensibilidades",s)
modelagem <- formulas[[3]]
b_completo <- coef(lm(modelagem, data = d))[["F"]]
sem_uf <- vapply(levels(d$uf), function(u) {
  coef(lm(modelagem, data = d[d$uf != u, ]))[["F"]]
}, numeric(1))
salvar("sem_uf", data.frame(M = 3, uf_omitida = names(sem_uf), b = as.vector(sem_uf),
                           diferenca = as.vector(sem_uf) - b_completo))

# Figuras ----
abrir_pdf <- function(nome,w,h) {
  pdf(file.path("resultados",paste0(nome,".pdf")),width=w,height=h,useDingbats=FALSE)
  par(bty="l",las=1)
}
regioes <- list(Norte=c("AC","AM","AP","PA","RO","RR","TO"),Nordeste=c("AL","BA","CE","MA","PB","PE","PI","RN","SE"),
                "Centro-Oeste"=c("DF","GO","MS","MT"),Sudeste=c("ES","MG","RJ","SP"),Sul=c("PR","RS","SC"))
caixas <- function(var,titulo,unidade) {
  grupos <- lapply(rev(regioes),function(ufs) base[[var]][base$uf %in% ufs])
  boxplot(grupos,horizontal=TRUE,outline=FALSE,col="#bdd6df",main=titulo,xlab=unidade,cex.axis=.8)
}
abrir_pdf("figura_descritiva",9,7)
par(mfrow=c(2,2),mar=c(4,5.6,2.5,.8),cex=.86)
caixas("folha","A. Famílias por região","Famílias na folha por mil habitantes")
caixas("Y","B. Votação por região","Votos válidos em Lula (%)")
plot(base$folha,base$Y,col="#26677940",pch=16,cex=.35,main="C. Associação municipal",
     xlab="Famílias na folha por mil habitantes",ylab="Votos válidos em Lula (%)")
cores <- c("#c6dce5","#8ebcc6","#407f87","#163d54")
plot(base$longitude,base$latitude,asp=1/cos(15*pi/180),axes=FALSE,pch=16,cex=.25,
     col=cores[faixas(base$folha,4)],main="D. Folha por sede municipal",xlab="",ylab="")
legend("bottomleft",legend=paste0("Q",1:4),col=cores,pch=16,bty="n",cex=.8)
dev.off()
abrir_pdf("ajustes",7,4)
plot(base$folha,base$Y,col="#63809b40",pch=16,cex=.35,ylim=range(0,100,fitted(reta),fitted(quad$m)),xlab="Famílias na folha por mil habitantes",ylab="Votos válidos em Lula (%)")
grade <- seq(min(base$F),max(base$F),length.out=400)
lines(10*grade,predict(reta,data.frame(F=grade)),lwd=2)
lines(10*grade,predict(quad$m,data.frame(F=grade)),col="#c65719",lwd=2)
legend("topleft",c("Linear","Quadrático"),col=c("black","#c65719"),lty=1,bty="n")
dev.off()
# Conferência dos resultados ----
referencias <- readRDS("conferencia.rds")
arquivos <- names(referencias)
comparar <- function(arq) {
  atual <- read.csv(file.path("resultados",arq),check.names=FALSE,fileEncoding="UTF-8")
  ref <- referencias[[arq]]
  stopifnot(identical(dim(atual),dim(ref)),identical(names(atual),names(ref)))
  for (v in names(ref)) {
    if (is.numeric(ref[[v]])) {
      if (!isTRUE(all.equal(as.numeric(atual[[v]]),ref[[v]],tolerance=2e-7,check.attributes=FALSE))) {
        stop("Diferença numérica em ",arq,", coluna ",v)
      }
    } else if (!identical(atual[[v]],ref[[v]])) stop("Identificação divergente em ",arq,", coluna ",v)
  }
  message("Conferido: ",arq," (",nrow(ref)," linhas)")
}
invisible(lapply(arquivos,comparar))
