# Abra o projeto no RStudio e execute este arquivo.
options(encoding = "UTF-8", na.action = "na.fail")
dir.create("resultados", showWarnings = FALSE)

# Base municipal, equações e inferência ----
ler <- function(nome) read.csv(paste0(nome, ".csv"), check.names = FALSE, fileEncoding = "UTF-8")
salvar <- function(nome, tabela) {
  if (!is.data.frame(tabela)) tabela <- do.call(rbind, tabela)
  write.csv(tabela, file.path("resultados", paste0(nome, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
}
base <- ler("base_municipal")
stopifnot(nrow(base) == 5570, !anyDuplicated(base$codigo_ibge), !anyNA(base),
          sum(base$nis21) == 14654782, sum(base$nis) == 20065780,
          all(base$renda > 0 & base$pop > 0))
base <- transform(base, Y = 100*pt2022, H = 100*pt2018, F = folha/10,
                  primeiro = 100*haddad1, ln_renda = log(renda), ln_pop = log(pop),
                  inicial = 100*nis21/pop, saldo = 1000*(nis-nis21)/pop,
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
  list(m = m, A = A, se = sqrt(pmax(diag(V), 0)), G = G, n = N, correcao = correcao)
}
linha <- function(a, termo = "F", escala = 1, gl = a$G-1) {
  b <- unname(coef(a$m)[termo]) * escala
  ep <- unname(a$se[termo]) * escala
  margem <- qt(.975, gl) * ep
  data.frame(b, ep, lo = b-margem, hi = b+margem)
}
faixas <- function(x, n) cut(x, quantile(x, seq(0, 1, length.out = n+1)), include.lowest = TRUE, labels = FALSE)
p_moran <- function(z, indice) {
  set.seed(202211)
  E <- -1/(length(z)-1)
  obs <- abs(indice(z)-E)
  (1 + sum(replicate(999, abs(indice(sample(z))-E) >= obs)))/1000
}

# Resultados de referência ----
vars <- c("Y", "folha", "edu", "renda", "H", "dif", "saldo")
desc <- lapply(vars, function(v) {
  x <- base[[v]]
  r <- data.frame(v, length(x), mean(x), sd(x), t(quantile(x, c(0,.25,.5,.75,1))), check.names = FALSE)
  names(r) <- c("variavel","count","mean","std","min","25%","50%","75%","max")
  r
})
salvar("descritivas", desc)
refs <- diags <- list()
equacoes <- list(Y ~ F, Y ~ F + edu, Y ~ F + edu + ln_renda)
grupo <- as.integer(factor(base$rgi))
tam <- tabulate(grupo)
indice_rgi <- function(z) {
  soma <- as.vector(rowsum(z, grupo))
  length(z)/sum(tam[grupo]>1) * sum(z*(soma[grupo]-z)/pmax(1,tam[grupo]-1))/sum(z^2)
}
for (j in seq_along(equacoes)) {
  for (tipo in c("UF", "RGI", "HC1")) {
    g <- c(UF="uf", RGI="rgi", HC1="codigo_ibge")[[tipo]]
    a <- ajustar(equacoes[[j]], base, grupo = g)
    m <- a$m
    for (k in names(coef(m))) {
      est <- linha(a, k, if (k=="ln_renda") log(1.1) else 1,
                   if (tipo=="HC1") df.residual(m) else a$G-1)
      refs[[length(refs)+1]] <- data.frame(eq=j, termo=if (k=="(Intercept)") "intercepto" else k,
                                          erro=tipo, est, n=a$n,r2=summary(m)$r.squared)
    }
  }
  z <- resid(m)-mean(resid(m))
  # Breusch–Pagan clássico: metade da soma explicada na regressão
  # dos resíduos quadráticos normalizados pela variância estimada sob H0.
  e2 <- resid(m)^2
  bp_aux <- lm.fit(m$x, e2 / mean(e2))
  bp_lm <- sum((bp_aux$fitted.values - 1)^2) / 2
  bp_gl <- bp_aux$rank - 1L

  # White completo: termos lineares, quadrados e produtos cruzados.
  X <- m$x[, -1, drop = FALSE]
  produtos <- do.call(cbind, lapply(seq_len(ncol(X)), function(k) {
    X[, k] * X[, k:ncol(X), drop = FALSE]
  }))
  white_aux <- lm.fit(cbind(m$x, produtos), e2)
  white_lm <- a$n * (1 - sum(white_aux$residuals^2) / sum((e2 - mean(e2))^2))
  white_gl <- white_aux$rank - 1L
  dec <- tapply(z,faixas(base$F,10),mean)
  covars <- m$x[,-1,drop=FALSE]
  diags[[j]] <- data.frame(eq=j,moran=indice_rgi(z),p_perm=p_moran(z,indice_rgi),
    bp_lm, bp_gl, bp_p=pchisq(bp_lm,bp_gl,lower.tail=FALSE),
    white_lm, white_gl, white_p=pchisq(white_lm,white_gl,lower.tail=FALSE),
    acima_100=sum(fitted(m)>100),
    decis_acima5=sum(abs(dec)>5),max_decil=max(abs(dec)),cook=sum(cooks.distance(m)>4/a$n),
    vif=max(diag(solve(cor(covars)))))
}
salvar("referencia",refs)
salvar("diagnosticos_referencia",diags)
quad <- lm(Y ~ F + I(F^2), base)
q <- quantile(base$F,c(.25,.75))
salvar("contraste_quartis",data.frame(quartil=c("Q1","Q3"),familias_por_mil=10*q,
  reta=predict(lm(Y~F,base),data.frame(F=q)),quadratica=predict(quad,data.frame(F=q)),
  inclinacao_quadratica=coef(quad)[2]+2*coef(quad)[3]*q))

# Modelos complementares e sensibilidades ----
registro <- function(j,rotulo,form=formulas[[j]],dados=d,peso="um") {
  a <- ajustar(form,dados,peso)
  est <- linha(a)
  gamma <- if (j>=3) unname(coef(a$m)["H"]) else NA_real_
  data.frame(j,teste=rotulo,est[c("b","lo","hi")],r2=summary(a$m)$r.squared,gamma,n=a$n)
}
prog <- list()
for (j in 1:5) {
  for (p in names(pesos)) prog[[length(prog)+1]] <- registro(j,p,peso=pesos[[p]])
  prog[[length(prog)+1]] <- registro(j,"NIS",dados=transform(d,F=fav/10))
}
prog <- c(prog,list(registro(2,"antecedente",H~F+uf),registro(2,"diferença",dif~F+uf),
  registro(3,"1º turno",dados=transform(d,H=primeiro)),
  registro(3,"H quadrático",Y~F+uf+H+I((H-mean(H))^2)),
  registro(3,"H por UF",Y~F+uf*H),
  registro(3,"X quadrático",update(formulas[[3]],.~.+I((F-mean(F))^2))),
  registro(5,"X quadrático",update(formulas[[5]],.~.+I((F-mean(F))^2)))))
prog <- do.call(rbind,prog)
salvar("progressao",prog)
a <- ajustar(formulas[[5]],d)
termos <- c("H","ln_renda","edu","urb","jovens","idosos","ln_pop")
fatores <- c(1,log(1.1),1,10,10,10,log(1.1))
rotulos <- c("Haddad 2018 (+1 p.p.)","Renda (+10%)","Superior completo (+1 p.p.)",
             "Urbana (+10 p.p.)","Menores de 18 (+10 p.p.)","65 ou mais (+10 p.p.)","População (+10%)")
contrastes <- lapply(seq_along(termos),function(i) {
  e <- linha(a,termos[i],fatores[i])
  data.frame(contraste=rotulos[i],estimativa_pp=e$b,ic_inf=e$lo,ic_sup=e$hi)
})
salvar("coeficientes_m5",contrastes)
historico <- function(rotulo, form, dados, peso="um", grupo="uf") {
  a <- ajustar(form,dados,peso,grupo)
  data.frame(modelo=rotulo,n=a$n,G=a$G,linha(a))
}
hist <- list()
dem <- c("uf","renda","edu","urb","jovens","idosos","ln_pop")
for (exposicao in c("folha","fav")) {
  dados <- transform(d,F=d[[exposicao]]/10)
  formas <- c(formulas[1:4],list(Y~F+uf+renda+edu,reformulate(c("F",dem),"Y")))
  labs <- c("simples","UF","UF+2018","UF+2018+soc","UF + renda em nível e escolaridade (histórico)",
            "UF + bloco demográfico e renda em nível (histórico)")
  for (j in seq_along(formas)) hist[[length(hist)+1]] <- historico(paste(exposicao,labs[j]),formas[[j]],dados)
  m <- lm(formas[[6]],dados)
  hist[[length(hist)+1]] <- historico(paste(exposicao,"Cook"),formas[[6]],dados[cooks.distance(m)<=4/nrow(dados),])
  for (v in c("jovens","idosos")) hist[[length(hist)+1]] <- historico(paste(exposicao,"sem",v),reformulate(c("F",setdiff(dem,v)),"Y"),dados)
  for (p in c("pop","validos2022")) hist[[length(hist)+1]] <- historico(paste(exposicao,"peso",p),formas[[6]],dados,p)
}
for (g in c("uf","rgi")) {
  hist[[length(hist)+1]] <- historico(paste("referencia",g),equacoes[[1]],base,grupo=g)
  hist[[length(hist)+1]] <- historico(paste("equacao3",g),equacoes[[3]],base,grupo=g)
}
hist[[length(hist)+1]] <- historico("sem excesso",equacoes[[1]],subset(base,familias<=dom),grupo="rgi")
salvar("estimativas",hist)
sens <- function(j,dados,nome,peso="um") {
  a <- ajustar(formulas[[j]],dados,peso)
  data.frame(teste=nome,M=j,n=a$n,G=a$G,linha(a)[c("b","lo","hi")])
}
s <- list(sens(3,subset(d,familias<=dom),"sem excesso de famílias"))
maiores <- head(order(d$pop,d$codigo_ibge,decreasing=TRUE),ceiling(.01*nrow(d)))
for (j in 3:4) for (p in names(pesos)) s[[length(s)+1]] <- sens(j,d[-maiores,],paste0("sem 1% maior, ",p),pesos[[p]])
# A observação exatamente no segundo corte pertence ao tercil superior no recorte arquivado.
cortes <- quantile(d$H,c(1/3,2/3))
tercil <- ifelse(d$H<=cortes[1],1L,ifelse(d$H>=cortes[2],3L,2L))
for (t in 1:3) s[[length(s)+1]] <- sens(3,d[tercil==t,],paste("tercil",t,"de Haddad"))
salvar("sensibilidades",s)
expansao <- list()
for (nome in c("UF","UF e Haddad","UF, Haddad e nível inicial")) {
  vars <- switch(nome,"UF"=c("saldo","uf"),"UF e Haddad"=c("saldo","uf","H"),c("saldo","uf","H","inicial"))
  a <- ajustar(reformulate(vars,"Y"),transform(d,saldo=saldo/10))
  for (v in intersect(vars,c("saldo","inicial"))) {
    e <- linha(a,v)
    expansao[[length(expansao)+1]] <- data.frame(modelo=nome,termo=v,estimativa=e$b,ep=e$ep,lo=e$lo,hi=e$hi,n=a$n,grupos=a$G)
  }
}
salvar("expansao",expansao)

source("diagnosticos.R", encoding = "UTF-8")

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
abrir_pdf("regioes",9,3.5)
par(mfrow=c(1,2),mar=c(4,6,2,.5))
caixas("folha","Famílias por região","Famílias por mil habitantes")
caixas("Y","Votação por região","Votos válidos em Lula (%)")
dev.off()
abrir_pdf("ajustes",7,4)
plot(base$folha,base$Y,col="#63809b40",pch=16,cex=.35,ylim=range(0,100,predict(lm(Y~F,base)),fitted(quad)),xlab="Famílias na folha por mil habitantes",ylab="Votos válidos em Lula (%)")
grade <- seq(min(base$F),max(base$F),length.out=400)
lines(10*grade,predict(lm(Y~F,base),data.frame(F=grade)),lwd=2)
lines(10*grade,predict(quad,data.frame(F=grade)),col="#c65719",lwd=2)
legend("topleft",c("Linear","Quadrático"),col=c("black","#c65719"),lty=1,bty="n")
dev.off()
abrir_pdf("decis",7,3.5)
tab <- do.call(rbind,decis)
curvas <- sapply(c(2,3,5),function(j) tab$residuo_pp[tab$M==j])
matplot(1:10,curvas,type="b",lty=1,pch=1:3,col=1:3,xlab="Decil da folha",ylab="Resíduo médio (p.p.)")
abline(h=0,lty=2)
legend("bottomleft",c("M2","M3","M5"),col=1:3,lty=1,pch=1:3,horiz=TRUE,bty="n")
dev.off()
cr1 <- subset(prog,j==3 & teste=="igual")
boot <- subset(do.call(rbind,intervalos),M==3 & peso=="igual")
abrir_pdf("intervalos_m3",7,2.8)
par(mar=c(4,9,1,.7))
esq <- c(boot$lo,cr1$lo)
dir <- c(boot$hi,cr1$hi)
plot(rep(cr1$b,2),1:2,xlim=range(esq,dir),ylim=c(.6,2.4),yaxt="n",pch=16,
     xlab="Inclinação da folha (p.p. por dez famílias por mil)",ylab="")
axis(2,1:2,c("Bootstrap por UF","CR1 por UF"),tick=FALSE)
segments(esq,1:2,dir,1:2,lwd=2)
abline(v=0,lty=2)
dev.off()

# Conferência das 14 tabelas ----
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
