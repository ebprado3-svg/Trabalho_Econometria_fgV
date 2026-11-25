# Auxílio Brasil e voto em Lula em 2022

Reprodução em R do trabalho de Bruno Donadel Mainieri e Enrico Bertolini Prado. [Artigo em PDF](Artigo_Econometria_Vers%C3%A3oFinal.pdf) · [Conversa com a IA](https://econo-chat.vercel.app/conversa).

## Como executar

Abra `Trabalho_Econometria.Rproj` no RStudio e rode no console:

```r
source("reproduzir.R", encoding = "UTF-8")
```

Pelo terminal, use `Rscript reproduzir.R` na pasta do projeto. A reprodução usa a base municipal incluída e os pacotes que acompanham o R. O script carrega os testes de Breusch–Pagan e White de `diagnosticos.R`, gera as saídas em `resultados/` e compara seis tabelas com os valores de referência em `conferencia.rds`. A execução foi verificada com R 4.5.0.

As regressões de referência usam 5.570 municípios. Os modelos complementares excluem o Distrito Federal e usam 5.569 municípios em 26 UFs. A inclinação da folha mede a diferença em pontos percentuais de votos válidos em Lula por dez famílias adicionais na folha por mil habitantes. Os erros-padrão da Tabela 2 são agrupados por UF.

## Onde encontrar os resultados

| Item do artigo | Saída em `resultados/` |
| --- | --- |
| Tabela 1 e Figura 1 | `descritivas.csv` e `figura_descritiva.pdf` |
| Tabela 2, colunas (1) a (3) | Linhas com `erro = UF` em `referencia.csv`. O termo de renda nessa saída já incorpora a variação de 10% usada no texto. |
| Tabela 2, colunas (4) a (7) | `tabela_2_complementar.csv` |
| Seção 4.3, testes de Breusch–Pagan e White | `diagnosticos_referencia.csv`, com a estatística, os graus de liberdade e o valor de p |
| Figura 2 e seção 5.2 | `ajustes.pdf`, `ajuste_reta.csv`, `ajuste_por_decil.csv` e `contraste_quartis.csv` |
| Seções 5.3 e 5.4, comparações entre modelos e ponderações | `progressao.csv`. Os valores `j = 2, 3, 4, 5` correspondem às colunas (4), (5), (6) e (7) da Tabela 2. |
| Tabela 3, Apêndice A | `progressao.csv` para pesos, NIS, primeiro turno e diferença Lula − Haddad. `sensibilidades.csv` para recortes da amostra e terços de Haddad. `sem_uf.csv` para a retirada de uma UF por vez. |

## Reconstrução opcional da base

`reconstruir_base.R` refaz a base municipal a partir das fontes de `fontes.json`. Essa etapa requer `data.table`, `jsonlite` e `pdftools`. Use `Rscript reconstruir_base.R --baixar` para baixar os insumos ou `Rscript reconstruir_base.R` se eles já estiverem em `brutos/`. A base reconstruída fica em `reconstruidos/`.
