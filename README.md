# Auxílio Brasil e voto em Lula em 2022

Reprodução em R do trabalho de Bruno Donadel Mainieri e Enrico Bertolini Prado. [Artigo em PDF](Artigo_Econometria_versao_definitiva.pdf) · [Conversa com a IA](https://econo-chat.vercel.app/conversa).

## Executar no RStudio

Extraia o projeto e abra `Trabalho_Econometria.Rproj`. No console, execute:

```r
source("reproduzir.R", encoding = "UTF-8")
```

Pelo terminal, use `Rscript reproduzir.R` na pasta do projeto. A análise foi testada com R 4.5.0. As rotinas que partem da base municipal usam apenas os pacotes que acompanham o R. O programa usa os dados municipais incluídos, recalcula as tabelas e os gráficos em `resultados/` e compara os números com as 14 tabelas de referência guardadas em `conferencia.rds`. Uma divergência interrompe a conferência e identifica a tabela e a coluna.

## Organização do projeto

O PDF contém o artigo. A base reúne folha, votos, covariáveis, contagens de 2021 e coordenadas das sedes por código IBGE. O arquivo `reproduzir.R` gera as tabelas e figuras. Ele chama `diagnosticos.R` para os testes de influência, dependência espacial e bootstrap. Os insumos da conferência numérica e os sorteios do bootstrap acompanham o projeto. As saídas são criadas em `resultados/` na execução.

As regressões de referência usam 5.570 municípios. Os modelos complementares usam 5.569 municípios em 26 UFs. A unidade da inclinação da folha é o ponto percentual de votos válidos em Lula por dez famílias adicionais na folha por mil habitantes. A tabela de referência inclui intervalos por UF, por região imediata e HC1. Nas comparações centrais, a rotina calcula CR1 por UF com referência t de 25 graus de liberdade.

Para reproduzir os 9.999 sorteios do artigo nas três ponderações, o bootstrap em R lê a matriz de sinais em `sinais_bootstrap.csv.gz`. Cada linha do arquivo contém um sorteio e cada coluna identifica uma UF. A matriz registra os sinais de Rademacher do gerador usado na versão anterior (semente 202211). A rotina recalcula estatísticas e intervalos e verifica uma réplica por reestimação direta. Nos testes de Moran, as 999 permutações são geradas pelo próprio R com semente 202211.

## Localizar as saídas

| Item do artigo | Arquivos em `resultados/` |
| --- | --- |
| Tabela 1 e Figura 1 | `descritivas.csv` e `figura_descritiva.pdf` |
| Tabela 2, colunas (1) a (3) | Linhas com erro `UF` em `referencia.csv` |
| Tabela 2, colunas (4) a (7): folha e Haddad 2018 | Linhas de peso igual em `progressao.csv` (`j` de 2 a 5; Haddad na coluna `gamma`) |
| Figura 2 e contraste entre quartis | `ajustes.pdf` e `contraste_quartis.csv` |
| Tabela 3 (Apêndice A): pesos, NIS, primeiro turno e diferença Lula − Haddad | Linhas com `j` = 3 (e `diferença`, com `j` = 2) em `progressao.csv` |
| Tabela 3 (Apêndice A): amostras alternativas e terços de Haddad | `sensibilidades.csv` (`M` = 3) |
| Tabela 3 (Apêndice A): sem uma UF por vez | `sem_uf.csv` (`M` = 3) |

Os demais arquivos (`diagnosticos_referencia.csv`, `diagnosticos_centrais.csv`, `decis.csv`, `moran.csv`, `wild.csv`, `vif.csv`, `coeficientes_m5.csv`, `estimativas.csv`, `expansao.csv` e as figuras `decis.pdf`, `regioes.pdf` e `intervalos_m3.pdf`) vêm de análises de versões anteriores e não aparecem na versão definitiva do artigo. Eles continuam sendo gerados para que a conferência com `conferencia.rds` permaneça completa.

## Reconstrução opcional com microdados

A base municipal integrada acompanha o projeto. Para reconstruí-la, use `reconstruir_base.R`. Essa rotina adicional requer `data.table`, `jsonlite` e `pdftools`, acesso às fontes e memória para ler os arquivos administrativos. Os nomes e as URLs estão em `fontes.json`. Para baixar os arquivos que faltam e reconstruir, rode `Rscript reconstruir_base.R --baixar`. Se os insumos já estiverem em `brutos/`, use `Rscript reconstruir_base.R`. Os dados reconstruídos ficam em `reconstruidos/`. A rotina verifica a soma da folha, as 2.458 divergências de média monetária, as 17 exceções de nome e os totais dos NIS. Os registros individuais são usados apenas localmente.

