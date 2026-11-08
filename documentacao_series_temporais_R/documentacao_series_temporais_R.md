# Séries Temporais: do conceito ao primeiro modelo

**Um guia acessível — porém rigoroso — para iniciantes, com exemplos em R**

> Documento produzido em 05/11/2026. Fontes principais: *Forecasting: Principles and Practice* (Hyndman & Athanasopoulos, 3ª ed., [otexts.com/fpp3](https://otexts.com/fpp3/)), documentação oficial dos pacotes R (*forecast*, *fable/fabletools*, *prophet*, *rbcb*, *GetBCBData*), materiais didáticos da Duke University (Robert Nau), NCSS, ABS (agência de estatísticas da Austrália), além das séries oficiais do Banco Central do Brasil (BCB) consultadas diretamente. Todos os links estão na seção 15.

---

## Sumário

1. [O que é uma série temporal](#1-o-que-é-uma-série-temporal)
2. [Componentes e decomposição](#2-componentes-e-decomposição)
3. [Estacionariedade e transformações](#3-estacionariedade-e-transformações)
4. [Autocorrelação: ACF e PACF](#4-autocorrelação-acf-e-pacf)
5. [Modelos da família ARIMA e a metodologia Box–Jenkins](#5-modelos-da-família-arima-e-a-metodologia-boxjenkins)
6. [Previsão e validação (backtesting)](#6-previsão-e-validação-backtesting)
7. [Métricas de acurácia: como medir se o modelo é bom](#7-métricas-de-acurácia-como-medir-se-o-modelo-é-bom)
8. [Incerteza: intervalos de previsão e comunicação](#8-incerteza-intervalos-de-previsão-e-comunicação)
9. [Armadilhas, limitações e controvérsias](#9-armadilhas-limitações-e-controvérsias)
10. [Exemplos completos em R](#10-exemplos-completos-em-r)
11. [Estudo de caso: inflação e juros no Brasil](#11-estudo-de-caso-inflação-e-juros-no-brasil)
12. [Glossário](#12-glossário)
13. [Como continuar aprendendo](#13-como-continuar-aprendendo)
14. [Limitações explícitas deste documento](#14-limitações-explícitas-deste-documento)
15. [Referências](#15-referências)

---

## 1. O que é uma série temporal

### 1.1 Definição

**Série temporal** é um conjunto de observações de uma mesma variável, **ordenadas no tempo**. Segundo a definição usada pela agência de estatísticas da Austrália ([ABS – Time Series Analysis: The Basics](https://www.abs.gov.au/websitedbs/d3310114.nsf/home/time+series+analysis:+the+basics)), trata-se de “uma coleção de observações de itens bem definidos, obtidas por medições repetidas ao longo do tempo” — por exemplo, o valor das vendas do varejo a cada mês. Dados coletados apenas uma vez, ou sem uma estrutura temporal clara, **não** formam uma série temporal.

O que torna esse tipo de dado especial é que **a ordem importa**. Em um experimento clássico, embaralhar as linhas de uma tabela não muda nada; em uma série temporal, embaralhar destrói justamente a informação que queremos estudar: a dependência entre o valor de hoje e os de ontem, a sazonalidade de ano a ano, a tendência de longo prazo. É por isso que as técnicas de séries temporais têm regras próprias — inclusive para fazer validação de modelos (seção 6).

### 1.2 Exemplos

- Preço de fechamento do dólar a cada dia útil;
- IPCA (inflação) mês a mês;
- Número de passageiros transportados por ônibus a cada hora;
- Temperatura média diária de uma cidade;
- Focos de queimada por semana em cada município.

Repare que alguns desses exemplos são **diários**, outros **mensais**, **horários** etc. O intervalo entre observações é chamado de **período** ou **frequência** da série (m = 12 para dados mensais, m = 4 para trimestrais, m = 7 para padrões semanais em dados diários, e assim por diante).

### 1.3 Classificações úteis

Algumas distinções aparecem em praticamente todos os livros da área:

| Classificação | Significado | Exemplo |
|---|---|---|
| **Estoque × Fluxo** | Estoque mede um “instantâneo” numa data; fluxo mede atividade acumulada num período | População empregada (estoque) × produção de fábricas no mês (fluxo) |
| **Determinística × Estocástica** | Se o futuro é regido por uma fórmula exata ou por distribuições de probabilidade | Série calculada × série econômica real |
| **Regular × Irregular** | Observações em intervalos fixos ou não | Mensal fixo × ocorrências de eventos em datas aleatórias |

Uma consequência prática importante: **séries de fluxo** sofrem efeitos de calendário (número de dias úteis diferente a cada mês, feriados móveis como Páscoa), que precisam ser tratados antes de analisar sazonalidade ([ABS](https://www.abs.gov.au/websitedbs/d3310114.nsf/home/time+series+analysis:+the+basics)).

### 1.4 Por que analisar séries temporais?

Os objetivos típicos são três (não necessariamente simultâneos):

1. **Descrever** o comportamento da série (tendência? sazonalidade? ciclos?);
2. **Explicar/entender** a estrutura de dependência (autocorrelação, efeitos de calendário);
3. **Prever** valores futuros — o objetivo final na maioria das aplicações de negócio e política pública.

Todo o restante deste documento prepara o caminho para esses três objetivos, na ordem acima.

---

## 2. Componentes e decomposição

### 2.1 Os três componentes

Na tradição inaugurada pela escola clássica (registrada também no livro-texto *Forecasting: Principles and Practice*, [capítulo 3](https://otexts.com/fpp3/decomposition.html)), costuma-se pensar toda série como composta de três partes:

- **Tendência-ciclo (T):** o movimento de longo prazo, incluindo flutuações lentas (“ciclos” econômicos, por exemplo);
- **Sazonalidade (S):** padrões sistemáticos ligados ao calendário (Natal, verão, meses do ano);
- **Resto/irregular (R):** tudo o que sobra — ruído, choques, combinações de fatores não modelados.

A forma mais simples de combiná-los é a **decomposição aditiva**:

$$y_t = S_t + T_t + R_t$$

Quando a intensidade da sazonalidade **cresce junto com o nível** da série (muito comum em economia), usamos a **decomposição multiplicativa**:

$$y_t = S_t \times T_t \times R_t$$

Em linguagem acessível: se as altas de dezembro dobram de tamanho quando o patamar das vendas dobra, faça uma decomposição **multiplicativa**; se as altas de dezembro têm sempre cerca de “+100 unidades”, não importando o nível, use a **aditiva**. Um truque frequente: aplicar **logaritmo** transforma uma relação multiplicativa em aditiva ($\log y_t = \log S_t + \log T_t + \log R_t$), permitindo usar ferramentas aditivas em seguida ([fpp3, cap. 3](https://otexts.com/fpp3/decomposition.html)). Há ainda modelos intermediários (como o “pseudo-aditivo” usado quando a série tem valores muito pequenos ou zeros que impediriam dividir por zero — caso da decomposição multiplicativa) ([ABS](https://www.abs.gov.au/websitedbs/d3310114.nsf/home/time+series+analysis:+the+basics)).

> **Ponto sutil (e importante):** um ciclo **sem duração fixa** não é a mesma coisa que sazonalidade. O número de peles de lince vendidas pela Hudson Bay Company, famoso conjunto de dados britânico, apresenta ciclos fortes, mas com timing imprevisível — e essa série é considerada **estacionária** (veja seção 3). Ciclos aperiódicos devem ser tratados com cautela: eles confundem muita gente ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)).

### 2.2 Métodos de decomposição

O objetivo da decomposição é estimar $\hat{S}$, $\hat{T}$ e $\hat{R}$ a partir de $y$. Os métodos mais usados:

1. **Média móvel simples (m-MA):** para estimar a tendência no período $t$, tira-se a média das observações numa janela centrada em $t$:

   $$\hat{T}_t = \frac{1}{m}\sum_{j=-k}^{k} y_{t+j}, \quad m = 2k+1$$

   Funciona como um “filtro passa-baixa”: elimina parte do ruído e da sazonalidade (se $m$ for igual ao período sazonal) e deixa a curva suave da tendência ([fpp3, cap. 3](https://otexts.com/fpp3/decomposition.html)).

2. **Médias móveis ponderadas:** pesos que somam 1 e são simétricos (por exemplo, [1/8, 1/4, 1/4, 1/4, 1/8]) produzem curvas mais suaves que a média móvel simples ([fpp3](https://otexts.com/fpp3/decomposition.html)).

3. **Decomposição clássica (1920s):** receita em passos — estimar tendência com média móvel (2×m quando o período $m$ é par; m ímpar caso contrário), subtrair, tirar a média das estações (todos os janeiros, todos os fev…) para achar o componente sazonal, ajustar a soma dos sazonal a zero e obter o resto ([fpp3](https://otexts.com/fpp3/decomposition.html)). Simples, mas assume sazonalidade constante ao longo dos anos.

4. **Métodos modernos:** **STL** (“Seasonal-Trend decomposition using LOESS”) e **X-13ARIMA-SEATS** (usado por agências estatísticas para ajuste sazonal oficial). São métodos mais flexíveis e robustos que a decomposição clássica — e é neles que as ferramentas atuais se apoiam ([fpp3](https://otexts.com/fpp3/decomposition.html)).

> **Cuidado epistemológico:** a decomposição não revela “componentes verdadeiros” ocultos na natureza; ela depende do método escolhido e das hipóteses feitas. Duas decomposições diferentes da mesma série podem contar histórias diferentes. Use decomposição para **explorar** e **comunicar**, não como prova.

### 2.3 Decomposição no R

No R há várias portas de entrada. As mais usadas:

```r
# Caminho 1: função base do R (stats) — decomposição STL
decomp <- stl(AirPassengers, s.window = "periodic")
plot(decomp)

# Caminho 2: pacote {forecast} — STL com previsão e ajuste automático
library(forecast)
fit <- stlm(AirPassengers, modelfunction = ar)   # STL + modelo AR sobre o destendenciado
plot(forecast(fit, h = 24))

# Caminho 3: abordagem moderna (tidy) com {fable}/{feasts}
# (carregue seus dados como tsibble e use stl() / STLM() — ver seção 10)
```

O uso de `stlm()`/`stlf()` permite combinar decomposição STL com um modelo (ARIMA, ETS…) sobre a série ajustada, produzindo previsões que respeitam o padrão sazonal ([documentação do pacote forecast](https://pkg.robjhyndman.com/forecast/)).

**Figura didática (gerada por simulação, não dados reais):**

![Figura 1 — Decomposição clássica de uma série simulada](figuras/fig1_componentes.png)

Observe como o painel ① mostra a série bruta; o painel ② recupera a curva lenta de tendência; o painel ③ captura o padrão sazonal repetido; e o painel ④ — o “resto” — é o que resta para o acaso explicar. Em séries reais, o resto também contém eventos raros (greves, crises, pandemia).

---

## 3. Estacionariedade e transformações

### 3.1 Definição

Uma série é **estacionária** quando suas propriedades estatísticas (média, variância, autocorrelação) **não dependem do momento** em que olhamos para ela. Em termos práticos: o jogo de dados não muda com o passar do tempo; qualquer janela da série “se comporta” de forma parecida ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)).

Consequências imediatas:

- Séries **com tendência** ou **com sazonalidade** não são estacionárias (o nível muda sistematicamente);
- **Ruído branco** (valores aleatórios independentes, com média e variância constantes) é estacionária;
- Série **cíclica sem período fixo pode ser estacionária** (caso do lince) — o critério não é “ter onda”, e sim poder prever o padrão de longo prazo.

Por que isso importa tanto? Porque **a maior parte dos modelos clássicos de séries temporais (ARIMA e parentes) assume estacionariedade** para a série que modela. Sem essa propriedade, os parâmetros estimados perdem interpretação estável e as previsões tendem a simplesmente replicar o último valor observado ([quantinsti](https://blog.quantinsti.com/stationarity/) e [fpp3](https://otexts.com/fpp3/arima.html)).

### 3.2 Testes formais (e suas limitações)

Dois testes dominam a prática:

| Teste | Hipótese nula (H₀) | Como interpretar |
|---|---|---|
| **ADF** (Augmented Dickey–Fuller) | Há **raiz unitária** (série não estacionária) | Estatística **mais negativa** rejeita H₀ → indício de estacionariedade ([Wikipedia – ADF](https://en.wikipedia.org/wiki/Augmented_Dickey%E2%80%93Fuller_test)) |
| **KPSS** (Kwiatkowski–Phillips–Schmidt–Shin) | Série é **estacionária** | p-valor pequeno → evidência contra estacionariedade |

Atenção: os testes **não concordam sempre** — e ambos têm **poder baixo** em amostras pequenas. A recomendação do livro-texto é usá-los como *triagem*, complementando-os com análise visual (gráfico da série e da série diferenciada) e com bom senso sobre o contexto econômico ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html); [análise macro – estacionariedade](https://analisemacro.com.br/estatistica-e-econometria/estacionariedade-de-series-temporais/)).

Intuição do teste ADF: estima-se uma equação do tipo

$$\Delta y_t = \alpha + \beta t + \gamma y_{t-1} + \delta_1 \Delta y_{t-1} + \cdots + \delta_{p-1}\Delta y_{t-p+1} + \varepsilon_t$$

e testa-se $\gamma = 0$ (ou seja, “o nível anterior não ajuda a prever a variação”; $\gamma < 0$ indica reversão à média → estacionária). “Augmented” significa incluir defasagens da diferença para limpar a autocorrelação remanescente ([Wikipedia – ADF](https://en.wikipedia.org/wiki/Augmented_Dickey%E2%80%93Fuller_test)).

### 3.3 Transformações para atingir estacionariedade

**(a) Diferenciação.** A ferramenta nº 1. Define-se a diferença de ordem 1 como $y'_t = y_t - y_{t-1}$ (a “variação de um período para o outro”), e a diferença sazonal como $y'_t = y_t - y_{t-m}$ (a “variação em relação ao mesmo mês do ano anterior”). Aplicar diferenças sucessivas até obter série estacionária é o “I” do ARIMA (o “d”):

- d = 0 → série já estacionária (na média);
- d = 1 → primeira diferença;
- d = 2 → segunda diferença (raro; implica que a previsão seguirá uma **linha reta**).

**Regra de ouro:** diferencie **o mínimo necessário**. Diferenciar demais **induz dinâmicas falsas** — cria autocorrelações artificiais onde não existem e destrói sinal interpretável ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)). Quando houver sazonalidade forte, prefira começar pela **diferença sazonal**; depois, se necessário, aplicar diferença adicional (a ordem entre elas não altera o resultado final, mas a escolha afeta a interpretabilidade durante o processo) ([fpp3](https://otexts.com/fpp3/arima.html)).

**(b) Transformações de estabilização de variância.** Se a variabilidade cresce com o nível (típico em séries econômicas), use escala logarítmica ou Box-Cox com parâmetro λ escolhido por critérios (ex.: método de Guerrero, implementado no pacote *fable* como `guerrero`) ([fpp3, caps. 3 e 9](https://otexts.com/fpp3/decomposition.html)). Depois da modelagem, reverta a transformação para reportar previsões na escala original.

> **Erro comum:** transformar sem necessidade distorce relações e dificulta a interpretação das previsões; não transformar quando há heterocedasticidade corrói a confiabilidade dos intervalos. Documente o que você fez e faça o caminho inverso nos resultados.

### 3.4 No R

```r
# Diferenças com o pacote {forecast} (diference() e diff())
library(forecast)
dd <- diff(AirPassengers)         # primeiras diferenças
dds <- diff(diff(AirPassengers))  # segundas diferenças (evite por hábito!)

# Testes de raiz unitária (pacote {tseries})
library(tseries)
adf.test(AirPassengers)           # ADF: p pequeno => rejeita raiz unitária
kpss.test(AirPassengers)          # KPSS: p pequeno => rejeita estacionariedade

# O auto.arima() faz esse teste internamente (default: KPSS),
# escolhendo automaticamente o número de diferenças d e D.
```

---

## 4. Autocorrelação: ACF e PACF

### 4.1 Intuição

**Autocorrelação** é a correlação de uma série com ela mesma, deslocada no tempo (defasada). Formalmente, $\rho_k$ mede a correlação entre $y_t$ e $y_{t-k}$. A função de autocorrelação (**ACF**) é o gráfico de $\rho_1, \rho_2, \dots$ contra as defasagens; a função de autocorrelação **parcial** (**PACF**) mede a correlação entre $y_t$ e $y_{t-k}$ *depois de remover o efeito* das defasagens intermediárias ([Nau, Duke](https://people.duke.edu/~rnau/411arim3.htm)).

Essas duas funções são a “impressão digital” de um processo e servem para **identificar as ordens p (AR) e q (MA)** dos modelos da próxima seção. O guia clássico (Nau) resume:

- **AR(p):** a ACF se **extingue gradualmente** (cauda exponencial ou senoidal amortecida), enquanto a PACF **corta** depois do lag p (picos significativos até p);
- **MA(q):** o oposto — a ACF **corta** depois do lag q, e a PACF se extingue aos poucos;
- **ARMA misto:** os dois comportamentos aparecem juntos, exigindo leitura conjunta.

**Figura didática:**

![Figura 2 — Assinaturas típicas de ACF e PACF para AR(1) e MA(1)](figuras/fig2_acf_pacf.png)

Detalhes de leitura:

- As linhas tracejadas nos gráficos são as bandas de referência (±1,96/√T): picos fora delas sugerem autocorrelação significativa;
- Antes de olhar ACF/PACF, trate tendência e sazonalidade (diferencie ou ajuste)! A ACF de uma série com tendência quase sempre “explode” e engana o diagnóstico ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html));
- ACF/PACF informam sobre **estrutura linear** — não capturam causalidade nem relação não linear. Correlação, aqui, é só padrão, nunca causa (seção 9).

### 4.2 No R

```r
# ACF e PACF brutas
acf(AirPassengers, lag.max = 24)
pacf(AirPassengers, lag.max = 24)

# Comparação: série bruta (com tendência/sazonalidade) vs tratada
par(mfrow = c(2, 2))
acf(AirPassengers, main = "ACF bruta")
acf(log(AirPassengers), main = "ACF em log")
acf(diff(log(AirPassengers)), main = "ACF log-diferenciada")
pacf(diff(log(AirPassengers)), main = "PACF log-diferenciada")

# Bandas de confiança (manual)
T <- length(AirPassengers)
banda <- 1.96 / sqrt(T)
cat("Banda:", round(banda, 3), "\n")
```

---

## 5. Modelos da família ARIMA e a metodologia Box–Jenkins

### 5.1 Peças do Lego: AR, MA e integração

**(a) AR(p) — autorregressivo.** O valor de hoje é explicado por uma combinação linear dos p valores anteriores:

$$y_t = c + \phi_1 y_{t-1} + \phi_2 y_{t-2} + \cdots + \phi_p y_{t-p} + \varepsilon_t$$

onde $\varepsilon_t$ é ruído branco. Para a série ser estacionária, os coeficientes precisam respeitar restrições (para AR(1): $-1 < \phi_1 < 1$; para AR(2), um conjunto de desigualdades envolvendo $\phi_1$ e $\phi_2$) ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)). Casos notáveis do AR(1):

| Parâmetros | Processo resultante |
|---|---|
| $\phi_1 = 0$, c = 0 | Ruído branco |
| $\phi_1 = 1$, c = 0 | **Passeio aleatório** (não estacionário) |
| $\phi_1 = 1$, c ≠ 0 | Passeio aleatório com *drift* (tendência) |
| $\phi_1 < 0$ | Oscilação ao redor da média |

**(b) MA(q) — médias móveis de erros.** O valor de hoje depende de choques passados:

$$y_t = c + \varepsilon_t + \theta_1 \varepsilon_{t-1} + \theta_2 \varepsilon_{t-2} + \cdots + \theta_q \varepsilon_{t-q}$$

Não confunda **modelo MA** (equação acima, para prever) com **média móvel de suavização** (operador de filtro, para estimar tendência) — são conceitos distintos com nomes parecidos ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)).

**(c) ARMA(p,q).** Combina as duas peças:

$$y_t = c + \underbrace{\phi_1 y_{t-1} + \cdots + \phi_p y_{t-p}}_{\text{parte AR}} + \varepsilon_t + \underbrace{\theta_1 \varepsilon_{t-1} + \cdots + \theta_q \varepsilon_{t-q}}_{\text{parte MA}}$$

**(d) ARIMA(p,d,q).** O “I” = *integrated* = **diferenciada**. Aplicar $(1-B)^d$ à série (onde $B$ é o operador de retrocesso, $B y_t = y_{t-1}$) equivale a diferenciar d vezes; então o ARIMA(p,d,q) é um ARMA(p,q) aplicado à série diferenciada d vezes, com polinômios em $B$ escritos de forma compacta:

$$(1 - \phi_1 B - \cdots - \phi_p B^p)(1-B)^d\, y_t = (1 + \theta_1 B + \cdots + \theta_q B^q)\,\varepsilon_t$$

O termo constante $c$ pode aparecer ou não; sua presença e o valor de $d$ mudam radicalmente o comportamento das previsões de longo prazo ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)):

| c | d | Comportamento da previsão de longo prazo | Intervalos |
|---|---|---|---|
| 0 | 0 | convergem para 0 | — |
| 0 | 1 | convergem para uma constante | crescem rápido (random walk) |
| 0 | 2 | seguem reta | crescem ainda mais rápido |
| ≠0 | 0 | convergem para a média histórica | largura ~constante |
| ≠0 | 1 | seguem reta (com inclinação) | crescem com $h$ |
| ≠0 | 2 | seguem parábola (não recomendado) | — |

**(e) SARIMA — sazonalidade.** Extende-se a notação com componentes sazonais: $\text{ARIMA}(p,d,q)(P,D,Q)_m$. Por exemplo, $(0,1,1)(0,1,1)_{12}$ significa: primeira diferença + MA de ordem 1, e no lado sazonal (periodicidade 12): primeira diferença sazonal + MA sazonal de ordem 1. Essa estrutura explica padrões repetitivos anuais com poucos parâmetros ([NCSS – Box-Jenkins](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/NCSS/The_Box-Jenkins_Method.pdf)).

> **Nota histórica.** Existem convenções distintas de sinais para os coeficientes MA (uns escrevem $-\theta_1\varepsilon_{t-1}$, outros $+\theta_1\varepsilon_{t-1}$). O conteúdo matemático é o mesmo; sempre confira a definição adotada pela fonte/software — mudanças de sinal trocam a interpretação de “choque” que se transmite positiva ou negativamente.

### 5.2 Metodologia Box–Jenkins

Box & Jenkins (1970) propuseram um **ciclo iterativo** para montar modelos ARIMA de forma disciplinada — descrito em detalhes nas fontes ([NCSS](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/NCSS/The_Box-Jenkins_Method.pdf), [análise macro](https://analisemacro.com.br/econometria-e-machine-learning/modelos-arima-e-metodologia-box-jenkins/)):

1. **Identificação** — examinar gráfico da série, sazonalidade, tendência; aplicar testes de raiz unitária; computar ACF/PACF para propor ordens $(p,d,q)$;
2. **Estimação** — estimar os coeficientes por mínimos quadrados não lineares ou máxima verossimilhança;
3. **Verificação (diagnóstico)** — testar resíduos: devem parecer **ruído branco** (sem autocorrelação; testes como Ljung-Box; gráficos). Se ainda houver estrutura, volta-se ao passo 1;
4. **Previsão** — usar o modelo validado para projetar. Repetir o ciclo periodicamente, pois séries reais mudam de comportamento.

O ciclo é *iterativo*: diagnóstico ruim ⇒ nova identificação. Hoje, quase ninguém faz manualmente todas as etapas: funções como `auto.arima()` fazem busca automática sobre as ordens, usando critério de informação (AICc) como bússola ([Hyndman & Khandakar, JSS 2008](https://www.jstatsoft.org/index.php/jss/article/view/v027i03)). Mas conhecer o ciclo ajuda a **interpretar diagnósticos e falhas** — e é obrigatório entender o que a automação está fazendo.

### 5.3 Seleção de ordens e critério de informação

Para comparar candidatos, usamos **AIC**, **AICc** e **BIC**, todos com a mesma lógica: penalizam verossimilhança alta (melhor ajuste) pelo número de parâmetros (evitar overfitting). Regra prática: **menor valor ganha**; AICc é preferível em amostras pequenas ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html)). O `auto.arima()` considera dezenas de configurações (padrão: até 5º grau, sazonal até 2×2, busca *stepwise* rápida, IC = AICc) e devolve o melhor encontrado ([documentação da função](https://www.rdocumentation.org/packages/forecast/versions/9.0.1/topics/auto.arima)).

**Saída real de um `auto.arima`** (fonte: material didático de diagnose, [r-statistics.co](https://r-statistics.co/ARIMA-Diagnostics-in-R.html)):

```r
library(forecast)
fit_good <- auto.arima(WWWusage)
fit_good
#> Series: WWWusage
#> ARIMA(1,1,1)
#> Coefficients:
#>        ar1      ma1
#>      0.6504  0.5256
#> s.e.  0.0842  0.0896
#> sigma^2 = 9.995: log likelihood = -254.15
#> AIC=514.3 AICc=514.55 BIC=522.08
```

Interpretação: a série foi diferenciada uma vez (d = 1); o modelo tem um termo AR(1) com coeficiente 0,65 e um MA(1) com 0,53 — e os p-valores (não mostrados) e o diagnóstico de resíduos precisam confirmar que isso é suficiente.

**Exemplo de saída de um `Arima`** (reproduzido de [fpp3, cap. 9](https://otexts.com/fpp3/arima.html)):

```r
#> Model: ARIMA(2,0,1) w/ mean
#> Coefficients:
#>          ar1      ar2      ma1  constant
#>       1.6764  -0.8034  -0.6896    2.5623
#> sigma^2 estimated as 8.046:  log likelihood=-141.6
#> AIC=293.1   AICc=294.3   BIC=303.4
```

Isso corresponde à equação (arredondada):

$$y_t = 2{,}56 + 1{,}68\,y_{t-1} - 0{,}80\,y_{t-2} - 0{,}69\,\varepsilon_{t-1} + \varepsilon_t$$

### 5.4 No R — receitas básicas

```r
library(forecast)

# (1) Modelo ARIMA com ordens fixas
m1 <- Arima(AirPassengers, order = c(1, 1, 1))
report(m1)                     # imprime a equação estimada

# (2) Busca automática (inclui sazonalidade por default)
m2 <- auto.arima(AirPassengers)
report(m2)

# (3) Previsões + intervalos
f2 <- forecast(m2, h = 24, level = c(80, 95))
plot(f2)

# (4) Diagnóstico dos resíduos (deve parecer ruído branco)
checkresiduals(m2)
```

`checkresiduals()` reúne: gráfico dos resíduos no tempo, ACF dos resíduos e o **teste de Ljung-Box**, cujo p-valor “idealmente deve ser maior que 0,05” para não rejeitarmos a hipótese de resíduos sem autocorrelação ([r-statistics.co](https://r-statistics.co/ARIMA-Diagnostics-in-R.html)).

---

## 6. Previsão e validação (backtesting)

### 6.1 Por que a validação cruzada comum não serve

Em séries temporais, **não se pode embaralhar os dados**: o futuro não pode ajudar a prever o passado. Portanto, k-fold CV clássico (partição aleatória) geralmente **vaza informação futura** e superestima o desempenho real. A alternativa padrão é a **validação por origem rolante** (*rolling origin* / *walk-forward*):

- divide-se a série em blocos crescentes de treino;
- para cada bloco, ajusta-se o modelo e prevê-se **à frente** (um ou vários passos);
- compara-se a previsão com o valor que, de fato, ocorreu;
- repete-se para todas as origens; ao final, agregam-se os erros (média, RMSE etc.).

Diagrama mental (veja a Figura 3a): os segmentos azuis são treinos que crescem; cada segmento vermelho é o ponto de teste previsto a partir daquele treino; o processo “rola” para o futuro ([Hyndman, tscv](https://robjhyndman.com/hyndsight/tscv/)).

Para **múltiplos passos à frente**, a ideia é a mesma, mas mede-se o erro de cada horizonte separadamente (ex.: 1, 2, 3, 4 trimestres), rodando o procedimento para cada horizonte desejado ([Hyndman, tscv](https://robjhyndman.com/hyndsight/tscv/)).

### 6.2 No R — duas rotas

**Rota 1 — pacote {forecast}:**

```r
library(forecast)
# Erros de previsão 1-passo com origem rolante (modelo rwf = random walk com drift)
e <- tsCV(dj, rwf, drift = TRUE, h = 1)
sqrt(mean(e^2, na.rm = TRUE))   # RMSE por validação rolling-origin
```

Comparar esse RMSE com o RMSE *in-sample* dos resíduos do mesmo modelo é um excelente teste de otimismo excessivo (a diferença mostra quanta “sorte” o ajuste pegou) ([Hyndman, tscv](https://robjhyndman.com/hyndsight/tscv/)).

**Rota 2 — framework {fable} (tidy):**

```r
library(fpp3)

# Cria os folds (treino crescente a partir de 12 pontos)
beer <- aus_production %>%
  select(Beer) %>%
  stretch_tsibble(.init = 12, .step = 1)

# Ajusta ETS em cada fold e prevê 1 ano à frente
fc <- beer %>%
  model(ETS(Beer)) %>%
  forecast(h = "1 year") %>%
  group_by(.id) %>% mutate(h = row_number()) %>%
  ungroup() %>%
  as_fable(response = "Beer", distribution = Beer)

# Compara os erros por horizonte
fc %>% accuracy(aus_production, by = c("h", ".model")) %>%
  select(h, RMSE)
#> # A tibble: 4 × 2
#>       h  RMSE
#>   <int> <dbl>
#> 1     1  17.1
#> 2     2  16.7
#> 3     3  18.1
#> 4     4  19.2
```

(Saída reproduzida de [Hyndman, tscv-fable](https://robjhyndman.com/hyndsight/tscv-fable/).) Repare: os erros do **horizonte curto** também contam — e muitas vezes o melhor modelo para h=1 não é o melhor para h=4. Escolha o horizonte pela **decisão** que você precisa apoiar (seção 8).

### 6.3 Figura ilustrativa

![Figura 3 — Origem rolante e comportamento do erro conforme o horizonte](figuras/fig3_validacao.png)

A intuição central: (a) toda avaliação honesta de previsão ocorre **no futuro** de cada ponto de treino; (b) o erro tipicamente **cresce com o horizonte**, mas a taxa de crescimento diferencia os modelos — a curva do passeio aleatório sobe com $\sqrt{h}$, enquanto um modelo com sinal útil pode subir mais devagar.

---

## 7. Métricas de acurácia: como medir se o modelo é bom

### 7.1 Catálogo essencial

Considere $y_t$ o valor observado, $\hat{y}_t$ a previsão, $T$ o número de pontos avaliados. Em notação compacta (e seguindo a discussão de [Hyndman & Koehler, 2006](https://robjhyndman.com/papers/foresight.pdf) e do [capítulo 5 do fpp3](https://otexts.com/fpp3/toolbox.html)):

| Métrica | Fórmula | Escala | Melhor uso | Cuidados |
|---|---|---|---|---|
| **MAE** (erro absoluto médio) | $\frac{1}{T}\sum \|y_t - \hat{y}_t\|$ | unidade da série | Interpretação direta numa série única | Não comparável entre séries com escalas diferentes |
| **MSE / RMSE** | $\frac{1}{T}\sum (y_t-\hat{y}_t)^2$ ; $\sqrt{\text{MSE}}$ | unidade² / unidade | Penaliza erros grandes | Sensível a outliers; RMSE volta à escala original |
| **MAPE** | $\frac{100\%}{T}\sum \left|\frac{y_t-\hat{y}_t}{y_t}\right|$ | % | Séries de escala parecida | **Falha com zeros ou valores próximos de zero** (divisão por ~0); ignora direção |
| **MASE** (erro absoluto escalonado) | $\text{mean}(|q_j|)$ com $q_j = \frac{e_j}{\frac{1}{T-m}\sum_{t=m+1}^{T}\|y_t-y_{t-m}\|}$ | livre de escala | **Comparar séries diferentes**; robusto a zeros | Interpretação “vs. ingênuo”: < 1 melhor que naïve; > 1 pior |
| **sMAPE** | versão simétrica do MAPE | % | Evitar — apenas reconhecer | Pode dar valores infundados/negativos; especialistas recomendam não usar ([fpp3](https://otexts.com/fpp3/toolbox.html)) |
| **Viés (MPE / ME)** | $\frac{1}{T}\sum (y_t-\hat{y}_t)$ | unidade | Detectar erro sistemático (sub/superprever) | Pode cancelar positivo/negativo |

Sobre o **MASE** vale detalhar: é a métrica proposta por Hyndman & Koehler (2006) para superar os problemas do MAPE — escalona o erro pela *média absoluta das variações de um passo do próprio histórico* (o erro do método ingênuo). Assim, $|q_j| < 1$ significa “acertamos mais que a previsão naïve”; a média dessas razões dá o MASE. Ele funciona com zeros, com escalas arbitrárias e através de séries — por isso virou o padrão recomendado para comparações amplas ([Hyndman & Koehler](https://robjhyndman.com/papers/foresight.pdf); [fpp3 toolbox](https://otexts.com/fpp3/toolbox.html)).

**Escolha por situação:**

- Uma única série, comunicação fácil → **MAE**;
- Comparar modelos em várias séries de magnitudes diferentes → **MASE**;
- Medir magnitude de erros extremos → **RMSE** (mas combine com MAE);
- Sempre complementar com **viés** (será que o modelo erra sistematicamente para cima?) e, idealmente, **skill score** relativo ao método ingênuo.

### 7.2 No R

```r
# O pacote {forecast} traz accuracy() pronta; o {fabletools} idem,
# com colunas .model / .type (Treino/Teste) e medidas:
#   ME, RMSE, MAE, MPE, MAPE, MASE, RMSSE, ACF1
accuracy(modelo)                       # métricas in-sample (resíduos)
accuracy(previsao, dados_completos)     # métricas out-of-sample

# Exemplo típico de saída (fabletools):
#> .model  .type       ME   RMSE    MAE   MPE   MAPE   MASE  RMSSE   ACF1
#> ets     Training  0.321  15.8   12.0  0.103  2.87   0.752 0.795  -0.177
#> ets     Test      ...    ...    ...   ...    ...    ...   ...    ...
```

Note a linha **“Test”**: para avaliar de verdade, forneça ao `accuracy()` o dataset com os **valores futuros reais** — ele precisará deles para calcular os erros das previsões. Sem dados futuros, não há avaliação honesta ([fabletools – accuracy](https://fabletools.tidyverts.org/reference/accuracy.html)).

---

## 8. Incerteza: intervalos de previsão e comunicação

### 8.1 Conceito

Uma previsão pontual (“IPCA será 0,4% no próximo mês”) é quase certamente errada. Uma previsão **intervalar** entrega faixas com níveis de confiança — tipicamente 80% e 95% — e comunica honestamente a incerteza. Três fatos a saber ([fpp3, caps. 5 e 9](https://otexts.com/fpp3/toolbox.html)):

1. Os intervalos **crescem com o horizonte** (mais tempo = mais incerteza);
2. O quanto crescem depende de **d**: séries diferenciadas (d ≥ 1) têm intervalos que abrem mais rápido — para d = 0 eles tendem a uma largura constante;
3. Intervalos de 95% ≈ $\hat{y}_t \pm 1{,}96 \times \text{desvio-padrão do erro}$ em condições regulares; modelos modernos fornecem **distribuições preditivas** completas (pacote *fable* gera distribuições, não apenas intervalos).

### 8.2 Comunicação

- Prefira **informar o intervalo** ao gestor/leitor, não só o ponto;
- Diga explicitamente o nível (80%? 95%?) e mostre-o visualmente (fan chart);
- Lembre que intervalos **cobrem a incerteza estatística**, não a incerteza estrutural (mudanças de regime podem invalidar o modelo inteiro — veja seção 9);
- Ao comparar métodos, use métricas adequadas + backtesting; não troque o modelo por um alívio de 1% no MAPE em prejuízo de robustez.

---

## 9. Armadilhas, limitações e controvérsias

Esta é a parte que separa quem “usa modelos” de quem “confia demais em modelos”. Os problemas abaixo são clássicos e bem documentados.

### 9.1 Correlação (espúria) não é causação

Duas séries podem se mover juntas **sem qualquer relação causal**. O exemplo de manual: vendas de sorvete e mortes por afogamento crescem juntas no verão — o fator de confusão é a temperatura ([Wikipedia – Spurious relationship](https://en.wikipedia.org/wiki/Spurious_relationship)). Existe até um site dedicado a catálogos de correlações absurdas reais, construídas com dados públicos — da safra de queijo consumida per capita ao número de doutoramentos concedidos: tudo “correlaciona” se você procurar o bastante ([Tyler Vigen – Spurious Correlations](https://www.tylervigen.com/spurious-correlations)).

Versão mais traiçoeira para séries temporais: a **regressão espúria** (Granger & Newbold, 1974) — regredir duas séries **não estacionárias** e independentes tende a produzir R² alto, estatísticas significativas e conclusão errada ([Wikipedia](https://en.wikipedia.org/wiki/Spurious_relationship); [metricgate](https://metricgate.com/blogs/spurious-regression-nonstationary-series/)). Qualquer tendência comum basta para enganar. Remédio: tratar não estacionariedade (diferenciar, cointegração quando apropriado) e **pensar** antes de conectar causas.

> **Regra de ouro:** nenhuma correlação, por mais forte e “significativa” que seja, substitui um mecanismo explicativo. Antes de afirmar “X move Y”, pergunte: (1) X e Y compartilham um terceiro fator? (2) ambos apenas tiveram tendências coincidentes? (3) a relação existe dentro de um subsistema controlado?

### 9.2 Overfitting e overdiferenciação

Modelos com muitos parâmetros memorizam o ruído do treino e fracassam no futuro. Sintomas: ajuste quase perfeito *in-sample*, desempenho pífio *out-of-sample*. Práticas preventivas:

- Prefira modelos parcimoniosos; valide **fora da amostra** sempre (seção 6);
- **Não diferencie mais que o necessário** — a overdose cria autocorrelação artificial (“dinâmica falsa”) que você vai tentar modelar inutilmente ([fpp3, cap. 9](https://otexts.com/fpp3/arima.html));
- Desconfie de AICc empurrando ordens altas em séries curtas; restrinja `max.p/max.q` em aplicações críticas;
- Use validação por origem rolante como critério final, não o critério de informação.

### 9.3 Diagnóstico de resíduos: nuances do Ljung-Box

O teste de Ljung-Box verifica se os resíduos ainda contêm autocorrelação em um grupo de defasagens:

$$Q = T(T+2)\sum_{k=1}^{\ell}(T-k)^{-1}r_k^2$$

onde $r_k$ são as autocorrelações dos resíduos. A questão delicada são os **graus de liberdade**: tradicionalmente usa-se $\ell - p - q$ quando um ARMA(p,q) é ajustado; debate recente mostrou que para outros modelos (ETS, regressões, ajustes sazonais) não há resposta universal, e a implementação do `{forecast}` mudou seu cálculo em 2023 ([Hyndman, ljung_box_df](https://robjhyndman.com/hyndsight/ljung_box_df.html)). Leitura prática:

- p > 0,05 → não rejeitamos ausência de autocorrelação; p < 0,05 → modelo deixou estrutura para trás;
- p muito alto (ex.: 0,95) também merece desconfiança: com amostras pequenas, o teste tem **poder baixo** e aceita qualquer coisa;
- Sempre combine com inspeção visual dos resíduos e da ACF dos resíduos ([r-statistics.co](https://r-statistics.co/ARIMA-Diagnostics-in-R.html)).

```r
library(forecast)
checkresiduals(modelo)   # gráfico + ACF dos resíduos + Ljung-Box
lb<- Box.test(residuals(modelo), lag = 12, type = "Ljung-Box")
```

### 9.4 Quebras estruturais e mudanças de regime

Modelos estimados sobre décadas de dados assumem que o **mecanismo gerador** permaneceu estável. Crises mostram o contrário. No choque da COVID-19, por exemplo, economistas do Banco Central Europeu documentaram que: (i) parâmetros de curvas de Phillips e VARs mudaram notavelmente com a inclusão das observações da pandemia; (ii) permitir resíduos com **caudas gordas** (distribuição t de Student) melhora a robustez; (iii) condicionamento em informações externas (projeções oficiais, expectativas de analistas) ajuda a navegar cenários atípicos ([Bobeica & Hartwig, ECB WP 2558](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2558~22b223a7c6.en.pdf)). Em outras palavras: em tempos incomuns, **dados históricos dizem menos**; combine o modelo com informação externa e monitore quebras (CUSUM, testes de ruptura, controle de indicadores de regime).

**Tradução prática:** todo modelo é condicionado ao “regime” vigente no histórico. Se a economia muda (lei nova, pandemia, choque geopolítico), reavalie — não confie no modelo só porque “funcionou” em 2010–2019.

### 9.5 Checklist anti-ilusão

- [ ] Gráfico da série e das séries transformadas antes de qualquer modelo
- [ ] ACF/PACF lidos com cuidado (tendência/sazonalidade tratadas)
- [ ] Diferenciação mínima; log/Box-Cox só quando estabiliza variância
- [ ] Validação fora da amostra por origem rolante (nunca k-fold aleatório)
- [ ] Métricas complementares (MAE + MASE + viés), não uma só
- [ ] Resíduos testados (Ljung-Box + gráfico) e entendidos
- [ ] Intervalos de previsão relatados, com nível explicitado
- [ ] Possíveis quebras/regimes considerados; cenário externo incorporado
- [ ] Correlação jamais tratada como causalidade sem teoria
- [ ] Simplicidade preferida até que os dados peçam complexidade

---

## 10. Exemplos completos em R

A seguir, cinco fluxos de trabalho completos. Cada um foi baseado na documentação oficial dos pacotes; rode-os no seu ambiente R (versões ≥ 4.2 e pacotes recentes recomendados).

### 10.0 Instalação e dados de exemplo

```r
# Pacotes essenciais
install.packages(c("forecast", "fable", "tsibble", "fpp3", "prophet", "tseries"))
# Conexão com dados do Banco Central (duas opções):
install.packages(c("rbcb", "GetBCBData"))
```

O livro *Forecasting: Principles and Practice* acompanha o pacote *fpp3* (dados e código); qualquer série própria pode entrar no lugar dos datasets de exemplo.

### 10.1 Exemplo A — Do zero: simular, ajustar e prever AR(1)

```r
set.seed(1)
n <- 300
e <- rnorm(n)
y <- numeric(n)
for (t in 2:n) y[t] <- 0.7 * y[t-1] + e[t]   # AR(1) com phi = 0.7

plot(y, type = "l", main = "Série AR(1) simulada (phi = 0.7)")

library(forecast)
fit <- auto.arima(y)            # deve retornar algo próximo de ARIMA(1,0,0)
report(fit)                     # leia a equação estimada: o "phi" estará em ar1
par(mfrow = c(2, 1))
plot(fit)                       # ajuste
plot(forecast(fit, h = 50))     # previsões com intervalos
```

**O que esperar:** `ar1` próximo de 0,7, constante pequena; AICc perto do mínimo teórico; previsões de longo prazo convergindo para a média da série (efeito típico de d = 0 com constante ≠ 0).

### 10.2 Exemplo B — Dados reais do BCB: baixar, transformar e modelar o IPCA

```r
# --- Rota 1: pacote {rbcb} ---
ipca_df <- rbcb::get_series(c(ipca = 433), last = 120)   # tibble com colunas date e ipca

# --- Rota 2: pacote {GetBCBData} ---
# ipca_df <- GetBCBData::gbcbd_get_series(id = 433,
#                first_date = "2012-01-01", last_date = Sys.Date(),
#                format_data = "long")

# Montar um tsibble moderno (formato tidy) — padrão do {fpp3}
library(dplyr)
library(tsibble)
ipca_tbl <- ipca_df |>
  mutate(mes = yearmonth(date), valor = ipca) |>
  select(mes, valor) |>
  as_tsibble(index = mes)

library(fable)
modelo <- ipca_tbl |>
  model(ARIMA(log(valor)))       # log + ARIMA automático
report(modelo)                   # equação e AICc

prev <- forecast(modelo, h = "1 year")
prev |> autoplot(ipca_tbl) + labs(y = "IPCA (log)", title = "Previsão do IPCA")
```

(Observação: `autoplot()` funciona direto sobre objetos de previsão do *fable*; o fluxo “model() → forecast() → autoplot() é o canônico do framework.”) Ajuste o recorte temporal ao seu interesse; para dados mensais, o pacote cuida da estrutura temporal desde que o índice esteja como `yearmonth`.

### 10.3 Exemplo C — Decomposição STL + previsão (via {forecast})

```r
library(forecast)
fit_stl <- stlm(AirPassengers, modelfunction = ar)  # STL + modelo AR
plot(forecast(fit_stl, h = 24), main = "STL + AR: previsão")

# Alternativa pura do STL (sem previsão, só análise):
decomp <- stl(AirPassengers, s.window = "periodic")
plot(decomp)
```

**O que olhar:** a coluna “seasonal” do objeto STL; se o resto ainda mostrar padrão, tente outra janela (`s.window`) ou modelo mais flexível. O comando `stlf()` aplica a mesma ideia com previsão automática.

### 10.4 Exemplo D — Backtest sério com origem rolante ({fable})

Template genérico, pronto para adaptar (segue exatamente o fluxo documentado por [Hyndman](https://robjhyndman.com/hyndsight/tscv-fable/)):

```r
library(fpp3)

dados <- sua_base |>                      # tsibble com coluna `valor`
  stretch_tsibble(.init = 24, .step = 1)   # 24 pontos mínimos por treino

fc <- dados |>
  model(ARIMA(valor)) |>                  # ou ETS(valor), ou TBATS...
  forecast(h = "12 months") |>
  group_by(.id) |> mutate(h = row_number()) |>
  ungroup() |>
  as_fable(response = "valor", distribution = valor)

fc |>
  accuracy(sua_base_completa, by = c("h", ".model")) |>
  select(h, MAE, RMSE, MASE)
```

**Diagnóstico esperado:** tabela de erros por horizonte. Compare modelos (ETS vs ARIMA vs TBATS) na mesma tabela; escolha pela métrica que conversa com sua decisão (h! h! h!).

### 10.5 Exemplo E — {prophet}: modelo de componentes para dados diários

O *prophet* decompõe a série em **quatro partes aditivas**: tendência flexível (linear por partes, com pontos de mudança automáticos), sazonalidades periódicas (Fourier), efeitos de feriados/eventos e ruído ([material introdutório r-statistics.co](https://r-statistics.co/Prophet-in-R.html)):

```r
library(prophet)

# Formato obrigatório: coluna ds (datas) e y (numérica)
dados <- data.frame(
  ds = as.Date("2022-01-01") + 0:1094,
  y = rnorm(1095, mean = 200, sd = 15) + seq_len(1095) * 0.1
)

modelo <- prophet(dados)                    # tendência + sazonalidades automáticas
futuro <- extend_dataframe(model@df, period = 365, freq = "D")
prev   <- predict(modelo, future = futuro)

plot(modelo, prev)                           # previsão global
plot_components(modelo)                      # trend, sazonalidades, holidays

# Adicionar feriados brasileiros:
# carnaval <- data.frame(holiday = "Carnaval",
#                        ds = as.Date(c("2023-02-20","2024-02-12","2025-03-03")),
#                        lower = -7, upper = 7)
# modelo <- add_holidays(modelo, carnaval)
```

**Quando brilha:** séries diárias com sazonalidades fortes (dias da semana + anual), feriados relevantes e lacunas toleráveis. **Quando evita:** fenômenos de alta frequência intradiária, séries curtas (precisa de ~≥ 1 ano para anual), ou quando você deseja rigor estatístico bayesiano (prophet é “engenharia de previsão”, com inferência aproximada).

---

## 11. Estudo de caso: inflação e juros no Brasil

Dados reais consultados nas **séries SGS do Banco Central do Brasil** (base oficial; códigos SGS indicados) para ilustrar tudo o que discutimos — com o cuidado de distinguir **fato** (valor publicado) de **expectativa** (Focus) e de **inferência** (interpretação).

### 11.1 Dois fatos e duas expectativas (consulta de 05/11/2026)

**Fatos (séries oficiais):**

- **IPCA** (SGS 433): **+0,82%** em setembro/2026 (variação mensal); acumulado em 12 meses até setembro/2026: **4,58%** (SGS 13522);
- **Meta Selic** (SGS 432): **13,75%** ao ano, definida em setembro/2026 — depois de um ciclo de corte iniciado em março/2026; pico recente de **15,00%** em junho/2025.

**Expectativas (Focus, pesquisa semanal com analistas; consulta de 23/10/2026):**

- **IPCA 2026:** mediana **≈ 5,19%** (142 respondentes); Top-5: 5,20%;
- **IPCA 2027:** mediana ≈ 4,34% (todos) / 4,47% (Top-5);
- **Selic fim de 2026:** mediana ≈ 13,25% ao ano.

**Dólar comercial** (SGS 1): R$ 5,4778 em 04/11/2026.

> **Interpretação cautelosa:** Focus é **expectativa de mercado**, não previsão oficial do BCB; as séries acima retratam um instante (novembro/2026) e podem sofrer revisões metodológicas. Nada aqui constitui recomendação de investimento.

### 11.2 Panorama visual

![Figura 4 — IPCA mensal e acumulado 12 meses (painel superior); Meta Selic (painel inferior)](figuras/fig4_brasil.png)

O painel (a) mostra a natureza **volátil e ruidosa** de uma série mensal de inflação: meses negativos (agosto/2025: −0,11%; agosto/2026: −0,32%) convivem com saltos (fevereiro/2025: +1,31%; setembro/2026: +0,82%). A linha laranja — acumulado em 12 meses — é o número que a imprensa cita; note que ela suaviza os soluços mensais, exatamente como um filtro de média móvel faz com a série bruta. Já o painel (b) exibe a Meta Selic como uma “escada”: decisões discretas do Copom em reuniões ordinárias.

### 11.3 Como trazer esses dados para o R

**Via API oficial (padrão documentado):** a interface JSON do SGS segue o padrão

```
https://api.bcb.gov.br/dados/serie/bcdata.sgs.{codigo}/dados?formato=json&dataInicial=dd/mm/aaaa&dataFinal=dd/mm/aaaa
```

Desde 26/03/2025, consultas exigem filtro de datas limitado a **10 anos**; há também o endpoint de últimos N valores, **limitado a N ≤ 20** ([Portal de Dados Abertos do BCB – recurso SGS 11](https://dadosabertos.bcb.gov.br/dataset/11-taxa-de-juros---selic/resource/b73edc07-bbac-430c-a2cb-b1639e605fa8)). Para consultar sem escrever requisições HTTP, use os pacotes R:

```r
# Opção 1 — {rbcb}: séries como tibble/xts/ts
library(rbcb)
ipca <- rbcb::get_series(c(IPCA = 433), last = 60, as = "xts")
selic <- rbcb::get_series(c(Selic = 11), last = 60, as = "xts")
expectativas <- rbcb::get_market_expectations(type = "annual",
                                              variable = "IPCA", year = 2026)

# Opção 2 — {GetBCBData}: data.frame plano, ideal para juntar séries
ipca_df <- GetBCBData::gbcbd_get_series(
  id = 433, first_date = "2010-01-01", last_date = Sys.Date(),
  format_data = "long"
)
```

**Sugestão de exercício (completo):** baixe IPCA (433) e meta Selic (432), monte um `tsibble` com `yearmonth()`, estime um `ARIMA(log(ipca))` via `{fable}`, faça um backtest de 3 anos com `stretch_tsibble()` e compare com um `ETS`. Depois, comente o desempenho por horizonte (3, 6, 9, 12 meses) usando MAE e MASE. Você verá na prática os temas deste guia: transformação, sazonalidade, identificação, validação e métricas.

---

## 12. Glossário

| Termo | Significado |
|---|---|
| **ACF / PACF** | Funções de (parcial) autocorrelação; medem dependência de uma série com defasagens de si mesma |
| **AR(p)** | Modelo autorregressivo: valor atual explicado por p valores anteriores |
| **MA(q)** | Modelo de médias móveis: valor atual explicado por q choques anteriores |
| **ARIMA(p,d,q)** | ARMA aplicado a série diferenciada d vezes; “I” = integrated |
| **SARIMA** | Versão sazonal do ARIMA, com componentes sazonais (P,D,Q) de período m |
| **ETS** | Modelos de suavização exponencial com erro/tendência/sazonalidade (Aditivo/Multiplicativo) |
| **Estacionária** | Propriedades estatísticas constantes no tempo |
| **Raiz unitária** | Propriedade de séries não estacionárias onde choques não se dissipam (passeio aleatório) |
| **Passeio aleatório (random walk)** | Caso especial de AR(1) com φ₁ = 1 e sem constante; previsões = último valor |
| **Drift** | Tendência constante adicionada ao passeio aleatório |
| **Diferença (d, D)** | Operação yₜ − yₜ₋₁ (ou yₜ − yₜ₋ₘ para diferença sazonal); remove tendência/sazonalidade |
| **Transformação Box-Cox** | Família de transformações (inclui log) para estabilizar variância |
| **Decomposição** | Separação em tendência/sazonalidade/resto (clássica, STL, X-13) |
| **STL** | “Seasonal-Trend decomposition using LOESS” — método moderno de decomposição |
| **Backtesting / origem rolante** | Validação em que treino cresce e teste sempre fica no futuro |
| **Horizonte (h)** | Quantos passos à frente se está prevendo |
| **MAE / RMSE / MAPE / MASE** | Métricas de acurácia (ver seção 7) |
| **Skill score** | Ganho percentual vs. método ingênuo (naïve) |
| **Intervalo de previsão** | Faixa com nível de confiança (ex.: 80%, 95%) para o valor futuro |
| **Viés (bias)** | Tendência sistemática de subestimar/superestimar |
| **Quebra estrutural** | Mudança no mecanismo gerador da série (regime novo) |
| **Correlação espúria** | Associação aparente sem causalidade (confundidor ou coincidência) |
| **Pseudo R² / R² de regressão espúria** | R² enganosamente alto ao regredir séries não estacionárias independentes |
| **Naïve / Seasonal naïve** | Método ingênuo de referência: repetir o último valor (ou o mesmo mês do ano anterior) |
| **Whitening / resíduos brancos** | Resíduos sem autocorrelação — o alvo do diagnóstico |

---

## 13. Como continuar aprendendo

**Livros e cursos gratuitos (fontes verificadas):**

- **[Forecasting: Principles and Practice, 3ª ed.](https://otexts.com/fpp3/)** — Hyndman & Athanasopoulos. Gratuito e online, com todo código em R; é a espinha dorsal deste guia. Há edições anteriores ([fpp2](https://otexts.com/fpp2/)) e uma tradução/versão em português em desenvolvimento ([fpptr](https://otexts.com/fpptr/));
- **[Análise de Séries Temporais](https://users.ime.usp.br/~pam/ST.html)** — Morettin & Toloi (ABE/Blucher, 3ª ed., 2 vols.). Referência brasileira clássica, com conjuntos de dados e respostas disponíveis no site dos autores;
- **[Curso gratuito “Séries Temporais em R”](https://analisemacro.com.br/wp-content/uploads/2015/02/CursoCompleto.pdf)** (Ricardo Lima; hospedado na Análise Macro) — introdução em português ao estilo econométrico brasileiro. Nota: o próprio autor avisa que o material está em construção;
- **[CRAN Task View: Time Series Analysis](https://cran.r-project.org/view=TimeSeries)** — lista mantida por Hyndman & Killick com os pacotes organizados por tarefa (classes de dados, decomposição, testes, modelagem, previsão). Útil como mapa do ecossistema.

**Trilha sugerida (4 semanas, ~1h/dia):**

1. **Semana 1** — Conceitos e visualização: capítulos 1–3 do fpp3; refazer figuras de decomposição com seus dados;
2. **Semana 2** — Estacionariedade, ACF/PACF: capítulos 4 e 9; exercícios de identificação com séries simuladas;
3. **Semana 3** — ARIMA na prática: `auto.arima`, diagnóstico, previsões; começar backtesting com `stretch_tsibble`;
4. **Semana 4** — Projeto miniatura: pegar uma série real do SGS (IPCA, câmbio, Selic), tratar, modelar, backtestar, escrever conclusões com intervalos e métricas.

**Comunidades e suporte:** perguntas de forecasting no Cross Validated/Stats Stack Exchange; questões de R no Stack Overflow com tags `r` + `forecast`; blogs dos autores do {fable} ([Mitchell’s blog](https://mitchelloharawild.com/blog/fable/)) e [Hyndsight](https://robjhyndman.com/hyndsight/).

---

## 14. Limitações explícitas deste documento

1. **Sem execução R neste texto.** Os trechos de R foram montados a partir da documentação oficial dos respectivos pacotes e fontes citadas; não houve execução de código R na sessão em que este documento foi preparado. Antes de usar em produção, valide cada comando nas versões dos pacotes instaladas no seu ambiente. Funções e argumentos evoluem (especialmente {fable}, {forecast} e {prophet});
2. **Simulações didáticas.** As Figuras 1–3 foram geradas por simulação em Python (numpy/matplotlib) exclusivamente para fins educacionais; não representam dados reais. A Figura 4 usa dados oficiais do BCB (SGS), compilados no ambiente de consulta;
3. **Instante temporal.** Os números do estudo de caso (seção 11) refletem consultas de novembro/2026 às bases do BCB (séries e Focus). Séries são revisadas e expectativas mudam toda semana — não trate como “estado atual” em outro momento;
4. **Escopo.** Este guia cobre modelagem univariada e previsão clássica. Não aborda: modelos multivariados (VAR/VECM), aprendizado estatístico moderno (XGBoost, redes neurais), séries de alta frequência, espaço de estados avançado, Bootstrap, ou inferência causal em séries temporais — áreas em que se recomenda literatura específica;
5. **Links e versões.** URLs de documentação podem mudar; links quebrados são corrigidos no arquivo-fonte da documentação. Onde havia ambiguidade de versão (ex.: `forecast` 9.x, `fable` 0.x), priorizou-se a documentação oficial publicada nas páginas citadas.

---

## 15. Referências

**Materiais acadêmicos e institucionais**

- HYNDMAN, Rob J.; ATHANASOPOULOS, George. *Forecasting: Principles and Practice* (3ª ed.). OTexts, 2021 (atualizado continuamente). Online gratuito: [otexts.com/fpp3](https://otexts.com/fpp3/). Capítulos usados: [decomposição](https://otexts.com/fpp3/decomposition.html), [padrões/gráficos](https://otexts.com/fpp3/tspatterns.html), [ARIMA](https://otexts.com/fpp3/arima.html), [toolbox/acurácia](https://otexts.com/fpp3/toolbox.html);
- MORETTIN, Pedro A.; TOLOI, Clélia M. C. *Análise de Séries Temporais* (Vol. 1, 3ª ed.; Vol. 2). ABE-Projeto Fisher/Blucher, 2018–2020. Página dos autores com dados e erratas: [users.ime.usp.br/~pam/ST.html](https://users.ime.usp.br/~pam/ST.html);
- AUSTRALIAN BUREAU OF STATISTICS. *Time Series Analysis: The Basics*. [abs.gov.au](https://www.abs.gov.au/websitedbs/d3310114.nsf/home/time+series+analysis:+the+basics);
- NAU, Robert. *ARIMA models for time series forecasting* (notas de aula, Duke University) — [identificação de ordens via ACF/PACF](https://people.duke.edu/~rnau/411arim3.htm);
- NCSS. *The Box-Jenkins Method* (PDF, capítulo 470). [ncss.com](https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/NCSS/The_Box-Jenkins_Method.pdf);
- BOBEICA, Elena; HARTWIG, Benny. *The COVID-19 shock and challenges for time series models*. ECB Working Paper n. 2558, 2021. [ecb.europa.eu](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2558~22b223a7c6.en.pdf);
- WIKIPÉDIA. *Augmented Dickey–Fuller test* — [link](https://en.wikipedia.org/wiki/Augmented_Dickey%E2%80%93Fuller_test); *Spurious relationship* — [link](https://en.wikipedia.org/wiki/Spurious_relationship);
- VIGEN, Tyler. *Spurious Correlations*. [tylervigen.com](https://www.tylervigen.com/spurious-correlations).

**Documentação dos pacotes R e tutoriais**

- FORECAST PACKAGE (Hyndman). [pkg.robjhyndman.com/forecast](https://pkg.robjhyndman.com/forecast/); [`auto.arima` — documentação](https://www.rdocumentation.org/packages/forecast/versions/9.0.1/topics/auto.arima);
- HYNDMAN, Rob J.; KHANDAKAR, Yeasmin. *Automatic Time Series Forecasting: the forecast Package for R*. Journal of Statistical Software, 27(3), 2008. [jstatsoft.org](https://www.jstatsoft.org/index.php/jss/article/view/v027i03);
- FABLE/FABLETOOLS (tidyverts). [*Accuracy* — referência](https://fabletools.tidyverts.org/reference/accuracy.html); anúncio do pacote [mitchelloharawild.com/blog/fable](https://mitchelloharawild.com/blog/fable/);
- HYRMAN, Rob J. *Cross-validation for time series* (blog Hyndsight) — [tscv](https://robjhyndman.com/hyndsight/tscv/) e [tscv-fable](https://robjhyndman.com/hyndsight/tscv-fable/);
- HYNDMAN, Rob J. *Another look at forecast-accuracy metrics for intermittent demand* (Foresight, 2006) — [robjhyndman.com/papers/foresight.pdf](https://robjhyndman.com/papers/foresight.pdf);
- HYNDMAN, Rob J. *Degrees of freedom for the Ljung-Box test* — [hyndsight](https://robjhyndman.com/hyndsight/ljung_box_df.html);
- R-STATISTICS.CO. *Time Series Objects in R: ts, xts, zoo, or tsibble?* — [link](https://r-statistics.co/Time-Series-Objects-in-R.html); *Prophet in R* — [link](https://r-statistics.co/Prophet-in-R.html); *ARIMA Diagnostics in R: checkresiduals() and Ljung-Box* — [link](https://r-statistics.co/ARIMA-Diagnostics-in-R.html);
- ANALYSE MACRO (Brasil). *Modelos ARIMA e Metodologia Box-Jenkins* — [link](https://analisemacro.com.br/econometria-e-machine-learning/modelos-arima-e-metodologia-box-jenkins/); *Estacionariedade de séries temporais* — [link](https://analisemacro.com.br/estatistica-e-econometria/estacionariedade-de-series-temporais/); curso *[Séries Temporais em R](https://analisemacro.com.br/wp-content/uploads/2015/02/CursoCompleto.pdf)* (Ricardo Lima);
- METRICGATE. *Spurious Regression: When Trends Fake a Relationship* — [link](https://metricgate.com/blogs/spurious-regression-nonstationary-series/).

**Dados brasileiros (Banco Central)**

- PORTAL DE DADOS ABERTOS DO BCB. Recurso `json_serie-sgs-11` (Taxa Selic) com a especificação do endpoint do SGS — [link](https://dadosabertos.bcb.gov.br/dataset/11-taxa-de-juros---selic/resource/b73edc07-bbac-430c-a2cb-b1639e605fa8);
- WILSONFREITAS. *python-bcb — documentação SGS* (padrão de uso da API; válido como referência de formato das chamadas) — [link](https://wilsonfreitas.github.io/python-bcb/sgs.html);
- RBCB: pacote R de acesso ao BCB — [GitHub](https://github.com/wilsonfreitas/rbcb) (funções `get_series`, `get_market_expectations`);
- GETBCBDATA: pacote R de acesso ao SGS — [GitHub](https://github.com/msperlin/GetBCBData) (função `gbcbd_get_series`).

---

*Documento encerrado. Comentários, correções e sugestões são bem-vindos — especialmente se você rodar os exemplos e encontrar divergências de versão.*
