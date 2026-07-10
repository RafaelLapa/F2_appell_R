library(ggplot2)

grid_pontos <- expand.grid(x = seq(-1.5, 1.5, length.out = 400),
                           y = seq(-1.5, 1.5, length.out = 400))

grid_pontos$Convergencia <- with(grid_pontos, abs(x) + abs(y) < 1)

grafico_diamante <- ggplot(grid_pontos, aes(x = x, y = y)) +
  geom_tile(aes(fill = Convergencia)) +
  scale_fill_manual(values = c("FALSE" = "lightgray", "TRUE" = "dodgerblue"), 
                    labels = c("Diverge (Precisa de continuação)", "Converge (Série Original)")) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +
  theme_minimal(base_size = 14) +
  labs(title = "Região de Convergência da Série F2 de Appell",
       subtitle = "Condição Absoluta: |x| + |y| < 1",
       x = "Eixo x", y = "Eixo y", fill = "Status") +
  theme(legend.position = "bottom")

print(grafico_diamante)
ggsave("diamante_convergencia.pdf", plot = grafico_diamante, width = 7, height = 5)




library(microbenchmark)
library(cubature)
library(hypergeo)
library(BAS)
library(dplyr)
library(reticulate)


F2_fixa <- function(a, b1, b2, c1, c2, x, y, max_iter = 100) {
  if (x < 0 || y < 0) stop("variável negativa")
  soma <- 0
  for (m in 0:max_iter) {
    for (n in 0:max_iter) {
      log_termo <- (lgamma(a + m + n) - lgamma(a)) +
        (lgamma(b1 + m) - lgamma(b1)) +
        (lgamma(b2 + n) - lgamma(b2)) -
        (lgamma(c1 + m) - lgamma(c1)) -
        (lgamma(c2 + n) - lgamma(c2)) -
        lgamma(m + 1) - 
        lgamma(n + 1) +
        m * log(x) + 
        n * log(y)
      
      termo <- exp(log_termo)
      soma <- soma + termo
    }
  }
  return(soma)
}


F2_otimizada <- function(x, y, a, b1, b2, c1, c2, max_iter = 500, tol = 1e-12) {
  if (x < 0 || y < 0) stop("Variáveis x e y não podem ser negativas.")
  soma <- 0
  for (m in 0:max_iter) {
    termo_n_anterior <- Inf 
    for (n in 0:max_iter) {
      log_termo <- (lgamma(a + m + n) - lgamma(a)) +
        (lgamma(b1 + m) - lgamma(b1)) +
        (lgamma(b2 + n) - lgamma(b2)) -
        (lgamma(c1 + m) - lgamma(c1)) -
        (lgamma(c2 + n) - lgamma(c2)) -
        lgamma(m + 1) - lgamma(n + 1) +
        m * log(x) + n * log(y)
      
      termo <- exp(log_termo)
      if (is.finite(termo)) soma <- soma + termo
      if (abs(termo) > termo_n_anterior && n > 1) break 
      if (abs(termo) < tol) break
      termo_n_anterior <- abs(termo)
    }
    if (n == 0 && abs(termo) < tol) break
  }
  return(soma)
}


F2_integral <- function(a, b1, b2, c1, c2, x, y) {
  integrando <- function(v) {
    u1 <- v[1]; u2 <- v[2]
    if(u1 == 0 | u1 == 1 | u2 == 0 | u2 == 1) return(0)
    return((u1^(b1 - 1) * (1 - u1)^(c1 - b1 - 1) * u2^(b2 - 1) * (1 - u2)^(c2 - b2 - 1)) / ((1 - u1*x - u2*y)^a))
  }
  const <- (gamma(c1)*gamma(c2)) / (gamma(b1)*gamma(b2)*gamma(c1-b1)*gamma(c2-b2))
  res <- hcubature(integrando, lowerLimit = c(0,0), upperLimit = c(1,1), tol = 1e-6)
  return(const * res$integral)
}


x_bench <- 0.25; y_bench <- 0.30


F2_fixa(1, 1, 1, 2, 2, 0.1, 0.2, max_iter = 100)
F2_otimizada(0.1, 0.2, 1, 1, 1, 2, 2)
F2_integral(1, 1, 1, 2, 2, 0.1, 0.2)


x_bench <- 0.25; y_bench <- 0.30

testes_tempo_R <- microbenchmark(
  `Grade Fixa O(N²)` = F2_fixa(1, 1, 1, 2, 2, x_bench, y_bench, max_iter = 100),
  `Série Otimizada` = F2_otimizada(x_bench, y_bench, 1, 1, 1, 2, 2),
  `Integração (hcubature)` = F2_integral(1, 1, 1, 2, 2, x_bench, y_bench),
  times = 50
)

grafico_bench1 <- ggplot(testes_tempo_R, aes(x = expr, y = time / 1e6, fill = expr)) +
  geom_violin(alpha = 0.8, trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.6) +
  coord_flip() +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  scale_fill_manual(values = c(
    "Grade Fixa O(N²)" = "tomato",
    "Série Otimizada" = "darkorange",
    "Integração (hcubature)" = "seagreen"
  )) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribuição do Tempo de Execução",
       subtitle = "Comparação da evolução algorítmica (50 iterações avaliadas)",
       x = "Método Implementado em R",
       y = "Tempo de Execução (Milissegundos - escala logarítmica)") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

print(grafico_bench1)
ggsave("benchmark_plot_metodos.pdf", plot = grafico_bench1, width = 8, height = 5)


tabela_R <- summary(testes_tempo_R, unit = "ms")
print(tabela_R)
write.csv2(tabela_R, "tabela_benchmark_R.csv", row.names = FALSE)



calc_erro <- function(estimado, real) {
  erro <- abs(estimado - real) / abs(real) * 100
  return(round(erro, 4))
}

py_run_string("
import mpmath
mpmath.mp.dps = 15

def appellf2_custom(a, b1, b2, c1, c2, x, y, max_iter=150):
    soma = mpmath.mpf(0)
    for m in range(max_iter):
        for n in range(max_iter):
            num = mpmath.rf(a, m+n) * mpmath.rf(b1, m) * mpmath.rf(b2, n)
            den = mpmath.rf(c1, m) * mpmath.rf(c2, n) * mpmath.fac(m) * mpmath.fac(n)
            termo = (num / den) * (mpmath.mpf(x)**m) * (mpmath.mpf(y)**n)
            
            soma += termo
            if abs(termo) < 1e-15:
                break
    return float(soma)
")

set.seed(42)
pontos_teste <- data.frame(
  id = 1:20,
  x = runif(20, -0.4, 0.4), 
  y = runif(20, -0.4, 0.4)
)

a <- 1; b1 <- 1; b2 <- 1; c1 <- 2; c2 <- 2

resultados <- pontos_teste %>%
  rowwise() %>%
  mutate(
    Valor_R = F2_integral(a, b1, b2, c1, c2, x, y),
    Valor_Python = py$appellf2_custom(a, b1, b2, c1, c2, x, y),
    Erro_Absoluto = abs(Valor_R - Valor_Python),
    Erro_Relativo = Erro_Absoluto / abs(Valor_Python)
  ) %>%
  ungroup()

grafico_precisao <- ggplot(resultados, aes(x = Valor_Python, y = Valor_R)) +
  geom_point(color = "dodgerblue", size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "Validação de Precisão Numérica: R vs. Python",
       subtitle = "A linha tracejada a vermelho representa a igualdade perfeita (Erro = 0)",
       x = "Estimativa Python (mpmath)",
       y = "Estimativa R (Integração hcubature)") +
  annotate("text", x = min(resultados$Valor_Python), y = max(resultados$Valor_R), 
           label = paste("Erro Máximo Observado:", format(max(resultados$Erro_Relativo), scientific = TRUE)), 
           hjust = 0, vjust = 1, size = 4, fontface = "italic")

print(grafico_precisao)
ggsave("precisao_python.pdf", plot = grafico_precisao, width = 7, height = 5)


print(head(resultados))
write.csv2(resultados, "tabela_erros_precisao.csv", row.names = FALSE)



testes_tempo_Py <- microbenchmark(
  `Integração (Nativo R)` = F2_integral(a, b1, b2, c1, c2, x_bench, y_bench),
  `Série mpmath (Python)` = py$appellf2_custom(a, b1, b2, c1, c2, x_bench, y_bench),
  times = 30 
)

grafico_bench2 <- ggplot(testes_tempo_Py, aes(x = expr, y = time / 1e6, fill = expr)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.5) +
  coord_flip() + 
  scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
  scale_fill_manual(values = c("Integração (Nativo R)" = "seagreen", "Série mpmath (Python)" = "darkorange")) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribuição do Tempo de Execução",
       subtitle = "Comparação da rapidez computacional (30 iterações avaliadas)",
       x = "Método Algorítmico",
       y = "Tempo de Execução (Milissegundos - escala logarítmica)") +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none" 
  )

print(grafico_bench2)
ggsave("benchmark_plot_python.pdf", plot = grafico_bench2, width = 8, height = 5)


tabela_Py <- summary(testes_tempo_Py, unit = "ms")
print(tabela_Py)
write.csv2(tabela_Py, "tabela_benchmark_python.csv", row.names = FALSE)





F1_otimizada <- function(x, y, a, b1, b2, c, max_iter = 500, tol = 1e-12) {
  if (x < 0 || y < 0) stop("Variáveis x e y não podem ser negativas.")
  soma <- 0
  for (m in 0:max_iter) {
    termo_n_anterior <- Inf 
    for (n in 0:max_iter) {
      log_termo <- (lgamma(a + m + n) - lgamma(a)) +
        (lgamma(b1 + m) - lgamma(b1)) +
        (lgamma(b2 + n) - lgamma(b2)) -
        (lgamma(c + m + n) - lgamma(c)) -
        lgamma(m + 1) - lgamma(n + 1) +
        m * log(x) + n * log(y)
      
      termo <- exp(log_termo)
      if (is.finite(termo)) soma <- soma + termo
      if (abs(termo) > termo_n_anterior && n > 1) break 
      if (abs(termo) < tol) break
      termo_n_anterior <- abs(termo)
    }
    if (n == 0 && abs(termo) < tol) break
  }
  return(soma)
}

F1_integral <- function(a, b1, b2, c, x, y) {
  integrando <- function(u) {
    if (u == 0 | u == 1) return(0)
    return( u^(a - 1) * (1 - u)^(c - a - 1) * (1 - u*x)^(-b1) * (1 - u*y)^(-b2) )
  }
  const <- gamma(c) / (gamma(a) * gamma(c - a))
  res <- integrate(Vectorize(integrando), lower = 0, upper = 1, rel.tol = 1e-6)
  return(const * res$value)
}



F3_otimizada <- function(x, y, a1, a2, b1, b2, c, max_iter = 500, tol = 1e-12) {
  if (x < 0 || y < 0) stop("Variáveis x e y não podem ser negativas.")
  soma <- 0
  for (m in 0:max_iter) {
    termo_n_anterior <- Inf 
    for (n in 0:max_iter) {
      log_termo <- (lgamma(a1 + m) - lgamma(a1)) +
        (lgamma(a2 + n) - lgamma(a2)) +
        (lgamma(b1 + m) - lgamma(b1)) +
        (lgamma(b2 + n) - lgamma(b2)) -
        (lgamma(c + m + n) - lgamma(c)) -
        lgamma(m + 1) - lgamma(n + 1) +
        m * log(x) + n * log(y)
      
      termo <- exp(log_termo)
      if (is.finite(termo)) soma <- soma + termo
      if (abs(termo) > termo_n_anterior && n > 1) break 
      if (abs(termo) < tol) break
      termo_n_anterior <- abs(termo)
    }
    if (n == 0 && abs(termo) < tol) break
  }
  return(soma)
}

F3_integral <- function(a1, a2, b1, b2, c, x, y) {
  integrando <- function(t) {
    t1 <- t[1]; t2 <- t[2]
    if(t1 == 0 | t1 == 1 | t2 == 0 | t2 == 1) return(0)
    u <- t1
    v <- (1 - t1) * t2
    jacobiano <- 1 - t1
    kernel <- u^(b1 - 1) * v^(b2 - 1) * (1 - u - v)^(c - b1 - b2 - 1) * 
      (1 - u*x)^(-a1) * (1 - v*y)^(-a2)
    return(kernel * jacobiano)
  }
  const <- gamma(c) / (gamma(b1) * gamma(b2) * gamma(c - b1 - b2))
  res <- hcubature(integrando, lowerLimit = c(0,0), upperLimit = c(1,1), tol = 1e-6)
  return(const * res$integral)
}



F4_otimizada <- function(x, y, a, b, c1, c2, max_iter = 500, tol = 1e-12) {
  if (x < 0 || y < 0) stop("Variáveis x e y não podem ser negativas.")
  soma <- 0
  for (m in 0:max_iter) {
    termo_n_anterior <- Inf 
    for (n in 0:max_iter) {
      log_termo <- (lgamma(a + m + n) - lgamma(a)) +
        (lgamma(b + m + n) - lgamma(b)) -
        (lgamma(c1 + m) - lgamma(c1)) -
        (lgamma(c2 + n) - lgamma(c2)) -
        lgamma(m + 1) - lgamma(n + 1) +
        m * log(x) + n * log(y)
      
      termo <- exp(log_termo)
      if (is.finite(termo)) soma <- soma + termo
      if (abs(termo) > termo_n_anterior && n > 1) break 
      if (abs(termo) < tol) break
      termo_n_anterior <- abs(termo)
    }
    if (n == 0 && abs(termo) < tol) break
  }
  return(soma)
}



x_bench <- 0.25; y_bench <- 0.30


a_F1 <- 1.5; b1_F1 <- 1.2; b2_F1 <- 1.2; c_F1 <- 3.0
a1_F3 <- 1.2; a2_F3 <- 1.2; b1_F3 <- 1.5; b2_F3 <- 1.5; c_F3 <- 4.0
a_F4 <- 1.5; b_F4 <- 1.5; c1_F4 <- 2.5; c2_F4 <- 2.5


bench_F1 <- microbenchmark(
  `Série Otimizada` = F1_otimizada(x_bench, y_bench, a_F1, b1_F1, b2_F1, c_F1),
  `Integração (1D)` = F1_integral(a_F1, b1_F1, b2_F1, c_F1, x_bench, y_bench),
  times = 30
)
print(summary(bench_F1, unit = "ms"))


bench_F3 <- microbenchmark(
  `Série Otimizada` = F3_otimizada(x_bench, y_bench, a1_F3, a2_F3, b1_F3, b2_F3, c_F3),
  `Integração (hcubature)` = F3_integral(a1_F3, a2_F3, b1_F3, b2_F3, c_F3, x_bench, y_bench),
  times = 30
)
print(summary(bench_F3, unit = "ms"))


bench_F4 <- microbenchmark(
  `Série Otimizada` = F4_otimizada(x_bench, y_bench, a_F4, b_F4, c1_F4, c2_F4),
  times = 30
)
print(summary(bench_F4, unit = "ms"))



py_run_string("
import mpmath
mpmath.mp.dps = 15

def appellf1_custom(a, b1, b2, c, x, y, max_iter=150):
    soma = mpmath.mpf(0)
    for m in range(max_iter):
        for n in range(max_iter):
            num = mpmath.rf(a, m+n) * mpmath.rf(b1, m) * mpmath.rf(b2, n)
            den = mpmath.rf(c, m+n) * mpmath.fac(m) * mpmath.fac(n)
            termo = (num / den) * (mpmath.mpf(x)**m) * (mpmath.mpf(y)**n)
            soma += termo
            if abs(termo) < 1e-15:
                break
    return float(soma)

def appellf3_custom(a1, a2, b1, b2, c, x, y, max_iter=150):
    soma = mpmath.mpf(0)
    for m in range(max_iter):
        for n in range(max_iter):
            num = mpmath.rf(a1, m) * mpmath.rf(a2, n) * mpmath.rf(b1, m) * mpmath.rf(b2, n)
            den = mpmath.rf(c, m+n) * mpmath.fac(m) * mpmath.fac(n)
            termo = (num / den) * (mpmath.mpf(x)**m) * (mpmath.mpf(y)**n)
            soma += termo
            if abs(termo) < 1e-15:
                break
    return float(soma)

def appellf4_custom(a, b, c1, c2, x, y, max_iter=150):
    soma = mpmath.mpf(0)
    for m in range(max_iter):
        for n in range(max_iter):
            num = mpmath.rf(a, m+n) * mpmath.rf(b, m+n)
            den = mpmath.rf(c1, m) * mpmath.rf(c2, n) * mpmath.fac(m) * mpmath.fac(n)
            termo = (num / den) * (mpmath.mpf(x)**m) * (mpmath.mpf(y)**n)
            soma += termo
            if abs(termo) < 1e-15:
                break
    return float(soma)
")


set.seed(42)
pontos_teste <- data.frame(
  id = 1:20,
  x = runif(20, 0.05, 0.35),
  y = runif(20, 0.05, 0.35)
)


resultados_F1 <- pontos_teste %>%
  rowwise() %>%
  mutate(
    Valor_R = F1_integral(a_F1, b1_F1, b2_F1, c_F1, x, y),
    Valor_Python = py$appellf1_custom(a_F1, b1_F1, b2_F1, c_F1, x, y),
    Erro_Relativo = abs(Valor_R - Valor_Python) / abs(Valor_Python)
  ) %>% ungroup()

max(resultados_F1$Erro_Relativo)

grafico_F1 <- ggplot(resultados_F1, aes(x = Valor_Python, y = Valor_R)) +
  geom_point(color = "darkgrey", size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "Validação de Precisão F1: R vs. Python",
       subtitle = "Cruzamento da Integral (R) contra Série Arbitrária (mpmath)",
       x = "Estimativa Python (mpmath)", y = "Estimativa R (Integral)") +
  annotate("text", x = min(resultados_F1$Valor_Python), y = max(resultados_F1$Valor_R), 
           label = paste("Erro Máx:", format(max(resultados_F1$Erro_Relativo), scientific = TRUE)), 
           hjust = 0, vjust = 1, fontface = "italic")
print(grafico_F1)
ggsave("grafico_F1.pdf", plot = grafico_F1, width = 7, height = 5)



resultados_F3 <- pontos_teste %>%
  rowwise() %>%
  mutate(
    Valor_R = F3_integral(a1_F3, a2_F3, b1_F3, b2_F3, c_F3, x, y),
    Valor_Python = py$appellf3_custom(a1_F3, a2_F3, b1_F3, b2_F3, c_F3, x, y),
    Erro_Relativo = abs(Valor_R - Valor_Python) / abs(Valor_Python)
  ) %>% ungroup()

max(resultados_F3$Erro_Relativo)

grafico_F3 <- ggplot(resultados_F3, aes(x = Valor_Python, y = Valor_R)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "Validação de Precisão F3: R vs. Python",
       subtitle = "Cruzamento da Quadratura Adaptativa (R) contra Série (mpmath)",
       x = "Estimativa Python (mpmath)", y = "Estimativa R (hcubature)") +
  annotate("text", x = min(resultados_F3$Valor_Python), y = max(resultados_F3$Valor_R), 
           label = paste("Erro Máx:", format(max(resultados_F3$Erro_Relativo), scientific = TRUE)), 
           hjust = 0, vjust = 1, fontface = "italic")
print(grafico_F3)
ggsave("grafico_F3.pdf", plot = grafico_F3, width = 7, height = 5)




set.seed(42)
pontos_teste_F4 <- data.frame(
  id = 1:20,
  x = runif(20, 0.01, 0.15),
  y = runif(20, 0.01, 0.15)
)

resultados_F4 <- pontos_teste_F4 %>%
  rowwise() %>%
  mutate(
    Valor_R = F4_otimizada(x, y, a_F4, b_F4, c1_F4, c2_F4),
    Valor_Python = py$appellf4_custom(a_F4, b_F4, c1_F4, c2_F4, x, y),
    Erro_Relativo = abs(Valor_R - Valor_Python) / abs(Valor_Python)
  ) %>% ungroup()

max(resultados_F4$Erro_Relativo)

grafico_F4 <- ggplot(resultados_F4, aes(x = Valor_Python, y = Valor_R)) +
  geom_point(color = "seagreen", size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "Validação de Precisão F4: R vs. Python",
       subtitle = "Cruzamento da Série Otimizada (R) contra Série Arbitrária (mpmath)",
       x = "Estimativa Python (mpmath)", y = "Estimativa R (Série Otimizada)") +
  annotate("text", x = min(resultados_F4$Valor_Python), y = max(resultados_F4$Valor_R), 
           label = paste("Erro Máx:", format(max(resultados_F4$Erro_Relativo), scientific = TRUE)), 
           hjust = 0, vjust = 1, fontface = "italic")
print(grafico_F4)
ggsave("grafico_F4.pdf", plot = grafico_F4, width = 7, height = 5)





cenarios <- list(
  `Cenario 1: (2.0, 3.0, 0.0)`    = c(2, 3, 0),
  `Cenario 2: (2.5, 3.5, 0.2)`    = c(2.5, 3.5, 0.2),
  `Cenario 3: (1.0, 1.0, 0.0)`    = c(1, 1, 0),
  `Cenario 4: (0.5, 0.5, 0.8)`    = c(0.5, 0.5, 0.8),
  `Cenario 5: (5.0, 5.0, -0.8)`   = c(5, 5, -0.8)
)


f_Z_original <- function(par, x) {
  alpha <- par[1]; beta <- par[2]; rho <- par[3]
  Ki <- function(s, tipo) {
    C   <- (1 + s)^2 * s^(beta - 1) / (gamma(alpha) * gamma(beta))
    nuc <- function(u) u^(alpha + beta - 1) * exp(-(1 + s) * u)
    integrando <- switch(tipo,
                         K1 = function(u) nuc(u),
                         K2 = function(u) nuc(u) * pgamma(u, shape = alpha),
                         K3 = function(u) nuc(u) * pgamma(s * u, shape = beta),
                         K4 = function(u) nuc(u) * pgamma(u, shape = alpha) * pgamma(s * u, shape = beta))
    fat <- c(K1 = 1, K2 = 2, K3 = 2, K4 = 4)[tipo]
    val <- tryCatch(integrate(integrando, 0, Inf, rel.tol = 1e-8, abs.tol = 1e-12)$value, error = function(e) NA_real_)
    return(C * fat * val)
  }
  sapply(x, function(zi) {
    if (zi <= 0 || zi >= 1) return(0)
    s  <- 1 / zi - 1
    k1 <- Ki(s, "K1"); k2 <- Ki(s, "K2"); k3 <- Ki(s, "K3"); k4 <- Ki(s, "K4")
    val <- (1 + rho) * k1 - rho * k2 - rho * k3 + rho * k4
    if (!is.finite(val)) NA_real_ else val
  })
}



f_Z_Appell <- function(par, x) {
  alpha <- par[1]; beta <- par[2]; rho <- par[3]
  Ki_simples <- function(s, tipo) {
    C <- (1 + s)^2 * s^(beta - 1) / (gamma(alpha) * gamma(beta))
    nuc <- function(u) u^(alpha + beta - 1) * exp(-(1 + s) * u)
    integrando <- switch(tipo,
                         K1 = function(u) nuc(u),
                         K2 = function(u) nuc(u) * pgamma(u, shape = alpha),
                         K3 = function(u) nuc(u) * pgamma(s * u, shape = beta))
    fat <- c(K1 = 1, K2 = 2, K3 = 2)[tipo]
    val <- tryCatch(integrate(integrando, 0, Inf, rel.tol = 1e-8, abs.tol = 1e-12)$value, error = function(e) NA_real_)
    return(C * fat * val)
  }
  K4_Appell <- function(s) {
    a_val <- 2 * alpha + 2 * beta
    x_val <- 1 / (2 * (1 + s))
    y_val <- s / (2 * (1 + s))
    log_termo <- beta * log(s) - lgamma(alpha + 1) - lgamma(beta + 1) + lgamma(a_val) - a_val * log(2 * (1 + s))
    log_C <- 2 * log(1 + s) + (beta - 1) * log(s) - lgamma(alpha) - lgamma(beta)
    f2_res <- F2_integral(a = a_val, b1 = 1, b2 = 1, c1 = alpha + 1, c2 = beta + 1, x = x_val, y = y_val)
    return(exp(log_C + log_termo) * 4 * f2_res)
  }
  sapply(x, function(zi) {
    if (zi <= 0 || zi >= 1) return(0)
    s  <- 1 / zi - 1
    k1 <- Ki_simples(s, "K1"); k2 <- Ki_simples(s, "K2"); k3 <- Ki_simples(s, "K3")
    k4 <- tryCatch(K4_Appell(s), error = function(e) NA_real_)
    val <- (1 + rho) * k1 - rho * k2 - rho * k3 + rho * k4
    if (!is.finite(val)) NA_real_ else val
  })
}


z_grid <- seq(0.02, 0.98, length.out = 40)
df_total <- data.frame()

for (nome in names(cenarios)) {
  cat("Processando", nome, "...\n")
  par_atual <- cenarios[[nome]]
  
  pdf_orig <- f_Z_original(par_atual, z_grid)
  pdf_appell <- f_Z_Appell(par_atual, z_grid)
  
  df_temp <- data.frame(
    Cenario = nome,
    Z = rep(z_grid, 2),
    Densidade = c(pdf_orig, pdf_appell),
    Metodo = factor(rep(c("Original (pgamma)", "Appell F2 (hcubature)"), each = length(z_grid)))
  )
  df_total <- bind_rows(df_total, df_temp)
}


df_total$Metodo <- factor(df_total$Metodo, levels = c("Original (pgamma)", "Appell F2 (hcubature)"))

grafico_grid <- ggplot(df_total, aes(x = Z, y = Densidade, 
                                     color = Metodo, linetype = Metodo, 
                                     size = Metodo, alpha = Metodo)) +
  geom_line() +
  facet_wrap(~ Cenario, scales = "free_y", ncol = 2) + 
  scale_color_manual(values = c("Original (pgamma)" = "tomato", "Appell F2 (hcubature)" = "dodgerblue")) +
  scale_linetype_manual(values = c("Original (pgamma)" = "solid", "Appell F2 (hcubature)" = "dashed")) +
  scale_size_manual(values = c("Original (pgamma)" = 2.0, "Appell F2 (hcubature)" = 1.2)) +
  scale_alpha_manual(values = c("Original (pgamma)" = 0.4, "Appell F2 (hcubature)" = 1.0)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Validação Analítica Múltipla: Curvas de Densidade (PDF)",
    subtitle = "Comparação para os 5 cenários iniciais estabelecidos na Tabela 5",
    x = "Domínio de Z (0 a 1)",
    y = "Densidade f(z)"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    legend.key.width = unit(2.5, "cm"),
    strip.background = element_rect(fill = "gray90", color = NA),
    strip.text = element_text(face = "bold")
  )

print(grafico_grid)
ggsave("validacao_pdf_cenarios.pdf", plot = grafico_grid, width = 10, height = 8)
