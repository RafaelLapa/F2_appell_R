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
