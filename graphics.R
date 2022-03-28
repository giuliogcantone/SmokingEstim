### PMF

library(VGAM)
library(tidyverse)

ggplot(tibble(
  z = seq(.05,.95,.1) %>% as.character() %>% str_remove("0"),
  x = 0:9,
  y = dbetabinom(0:9,9,.225,.3)),
  aes(x=z,y=y)) +
  geom_histogram(stat="identity") +
  theme_classic(base_size = 20) +
  ylab("Probability")+
  xlab("Propensity score")
  geom_vline(aes(xintercept = 2.75), color = "red")
  annotate("text",x=4,y=.32,label="Expected Value",
           color = "red", size = 10)
ggsave("PMF.eps")
  
(rbetabinom(10000,9,.225,.3)+.05) %>% mean()
(rbetabinom(10000,9,.225,.3)+.05) %>% var()

### MIXING MATRIX
library(Matrix)

map_dfc(
  seq(1,10,1),
  function(.x) {
    (dnorm(seq(1,10,1),.x,.x))}
) |> as.matrix() |>
  forceSymmetric(uplo = "U") |>
  as.matrix() -> B

B/sum(B) -> B
colnames(B) <- (seq(.05,.95,.1))
rownames(B) <- (seq(.05,.95,.1))

(B^8) / sum(B^8) -> G
(B^2) / sum(B^2) -> G

sum(G)
               
colnames(G) <- (seq(.05,.95,.1))
rownames(G) <- (seq(.05,.95,.1))

superheat::superheat(G,
                     X.text = round(G, 3),
                     X.text.size = 5,
                     legend = FALSE,
                     X.text.col = "white")
superheat::superheat(B,
                     X.text = round(B, 3),
                     X.text.size = 5,
                     legend = FALSE,
                     X.text.col = "black")
ggsave("B.png")