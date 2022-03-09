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
  ylab("PMF")+
  xlab("Propensity score") +
  geom_vline(aes(xintercept = 2.75), color = "red") +
  annotate("text",x=4,y=.32,label="Expected Value",
           color = "red", size = 10)

### MIXING
library(Matrix)

matrix(rep(.01,100),10) -> N

map_dfc(
  seq(.05,.95,.1),
  function(.x) {
    (pbeta(seq(.05,.95,.1),.x*10,10-.x*10))}
) %>% as.matrix() -> B

forceSymmetric(B) %>% as.matrix() -> B

B/sum(B) -> B
(B^.5)/sum(B^.5) -> BN

mean(B)

mtcars.col <- scale(mtcars) < -0.3
# set all values that satisfy the condition to "white"
mtcars.col <- gsub("TRUE", "white", mtcars.col)
# set all values that do not satisfy the condition to "black"
mtcars.col <- gsub("FALSE", "black", mtcars.col)
# convert to matrix
mtcars.col <- matrix(mtcars.col, ncol = ncol(mtcars))

names(BN)
superheat::superheat(BN,
                     X.text = round(BN, 3),
                     X.text.size = 4,
                     order.cols = 1:10,
                     order.rows = 10:1,
                     legend = FALSE,
                     X.text.col = "goldenrod")

superheat::superheat(B,
                     X.text = round(B, 3),
                     order.cols = 1:10,
                     order.rows = 10:1,
                     X.text.size = 4,
                     X.text.col = "goldenrod")

