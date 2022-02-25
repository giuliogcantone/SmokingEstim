library(tidyverse)
library(fastRG)
library(tidygraph)
library(igraph)


### Small component

do.call(
  graph.disjoint.union,
  lapply(
    rpois(5000,1.2)+1,
    function(k) {
      make_full_graph(k) %>% set_vertex_attr(
        name = "nucleus_size", value = k)
    }
  )
) %>%
  set_vertex_attr(name = "nucleus", value = membership(components(.))) %>% 
  set_vertex_attr(name = "name", value = seq(vcount(.))) -> graph
  
graph %>%
as_tbl_graph(directed = FALSE) %>% activate(nodes) %>%
  as_tibble() %>% group_by(nucleus_size) %>% count()

graph %>%
  as_tbl_graph(directed = FALSE) %>% activate(nodes) %>%
  as_tibble() %>% count()

###
library(VGAM)

ggplot() + xlim(1,9) +
  geom_function(fun = dbetabinom,
                args = list(size = 8,
                            prob = .2,
                            rho = .9))

dbetabinom(c(0:9),8,.2) %>% plot()
dbetabinom(c(0:9),8,.2,.1) %>% plot()

rbetabinom(100,8,.2,.05) %>% qplot()



ggplot() + geom_function(fun = dbeta,
                  args = list(
                    shape1 =5.1,
                    shape2 =4.9),
                  color = "red")

dbeta(1)

?dbetabinom

ggplot() + geom_function(fun = dbeta,
                         args = list(
                           shape1 =8,
                           shape2 =8))

ggplot() + geom_function(fun = dnorm,
                         args = list(
                           mean = .5,
                           sd = .5))

dbeta(seq(.01,.99,.01),.5,.5)

map_dfc(
  seq(.01,.99,.01),
  function(.x) {
    (dbeta(seq(.01,.99,.01),.x*10,10-.x*10))/
      sum(dbeta(seq(.01,.99,.01),.x*10,10-.x*10))}
) %>% as.matrix() -> B

map_dfc(
  seq(.01,.99,.01),
  function(.x) {
    (pbeta(seq(.01,.99,.01),.x*10,10-.x*10))/
      sum(dbeta(seq(.01,.99,.01),.x*10,10-.x*10))}
) %>% as.matrix() -> p

map_dfc(
  seq(.1,.9,.1),
  function(.x) {
    (pbeta(seq(.1,.9,.1),.x*10,10-.x*10))/
      sum(dbeta(seq(.1,.9,.1),.x*10,10-.x*10))}
) %>% as.matrix() %>% forceSymmetric() %>% as.matrix() -> P2

View(B)

superheat::superheat(forceSymmetric(B) %>% as.matrix())

B %>% View()
forceSymmetric(p) %>% as.matrix() %>% View()

B %*%t(B) %>% sum()
B * t(B) %>% sum()

B*t(B) %>% View()

B == C

View(B)
View(C)

(B*t(B))/2 -> C

(B+t(B)) %>% View()

(B+t(B)) %>% sum()

B %>% sum()
t(B) %>% sum()

transp

?forceSymmetric

superheat::superheat(forceSymmetric(B) %>% as.matrix())
superheat::superheat(B)
superheat::superheat(B %*%t (B))
superheat::superheat(B *t(B))
superheat::superheat(forceSymmetric(p) %>% as.matrix())
superheat::superheat(P2)

forceSymmetric(p) %>% as.matrix() %>% View()

sum(C)

map_dfc(
  seq(.01,.98,.01),
  function(.x) {
    dbeta(seq(.01,.99,.01),.x*10,10-.x*10)}
) %>% as.matrix() %>% View()

isSymmetric(B)

B %>% View()



latent <- fastRG::sbm(
  n = 10090,
  pi = dbetabinom(1:99,99,.2,.165),
  B = forceSymmetric(p),
  expected_degree = 20
)


dbetabinom(c(1:9),9,.2,.165) %>% plot()

dbetabinom(c(1:9),9,.35,.165) %>% plot()
rbetabinom(10090,9,.2,.165) %>% .[. %in% c(1:9)] %>% mean()

(rbetabinom(10090,8,.15,.2)+1) %>% qplot()

(rbetabinom(10090,8,.15,.5)+1) %>% mean()

rbetab

%>% mean()

latent <- fastRG::sbm(
  n = 10090,
  pi = dbetabinom(0:8,8,.2,.2),
  B = forceSymmetric(P2),
  expected_degree = 120
)

latent$z %>% table()
latent$pi

graph2 <- sample_tidygraph(latent)

graph2 %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble() -> nodes2

nodes2$name

graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes2, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
         ) %>% select(-block) %>% inner_join(nodes2, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% count() %>% View()

  
latent$pi

graph2 %>%
  activate(nodes) %>%
  as_tibble() %>% count()


library(ggplot2)
ggplot(data = data.frame(x = c(0, 1)), mapping = aes(x = x)) +
  stat_function(fun = dbeta, args = c(3,1), n = 100) +
  geom_histogram(aes(x = rbeta(100,3,1)))

ggplot() +
  geom_histogram(aes(x = rbeta(100,3,1), y = ..density..), bins = 30) +
  stat_function(fun = dbeta, args = c(6,1), n = 100)

mvtnorm::rmvnorm()

rmvnorm()
rmvnorm(100,c(1,1))

dmvnorm(c(1,1),c(1,1))

rmvnorm(1000,c(0,0),
        matrix(c(1,2,2,5), nrow=2)) %>% as_tibble() %>%
  ggplot() +
  stat_density2d_filled(aes(x = V1, y = V2))


dmvnorm(c(1,1),c(5,5),
        matrix(c(2.5,10,10,2.5), nrow=2))



matrix(c(2.5,1,2.5,1), nrow=2)
