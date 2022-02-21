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

###

latent <- sbm(
  n = 100,
  B = diag(30),
  pi= (rep(1,30))
)

Pmatrix(c(rep(c(1,rep(0,9))),9),1, nrow = 10) %>% View()

sbm

RGfas

graph %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble %>% filter(name == "1")

graph %>%
  activate(edges) %>%
  filter(
    from == 2 | to == 2
  ) %>% as_tibble()


library(ggplot2)
ggplot(data = data.frame(x = c(0, 1)), mapping = aes(x = x)) +
  stat_function(fun = dbeta, args = c(3,1), n = 100) +
  geom_histogram(aes(x = rbeta(100,3,1)))

ggplot() +
  geom_histogram(aes(x = rbeta(100,3,1), y = ..density..), bins = 30) +
  stat_function(fun = dbeta, args = c(6,1), n = 100)
