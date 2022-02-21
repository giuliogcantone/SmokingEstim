library(tidyverse)
library(fastRG)
library(tidygraph)
library(igraph)


make_full_graph(10) %>%
  set_vertex_attr(name = "name", value = 1) -> graph

?set_vertex_attr

V(graph)
E(graph)

### Small component

make_full_graph(sample(1:5)) %>% set_vertex_attr(
  name = "nucleus",
  value = 1) -> graph

for (i in c(1:3000)){
 graph.disjoint.union(graph,
   make_full_graph(sample(1:5)) %>% set_vertex_attr(
     name = "nucleus",
     value = i+1)) -> graph
}

graph %>% set_vertex_attr(name = "node", value = seq(vcount(.))) %>%
  as_tbl_graph(directed = FALSE)


  
set.seed(1)
do.call(
  graph.disjoint.union,
  lapply(
    1:10,
    make_full_graph
  )
) %>%
  set_vertex_attr(name = "names", value = seq(vcount(.))) %>%
  as_tbl_graph(directed = FALSE)

vcou

sample(c(0:1),replace = T, prob =(.3,.7))

do.call(
  graph.disjoint.union,
  lapply(
    1:3,
    (make_full_graph %>% set_vertex_attr(name = "name",
                                        value = "a"))
  )
) -> small

small %>% as_tbl_graph() 

runif()

set.seed(1)
do.call(
  graph.disjoint.union,
  lapply(
    3,
    make_full_graph
  )
) %>%
  set_vertex_attr(name = "names", value = seq(vcount(.))) %>% 
  as_tbl_graph(directed = FALSE) -> small

E(small) %>% View()

small %>% plot()


?make_full_graph()

set_vert

small$laten

small %>% activate(nodes) %>% as_tibble()

small %>% activate(nodes) %>% as_tibble() %>%
  nrow()
small %>% activate(edges) %>% as_tibble() %>%
  nrow()

small %>% plot()

small %>% activate(edges) %>% as_tibble() %>% View()
small %>% activate(nodes) %>% as_tibble() %>% View()
small %>% plot()


small %>% activate(edges) %>% as_tibble() 
small %>% activate(edges) %>% as_tibble() %>% group_by(from) %>%
  count() %>% group_by(n) %>% count() %>% ungroup() %>% summarize(sum(nn))

###


latent <- sbm(
  n = 100,
  B = diag(30),
  pi= (rep(1,30))
)

Pmatrix(c(rep(c(1,rep(0,9))),9),1, nrow = 10) %>% View()

re

graph <- sample_tidygraph(latent)

graph

small %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble %>% group_by(block) %>% count() %>% group_by(n) %>% count()


?sbm

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
