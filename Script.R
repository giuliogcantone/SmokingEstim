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


### BAD VISUAL
ggplot() + xlim(1,9) +
  geom_function(fun = dbetabinom,
                args = list(size = 8,
                            prob = .2,
                            rho = .9))

### GOOD VISUALS ONLY

dbetabinom(c(0:9),8,.2,.1) %>% plot()

ggplot() + geom_function(fun = dbeta,
                  args = list(
                    shape1 =5.1,
                    shape2 =4.9),
                  color = "red")

ggplot() + geom_function(fun = dbeta,
                         args = list(
                           shape1 =8,
                           shape2 =8))

ggplot() + geom_function(fun = dnorm,
                         args = list(
                           mean = .5,
                           sd = .5))


### GENERAZIONE DELLA MATRICE DI COMMISTIONE

library(superheat)
options(scipen = 999)

map_dfc(
  seq(.1,.9,.1),
  function(.x) {
    (dbeta(seq(.1,.9,.1),.x*10,10-.x*10))/
      sum(dbeta(seq(.1,.9,.1),.x*10,10-.x*10))}
) %>% as.matrix() -> B

B+t(B) -> B
(B+mean(B))/sum((B+mean(B))) -> B

# on B #
View(B)

B[1,] %>% sum()
B[3,] %>% sum()
B[5,] %>% sum()
B[7,] %>% sum()
B[9,] %>% sum()

### latent creation

latent <- fastRG::sbm(
  n = V(graph) %>% length(),
  pi = dbetabinom(0:8,8,.2,.2),
  B = B,
  expected_degree = 200,
  sort_nodes = T
)

latent$z
latent$z %>% table()
latent$pi

graph2 %>% active(nodes) %>%
  mutate_all(slice_sample(prop = .2))

graph2 <- sample_tidygraph(latent)
graph2 %>% activate(edges) %>% as_tibble() %>% slice_sample(prop = .2) -> edges

graph2 %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble() -> nodes

graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes2, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
         ) %>% select(-block) %>% inner_join(nodes2, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% count() -> links

links %>% View()

graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes2, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
  ) %>% select(-block) %>% inner_join(nodes2, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% group_by(from_block) %>% count() %>%
  rename(to = n) %>%
  add_column(graph2 %>%
               activate(edges) %>%
               as_tibble() %>% mutate(name = as.character(from)) %>%
               inner_join(nodes2, by = "name") %>%
               mutate(from = name,
                      from_block = block,
                      name = as.character(to)
               ) %>% select(-block) %>% inner_join(nodes2, by = "name") %>%
               rename(to_block = block) %>%
               group_by(from_block, to_block) %>% group_by(to_block) %>%
               count()) %>% ungroup %>% transmute(block = from_block,
                                      n_links = (to+n)/2) %>%
  inner_join(nodes2 %>%
               group_by(block) %>% count(), by = "block")


graph2 %>% activate(nodes) %>% as_tibble()

E(graph2) %>% length()
