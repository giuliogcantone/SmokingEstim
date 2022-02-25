library(tidyverse)
library(fastRG)
library(tidygraph)
library(igraph)

### Nucleotic component sampling
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
  set_vertex_attr(name = "name", value = seq(vcount(.))) -> graph1

### B

map_dfc(
  seq(.1,.9,.1),
  function(.x) {
    (dbeta(seq(.1,.9,.1),.x*10,10-.x*10))/
      sum(dbeta(seq(.1,.9,.1),
                .x*10,
                10-.x*10))}
) %>% as.matrix() -> B

###repeat 3 times
(B+mean(B))/sum((B+mean(B))) -> B

###
forceSymmetric(B) %>% as.matrix() -> B

# checks

superheat:superheat(B)
View(B)

B[1,] %>% sum()
B[3,] %>% sum()
B[5,] %>% sum()

### latent sbm
library(VGAM)

latent <- fastRG::sbm(
  n = V(graph1) %>% length(),
  pi = dbetabinom(0:8,8,.2,.2),
  B = B,
  expected_degree = 20,
  sort_nodes = T
)

### Block component sampling
graph2 <- sample_tidygraph(latent)

graph2 %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble() -> nodes

### Checks

graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
  ) %>% select(-block) %>% inner_join(nodes, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% count() -> links

links %>% View()

graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
  ) %>% select(-block) %>% inner_join(nodes, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% group_by(from_block) %>% count() %>%
  rename(to = n) %>%
  add_column(graph2 %>%
               activate(edges) %>%
               as_tibble() %>% mutate(name = as.character(from)) %>%
               inner_join(nodes, by = "name") %>%
               mutate(from = name,
                      from_block = block,
                      name = as.character(to)
               ) %>% select(-block) %>% inner_join(nodes, by = "name") %>%
               rename(to_block = block) %>%
               group_by(from_block, to_block) %>% group_by(to_block) %>%
               count()) %>% ungroup %>% transmute(block = from_block,
                                                  n_links = (to+n)/2) %>%
  inner_join(nodes %>%
               group_by(block) %>% count(), by = "block")


graph2 %>%
  activate(edges) %>%
  as_tibble() %>% mutate(name = as.character(from)) %>%
  inner_join(nodes, by = "name") %>%
  mutate(from = name,
         from_block = block,
         name = as.character(to)
  ) %>% select(-block) %>% inner_join(nodes, by = "name") %>%
  rename(to_block = block) %>%
  group_by(from_block, to_block) %>% count() %>% View()

# making the join

graph1 %>% as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>%
  mutate(name = as.character(name)) %>%
  inner_join(graph2 %>% activate(nodes) %>%
               as_tibble(), by = "name")

graph2 