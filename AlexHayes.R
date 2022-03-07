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
  seq(.05,.95,.1),
  function(.x) {
    (dbeta(seq(.05,.95,.1),.x*10,10-.x*10))/
      sum(dbeta(seq(.05,.95,.1),
                .x*10,
                10-.x*10))}
) %>% as.matrix() -> B

map_dfc(
  seq(.05,.95,.1),
  function(.x) {
    (pbeta(seq(.05,.95,.1),.x*10,10-.x*10))}
) %>% as.matrix() -> B

forceSymmetric(B) %>% as.matrix() -> B
B/sum(B) -> B
superheat::superheat(B)

View(B)
  
### latent sbm
library(VGAM)
library(tidyverse)

(18+36) / 2

dbetabinom(0:9,9,.27,.5) %>% barplot()

dbetabinom(0:9,9,.27,.3) %>% barplot()
dbetabinom(0:9,9,.28,.3) %>% barplot()
dbetabinom(0:9,9,.35,.3) %>% barplot()

(rbetabinom(10000,9,.22,.3)+.5) %>% mean()

(rbetabinom(10000,9,.22,.3)+.5) %>% var()

param=list()
param$p = .35

mean(((rbetabinom(10000,9,param$p,.3)+.5) +
    (rbetabinom(10000,9,param$p,.3)+.5))/
  2)

var(((rbetabinom(10000,9,param$p,.3)+.5) +
        (rbetabinom(10000,9,param$p,.3)+.5))/
       2)

latent <- fastRG::sbm(
  n = V(graph1) %>% length(),
  pi = dbetabinom(0:9,9,.22,.18),
  B = B,
  expected_degree = 500,
  sort_nodes = T
)

### Block component sampling
graph2 <- sample_tidygraph(latent)

graph2 %>%
  activate(nodes) %>%
  mutate(
    block = latent$z
  ) %>% as_tibble() -> nodes

graph2 %>%activate(edges) %>%
  sample_n(length(E(graph2)) / 25) -> graph2

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
                                                  n_links = (to+n)) %>%
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
               as_tibble(), by = "name") %>%
  bind_edges(graph2%>%activate(edges)%>%as_tibble) %>%
  activate(nodes) %>%
  select(-nucleus_size) %>%
  mutate(
    nucleus = sample(nucleus, replace = F),
    k = centrality_degree(),
    neighs = local_members()) %>%
  group_by(nucleus) %>%
  mutate(nuclear_prop = (rbetabinom(1,9,.2,.2)+.5)/10) %>%
  group_by(block) %>%
  mutate(block_prop = str_remove(block,"block") %>%
           as.integer(),
         block_prop = abs(block_prop - 10.5)/10) %>%
  ungroup() %>% group_by(name) %>%
  mutate(
    p = mean(c(nuclear_prop,block_prop)),
    smoker = rbinom(1,1,p)) -> pop

pop %>% as_tibble() %>% View()

### DIAGNOSTIC: Homophily in the pop

pop %>% activate(nodes) %>% mutate(name = as.character(name)) %>%
  as_tibble() -> pop_nodes

pop %>% activate(edges) %>%
  as_tibble() %>%
  mutate(from = as.character(from),
         to = as.character(to)) -> pop_edges

pop_nodes %>% group_by(smoker) %>% summarize(k = mean(k))

pop_edges %>% rename(name = from) %>%
  inner_join(pop_nodes %>% select(name,smoker), by = "name") %>%
  rename(from = name,
         smoker_from = smoker,
         name = to) %>%
  inner_join(pop_nodes %>% select(name,smoker), by = "name") %>%
  rename(to = name,
         smoker_to = smoker) %>%
  mutate(smoke_sum = smoker_to + smoker_from) %>%
  group_by(smoke_sum) %>% count() %>%
  ungroup %>% mutate(f = n/sum(n))

# Diagnostic: Degree connectivity

pop %>% activate(nodes) %>% as_tibble() %>%
  mutate(p = as.character(p)) %>% group_by(p) %>% summarise(
    n(),
    mean(k)
    ) %>% View()

### pop

pop %>% activate(nodes) %>% ungroup() %>%
  mutate(smok.neighs = sum(
           (name %in% neighs))) %>% as_tibble() %>%
  View()