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
  xlab("Propensity score")
  geom_vline(aes(xintercept = 2.75), color = "red")
  annotate("text",x=4,y=.32,label="Expected Value",
           color = "red", size = 10)

(rbetabinom(10000,9,.225,.3)+.05) %>% mean()
(rbetabinom(10000,9,.225,.3)+.05) %>% var()

### MIXING
library(Matrix)

map_dfc(
  seq(1,10,1),
  function(.x) {
    (dpois(seq(1,10,1),.x))}
  ) |> as.matrix() |>
  forceSymmetric() |>
  as.matrix() -> B

B/sum(B) -> B

colnames(B) <- (seq(.05,.95,.1))
rownames(B) <- (seq(.05,.95,.1))

superheat::superheat(B,
                     X.text = round(B, 3),
                     X.text.size = 5,
                     legend = FALSE,
                     X.text.col = "black")


(B^.2 / sum(B^.2)) %>%
  superheat::superheat(.,
                       X.text = round(., 3),
                       X.text.size = 5,
                       legend = FALSE,
                       X.text.col = "black")


### TOY NETWORK
library(igraph)
library(tidyverse)

# Graph1
do.call(
  graph.disjoint.union,
  lapply(
    rpois(10,1.2)+1,
    function(k) {
      make_full_graph(k) %>% set_vertex_attr(
        name = "nucleus_size", value = k)
    }
  )
) %>%
  set_vertex_attr(name = "nucleus", value = membership(components(.))) %>% 
  set_vertex_attr(name = "name", value = seq(vcount(.))) %>%
  set_edge_attr(name = "color", value = "red") %>% as_tbl_graph() %>%
  activate(nodes) %>% mutate(name_ = sample(name)) -> graph1

tbl_graph(graph1 %>% activate(nodes)
          %>% as_tibble() %>% mutate(name = name_) %>%
            select(-name_),graph1 %>%
            activate(edges) %>%
            as_tibble() %>%
            mutate(name = from) %>%
            inner_join(graph1 %>% activate(nodes)
                       %>% as_tibble() %>%
                         select(name,name_), by = "name") %>%
            mutate(from = name_) %>% select(-name_) %>%
            mutate(name = to) %>%
            inner_join(graph1 %>% activate(nodes)
                       %>% as_tibble() %>%
                         select(name,name_), by = "name") %>%
            mutate(to = name_) %>% select(-c(name_,name)),
          directed = F
          ) -> graph1

graph1 %>% plot()


### Graph2
map_dfc(
  seq(1,10,1),
  function(.x) {
    (dpois(seq(1,10,1),.x))}
) |> as.matrix() |>
  forceSymmetric() |>
  as.matrix() -> B

B/sum(B) -> B

colnames(B) <- (seq(.05,.95,.1))
rownames(B) <- (seq(.05,.95,.1))

superheat::superheat(B,
                     X.text = round(B, 3),
                     X.text.size = 5,
                     legend = FALSE,
                     X.text.col = "black")

(B^.2 / sum(B^.2)) %>%
  superheat::superheat(.,
                       X.text = round(., 3),
                       X.text.size = 5,
                       legend = FALSE,
                       X.text.col = "black")

(B^.2 / sum(B^.2)) -> G

#
library(fastRG)
library(tidygraph)

latent <- fastRG::sbm(
  n = V(graph1) %>% length(),
  pi = dbetabinom(0:9,9,.225,.3),
  B = G,
  expected_degree = 20,
  sort_nodes = T
)

dbetabinom(0:9,9,.225,.3)

graph2 <- sample_tidygraph(latent) %>%
  activate(edges) %>%
  mutate(
    color = "blue"
  )

graph2 %>%
  activate(nodes) %>%
  mutate(
    beta = latent$z
  ) %>% as_tibble() %>% mutate(
    beta = abs((str_remove(beta,"block") %>%
      as.integer()) - 10)/10 + .05
  ) -> nodes

nodes %>% group_by(beta) %>% count()
nodes %>% count()

### Union

graph1 %>% as_tbl_graph() %>%
  activate(nodes) %>% as_tibble() %>% View()

graph1 %>% as_tbl_graph() %>%
  activate(nodes) %>%
  group_by(nucleus) %>%
  mutate(alpha = rbetabinom(1,9,.225,.3),
         alpha = (alpha/10) + .05) %>% ungroup() %>%
  mutate(
    beta = latent$z,
    beta = abs((str_remove(beta,"block") %>%
                         as.integer()) - 10)/10 + .05,
    name = as.character(name) %>%
      bind_edges(graph2%>%activate(edges)%>%as_tibble)
    
  ) %>%
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
