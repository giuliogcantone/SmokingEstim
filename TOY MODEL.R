library(tidyverse)
library(igraph)
library(tidygraph)
library(VGAM)

#g1
n = 46000

tibble(k = rpois(n,1.2)+1,
       nucleus = 1:n) -> g1

tibble(
  w = .4,
  nucleus = rep(g1$nucleus,g1$k)  %>% sample()) %>%
  mutate(name = cur_group_rows()) %>%
  group_by(nucleus) %>%
  mutate(alpha = rbetabinom(1,9,.225,.3),
         alpha = (alpha/10) + .05) %>%
  ungroup() -> g1

g1 %>% with(
    .,
    do.call(
      rbind,
      Map(
        function(v) {
          get.data.frame(
            set_vertex_attr(
              make_full_graph(length(v)),
              name = "name", value = v
            )
          )
        },
        split(name, nucleus)
      )
    )
  ) %>% as_tibble() %>% mutate(color = "red") %>%
  tbl_graph(g1,.,directed = F) %>%
  activate(nodes) %>% arrange(name) -> g1

# g2
library(Matrix)

map_dfc(
  seq(1,10,1),
  function(.x) {
    (dnorm(seq(1,10,1),.x,.x))}
) |> as.matrix() |>
  forceSymmetric(uplo = "U") |>
  as.matrix() -> B

B/sum(B) -> B
(B^.5 / sum(B^.5)) -> G

#
library(fastRG)
library(tidygraph)

latent <- fastRG::sbm(
  n = V(g1) %>% length(),
  pi = dbetabinom(0:9,9,.225,.3),
  B = G,
  expected_degree = 20,
  sort_nodes = T
)

sample_tidygraph(latent) %>%
  activate(nodes) %>% mutate(name = as.integer(name)) -> g2
                             
graph_join(g1,g2) %>% as_tbl_graph() %>%
  activate(edges) %>% mutate(
    color = ifelse(is.na(color),"blue","red")
  ) %>%
  activate(nodes) %>%
  mutate(
    beta = latent$z,
    beta = abs((str_remove(beta,"block") %>%
                  as.integer()) - 10)/10 + .05
  ) %>% group_by(name) %>%
  mutate(
    e = (alpha*w)+beta*(1-w),
    smoker = rbinom(1,1,e),
    color = ifelse(smoker==1,"green","pink")) %>%
  activate(edges) %>%
  distinct() %>%
  activate(nodes) %>% arrange(name) -> g

g %>% activate(nodes) %>% as_tibble() %>%
  group_by(smoker) %>% count() %>%
  ungroup %>% mutate(f = n/sum(n))

library(ggraph)
g %>% ggraph(layout = "stress") +
  geom_edge_link(aes(colour = color,
                    width = color)) +
  geom_node_point(aes(fill = color),
                  size = 7, shape = 21, stroke = 1, color = 'black') +
  scale_fill_manual(values = c("pink" = "white","green" = "forestgreen")) +
  scale_edge_color_manual(values = c("red" = "red","blue" = "orange"))+
  geom_node_text(aes(label = name))+
  scale_edge_width_manual(values = c(
    "red" = 1,
    "blue" = .5)) +
  theme_void() +
  theme(legend.position="none") 

ggsave("network.eps", device = "eps")


### Diagnostica

g %>% activate(nodes) %>%
  mutate(k = centrality_degree()) %>%
  as_tibble() %>% group_by(smoker) %>%
  summarize('E(k)' = mean(k))
