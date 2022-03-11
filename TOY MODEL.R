library(tidyverse)
library(igraph)
library(tidygraph)
library(VGAM)

#g1
n = 46
g1 = list()

tibble(k = rpois(n,1.2)+1,
       nucleus = 1:n) -> g1$n

tibble(
  w = .4,
  nucleus = rep(g1$n$nucleus,g1$n$k)) %>%
  mutate(name = cur_group_rows() %>% sample()) %>%
  group_by(nucleus) %>%
  mutate(alpha = rbetabinom(1,9,.225,.3),
         alpha = (alpha/10) + .05) %>%
  ungroup() -> g1$n

g1$n %>%
  with(
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
  ) %>% as_tibble() %>% mutate(color = "red") -> g1$e

tbl_graph(g1$n,g1$e,directed = F) -> g1


# g2
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

(B^.4 / sum(B^.4)) %>%
  superheat::superheat(.,
                       X.text = round(., 3),
                       X.text.size = 5,
                       legend = FALSE,
                       X.text.col = "black")

(B^.4 / sum(B^.4)) -> G

#
library(fastRG)

latent <- fastRG::sbm(
  n = V(g1) %>% length(),
  pi = dbetabinom(0:9,9,.225,.3),
  B = G,
  expected_degree = 7,
  sort_nodes = T
)

latent$z %>% length()

g2 <- sample_tidygraph(latent) %>%
  activate(nodes) %>%
  mutate(
    beta = latent$z
  ) %>% mutate(
    name = as.integer(name),
    beta = abs((str_remove(beta,"block") %>%
                  as.integer()) - 10)/10 + .05
  ) -> g2

graph_join(g1,g2) %>% as_tbl_graph() %>%
  activate(edges) %>% mutate(
    color = ifelse(is.na(color),"blue","red")
  ) %>%
  activate(nodes) %>% group_by(name) %>%
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
g %>% ggraph(layout = "kk") +
  geom_edge_hive(aes(colour = color,
                    width = color)) +
  geom_node_point(aes(fill = color),
                  size = 7, shape = 21, stroke = 1, color = 'black') +
  scale_fill_manual(values = c("pink" = "pink","green" = "green")) +
  scale_edge_color_manual(values = c("red" = "red","blue" = "blue"))+
  geom_node_text(aes(label = name))+
  scale_edge_width_manual(values = c(
    "red" = 1,
    "blue" = .5)) +
  theme_void()

?scale_edge_size

g %>% activate(nodes) %>% as_tibble() %>% filter(nucleus == nucleus[name == 57])
g %>% activate(edges) %>% as_tibble() %>% filter(to == 57 | from == 57)


g1 %>% activate(nodes) %>% as_tibble() %>% filter(nucleus == nucleus[name == 78])
g1 %>% activate(edges) %>% as_tibble() %>% filter(to == 78 | from == 78)


g1 %>% ggraph(layout = "stress") +
  geom_edge_link(aes(colour = color,
                     width = color)) +
  geom_node_point(size = 7, shape = 21, stroke = 1, color = 'black') +
  scale_edge_color_manual(values = c("red" = "red","blue" = "blue"))+
  geom_node_text(aes(label = name))+
  scale_edge_width_manual(values = c(
    "red" = 1,
    "blue" = .5)) +
  theme_void()

g1 %>% as.igraph()
g1 %>% activate(edges) %>% as_tibble() %>% filter(to == 18 | from == 18)

g1 %>% activate(nodes) %>% arrange(name) %>% as.igraph
g1 %>% activate(nodes)  %>% as.igraph

g1 %>% activate(nodes)
