library(tidyverse)
library(igraph)
library(tidygraph)
library(fastRG)
library(VGAM)
library(Matrix)
library(fixest)
library(lme4)
library(broom.mixed)

# params

params = list(n.nuclei = NA,
              E_peers = NA,
              w = NA,
              gamma = NA,
              r = NA,
              G = NA,
              B = map_dfc(
                seq(1,10,1),
                function(.x) {
                  (dnorm(seq(1,10,1),.x,.x))}
              ) |> as.matrix() |>
                forceSymmetric(uplo = "U") |>
                as.matrix())

stats = list(
  g1 = NA,
  g2 = NA,
  g = NA,
  N = NA,
  mean.e = NA,
  f.smokers = NA,
  n.gold = NA,
  f.smokers_gold = NA,
  n.0 = NA,
  f.smokers_0 = NA,
  n.hp = NA,
  f.smokers_hp = NA,
  f.smokers_hp_m = NA,
  n.hy = NA,
  f.smokers_hy = NA,
  f.smokers_hy_m = NA
)

for(i in 226:250){

params$n.nuclei[i] = sample(5000:15000,1,T)
params$E_peers[i] = sample(3:23,1,T)
params$w[i] = runif(1,.1,.5)
params$gamma[i] = runif(1,.2,.8)
params$r[i] = runif(1,0,.5)

params$G[i] = list(params$B^params$gamma[i] /
                     sum(params$B^params$gamma[i]))


# g1

tibble(k = rpois(params$n.nuclei[i],1.2)+1,
       nucleus = 1:params$n.nuclei[i]) -> stats$g1

tibble(
  w = params$w[i],
  nucleus = rep(stats$g1$nucleus,
                stats$g1$k)  %>% sample()) %>%
  mutate(name = cur_group_rows()) %>%
  group_by(nucleus) %>%
  mutate(alpha = rbetabinom(1,9,.225,.3),
         alpha = (alpha/10) + .05) %>%
  ungroup() -> stats$g1

stats$g1 %>% with(
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
  tbl_graph(stats$g1,.,directed = F) %>%
  activate(nodes) %>% arrange(name) -> stats$g1

# g2

latent <- fastRG::sbm(
  n = V(stats$g1) %>% length(),
  pi = dbetabinom(0:9,9,.225,.3),
  B = params$G[[i]],
  expected_degree = 20,
  sort_nodes = T
)

sample_tidygraph(latent) %>%
  activate(nodes) %>% mutate(name = as.integer(name)) -> stats$g2

# UNION

graph_join(stats$g1,stats$g2) %>% as_tbl_graph() %>%
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
  activate(nodes) %>% arrange(name) -> stats$g

stats$g %>% as_tibble %>%
  summarise(n = n()) %>%
  pull(n) -> stats$N[i]


stats$g %>% as_tibble %>%
  summarise(e = mean(e)) %>%
  pull(e) -> stats$mean.e[i]

stats$g %>% as_tibble %>%
  summarise(f = mean(smoker)) %>%
  pull(f) -> stats$f.smokers[i]

stats$g %>%
  mutate(
    m = local_size(mindist = 1),
    ego = local_members(mindist = 1)
  ) -> stats$g

### Recruitment

stats$mean.e[i] <- stats$g %>% as_tibble() %>%
  summarise(e = mean(e)) %>% pull(e)

t = 0

stats$g %>% as_tibble() %>% sample_n(1000/(1-params$r[i])) %>%
  rowwise() %>%
  mutate(
    r = rbinom(1,100, params$r[i] + (e/10 - stats$mean.e[i]/10)*.25),
    r = r/100,
    fill = rbinom(1,1,(1-r))
  ) %>% filter(fill == 1) %>% pull(name) -> stage

stats$g %>% group_by(name) %>%
  mutate(
    crawl_stage = if_else(name %in% stage, 0, NA_real_),
    parent = if_else(crawl_stage == 0, 0, NA_real_),
    root = if_else(crawl_stage == 0, name, NA_integer_),
    k = if_else(crawl_stage == 0, rpois(1,.5), NA_integer_),
  ) -> stats$g

stats$g %>% as_tibble() %>% filter(crawl_stage == 0) %>%
  summarise(n = n()) %>%
  pull(n) -> stats$n.gold[i]

stats$g %>% as_tibble() %>% filter(crawl_stage == 0) %>%
  summarise(f = mean(smoker)) %>%
  pull(f) -> stats$f.smokers_gold[i]

stats$g %>% group_by(name) %>%
  mutate(
    crawl_stage = ifelse(crawl_stage == 0,
                         sample(c(0,NA),1), NA)) -> stats$g

stats$g %>% as_tibble() %>% filter(crawl_stage == 0) %>%
  summarise(n = n()) %>%
  pull(n) -> stats$n.0[i]

stats$g %>% as_tibble() %>% filter(crawl_stage == 0) %>%
  summarise(f = mean(smoker)) %>%
  pull(f) -> stats$f.smokers_0[i]

# Poisson

while(stats$g %>% as_tibble %>%
      filter(k>0,
             crawl_stage == t) %>% count() %>%
      as.integer() != 0) {
  
  stats$g %>% as_tibble %>%
    filter(crawl_stage == t,
           k>0) %>%
    select(name,ego,k,m,root) %>%
    rowwise() %>%
    mutate(picks = list(sample(ego,min(k,m)))) %>%
    select(root, parent = name, name = picks) %>%
    unnest(cols = c(name)) -> stage
  
  stage %>% anti_join(stats$g %>%
                        as_tibble %>%
                        filter(crawl_stage <= t), by = "name") %>%
    left_join(stats$g %>% as_tibble() %>% select(name,e),
              by = "name") %>%
    group_by(name) %>%
    mutate(
      r = rbinom(1,100, params$r[i] + (e/10 - stats$mean.e[i]/10)*.25),
      r = r/100,
      fill = rbinom(1,1,(1-r))) %>%
    filter(fill == 1)-> stage
  
  t = t+1
  
  stats$g %>% group_by(name) %>%
    mutate(crawl_stage = ifelse(name %in% stage$name,t,crawl_stage),
           parent = ifelse(name %in% stage$name,
                           stage$parent[stage$name == name],
                           parent),
           root = ifelse(name %in% stage$name,
                         stage$root[stage$parent == parent],
                         root),
           k = ifelse(crawl_stage == t,
                      rpois(1,.5),
                      k)
    ) %>% ungroup() -> stats$g
}

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  summarise(n = n()) %>%
  pull(n) -> stats$n.hp[i]

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  summarise(f = mean(smoker)) %>%
  pull(f) -> stats$f.smokers_hp[i]

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  group_by(root) %>% mutate(
    n.root = n()
  ) %>% lmer(data = .,
             smoker ~ 1 + 1|n.root) %>% tidy() %>%
  pull(estimate) %>% .[1] -> stats$f.smokers_hp_m[i]


# Yule 

stats$g %>% mutate(
  crawl_stage = ifelse(crawl_stage > 0,NA,0)
) -> stats$g

t = 0

while(stats$g %>% as_tibble %>%
      filter(k>0,
             crawl_stage == t) %>% count() %>%
      as.integer() != 0) {
  
  stats$g %>% as_tibble %>%
    filter(crawl_stage == t,
           k>0) %>%
    select(name,ego,k,m,root) %>%
    rowwise() %>%
    mutate(picks = list(sample(ego,min(k,m-1)))) %>%
    select(root, parent = name, name = picks) %>%
    unnest(cols = c(name)) -> stage
  
  stage %>% anti_join(stats$g %>%
                        as_tibble %>%
                        filter(crawl_stage <= t), by = "name") %>%
    left_join(stats$g %>% as_tibble() %>% select(name,e),
              by = "name") %>%
    group_by(name) %>%
    mutate(
      r = rbinom(1,100, params$r[i] + (e/10 - stats$mean.e[i]/10)*.25),
      r = r/100,
      fill = rbinom(1,1,(1-r))) %>%
    filter(fill == 1)-> stage
  
  t = t+1
  
  stats$g %>% group_by(name) %>%
    mutate(crawl_stage = ifelse(name %in% stage$name,t,crawl_stage),
           parent = ifelse(name %in% stage$name,
                           stage$parent[stage$name == name],
                           parent),
           root = ifelse(name %in% stage$name,
                         stage$root[stage$parent == parent],
                         root),
           k = ifelse(crawl_stage == t,
                      ryules(1,3)-1,
                      k)
    ) %>% ungroup() -> stats$g
}

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  summarise(n = n()) %>%
  pull(n) -> stats$n.hy[i]

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  summarise(f = mean(smoker)) %>%
  pull(f) -> stats$f.smokers_hy[i]

stats$g %>% as_tibble() %>% filter(crawl_stage > -1) %>%
  group_by(root) %>% mutate(
    n.root = n()
  ) %>% lmer(data = .,
             smoker ~ 1 + 1|n.root) %>% tidy() %>%
  pull(estimate) %>% .[1] -> stats$f.smokers_hy_m[i]
}
