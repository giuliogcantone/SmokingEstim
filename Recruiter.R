### PROB SAMPLING
library(tidygraph)

stage = list()
sample = 1/

g %>% activate(nodes) %>%
  sample_n(10) %>% pull(name) -> stage[[1]]

g %>% activate(nodes) %>%
  sample_n(10) %>% pull(name) -> stage[[2]]

mutate(g %>% activate(nodes),
         Sample_Stage = ifelse(name %in% stage[[1]],1,0))

stage %>% as.numeric()


###

stats$g %N>% as_tibble() %>% View()

t = 0

stats$g %>% as_tibble %>%
  filter(crawl_stage == t,
         k>0) %>%
  select(name,ego,k,root) %>%
  rowwise() %>%
  mutate(picks = list(sample(ego,k))) %>%
  select(root, parent = name, name = picks) %>%
  unnest(cols = c(name)) -> stage

stage %>% anti_join(stats$g %>%
                      as_tibble %>%
                      filter(crawl_stage <= t), by = "name") %>%
  left_join(stats$g %>% as_tibble() %>% select(name,e),
       by = "name") %>%
  group_by(name) %>%
  mutate(
    r = rbinom(1,100, params$r[1] + (e/10 - stats$mean.e[1]/10)*.25),
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
                       root[name == parent],
                       root),
         k = ifelse(crawl_stage == t,
                    rpois(1,.5),
                    k)
  ) %>% ungroup() -> stats$g


###

g %>% convert(to_local_neighborhood,
        node = which(.N()$name == "1"),
        order = 1,
        mode = "all")


rtruncpareto(20,0,1,1)
rtruncpareto(1,1,1000,.5)

rpareto(1,1,1)

ryules(100,6) %>% mean()

30*.05

rpois(1000000,.5) %>% sum()
(ryules(10000,3)-1) %>% sum()
(ryules(10000,2)-1) %>% sum()


####
library(tidyverse)
library(tidygraph)
set.seed(1)
play_erdos_renyi(n = 100000,m =1000000) %>%
  mutate(name = cur_group_rows(),
         ego = local_members(mindist = 1),
         m = local_size()) -> g

g

g %>% as_tibble() %>% sample_n(1000) %>%
    pull(name) -> stage

g %>% rowwise() %>%
  mutate(
    attrition = rbinom(1,100,
                       .75 +
                         e/10 -
                         (g %>% as_tibble %>%
                         summarise(e = mean(e)) %>%
                           pull(e)
                         )
                       ))


,
    crawl_stage = ifelse(name %in% stage,
                              0,NA),
         parent = NA,
         k = ifelse(crawl_stage == 0,
                         rpois(1,.5),
                         0)
             ) %>% ungroup() -> g

t = 0

g %>% as_tibble %>% View()

while(g %>% as_tibble %>%
      filter(k>0,
             crawl_stage == t) %>% count() %>%
      as.integer() != 0) {

g %>% as_tibble %>%
  filter(crawl_stage == t,
         k>0) %>%
  select(name,ego,k) %>% rowwise() %>%
  mutate(picks = list(sample(ego,k))) %>%
  select(name = picks, parent = name) %>%
  unnest(cols = c(name)) -> stage
  
stage %>% anti_join(g %>%
                      as_tibble %>%
                      filter(crawl_stage <= t), by = "name") -> stage
  
t = t+1

g %>% group_by(name) %>%
  mutate(crawl_stage = ifelse(name %in% stage$name,t,crawl_stage),
         parent = ifelse(name %in% stage$name,
                         stage$parent[stage$name == name],
                         parent),
         k = ifelse(crawl_stage == t,
                    rpois(1,.5),
                    k)
  ) %>% ungroup() -> g}

g %>% as_tibble() %>% View()


g %>% as_tibble %>%
  filter(k>0,
         crawl_stage == t) %>% count() %>%
  as.integer()

t


###


tibble = tibble(a = rbinom(100000,100,.5)) %>%
  mutate(b = a/100 - mean(a/100),
         d = rbinom(100000,100,.08 + b*.1)
  )

tibble %>% summarise(mean.d = mean(d),
                     r = cor(a,d))

tibble
