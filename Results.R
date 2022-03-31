library(tidyverse)
library(broom)

results =
  tibble(
    n.nuclei = params$n.nuclei,
    k = params$E_peers,
    w = params$w,
    gamma = params$gamma,
    phi_y = stats$phi,
    phi_k = stats$phi_m,
    Q = stats$Q,
    r = params$r,
    N = stats$N,
    y = stats$f.smokers,
    n_gold = stats$n.gold,
    y_gold = stats$f.smokers_gold,
    n_0 = stats$n.0,
    y_0 = stats$f.smokers_0,
    n_p = stats$n.p,
    y_p = stats$f.smokers_p,
    n_hp = stats$n.hp,
    y_hp = stats$f.smokers_hp,
    n_yu = stats$n.y,
    y_yu = stats$f.smokers_y,
    n_hyu = stats$n.hy,
    y_hyu = stats$f.smokers_hy,
  ) %>%
  mutate(
    err_gold = y_gold - y,
    err_p = y_p - y,
    err_hp = y_hp - y,
    err_yu = y_yu - y,
    err_hyu = y_hyu - y,
    err_0 = y_0 - y,
    abs_gold = abs(err_gold),
    abs_p = abs(err_p),
    abs_hp = abs(err_hp),
    abs_yu = abs(err_yu),
    abs_hyu = abs(err_hyu),
    abs_0 = abs(y_0 - y),
    delta_p = (abs_gold - abs_p)/y,
    delta_yu = (abs_gold - abs_yu)/y,
    delta_hp = (abs_0 - abs_hp)/y,
    delta_hyu = (abs_0 - abs_hyu)/y
  ) %>% mutate(
    gamma_l = cut_number(gamma,4,
                         c("Low",
                           "Mid-Low",
                           "Mid-High",
                           "High")),
    ub_p = err_p + .01,
    ub_yu = err_yu + .01,
    ub_hp = err_hp + .01,
    ub_hyu = err_hyu + .01,
    delta_ubp = (abs_gold - abs(ub_p)) / y,
    delta_ubyu = (abs_gold - abs(ub_yu)) / y,
    delta_ubhp = (abs_gold - abs(ub_hp)) / y,
    delta_ubhyu = (abs_gold - abs(ub_hyu)) / y
  )

library(ggdensity)
library(cowplot)
ggdraw() +
    draw_plot(
      results %>% ggplot(aes(x=phi_y,
                             y=phi_k)) +
        geom_hdr() +
        ylab(bquote(varphi[k])) +
        xlab(bquote(varphi[y])) +
        theme_linedraw(base_size = 30)+
        theme(legend.position = "top")
    ) +
  draw_image(magick::image_read("gamma.png"),
             scale = .25, x = .4, y = -.1)
ggsave("phis.png")

results %>%
  select(k,w,gamma,phi_y,phi_k,r,N,y,Q,abs_gold) %>%
  names() %>%
  paste('scale(abs_p)~',.)%>%
  map_df(~tidy(lm(as.formula(.x), 
                  data= results %>% select(-gamma_l) %>% mutate_all(scale)
  ))) %>%
  filter(!term %>% str_detect("Intercept")) %>%
  arrange(p.value) %>% select(-statistic) %>%
  mutate_if(is.numeric,round,3) %>%
  kableExtra::kable("latex", booktabs = T) %>%
  kableExtra::kable_styling(latex_options = c("scale_down"))

results %>%
  select(k,w,gamma,phi_y,phi_k,Q,r,N,y,abs_gold) %>%
  names() %>%
  paste('scale(abs_hp)~',.)%>%
  map_df(~tidy(lm(as.formula(.x), 
                  data= results %>% select(-gamma_l) %>% mutate_all(scale)
  ))) %>%
  filter(!term %>% str_detect("Intercept"))%>%
  arrange(p.value) %>% select(-statistic) %>%
  mutate_if(is.numeric,round,3) %>%
  kableExtra::kable("latex", booktabs = T) %>%
  kableExtra::kable_styling(latex_options = c("scale_down"))

results %>%
  select(k,w,gamma,phi_y,phi_k,Q,r,N,y,abs_gold) %>%
  names() %>%
  paste('scale(abs_yu)~',.)%>%
  map_df(~tidy(lm(as.formula(.x), 
                  data= results %>% select(-gamma_l) %>% mutate_all(scale)
  ))) %>%
  filter(!term %>% str_detect("Intercept"))%>%
  arrange(p.value) %>% select(-statistic) %>%
  mutate_if(is.numeric,round,3) %>%
  kableExtra::kable("latex", booktabs = T) %>%
  kableExtra::kable_styling(latex_options = c("scale_down"))

results %>%
  select(k,w,gamma,phi_y,phi_k,Q,r,N,y,abs_gold) %>%
  names() %>%
  paste('scale(abs_hyu)~',.)%>%
  map_df(~tidy(lm(as.formula(.x), 
                  data= results %>% select(-gamma_l) %>%
                    mutate_all(scale)
  ))) %>%
  filter(!term %>% str_detect("Intercept")) %>%
  arrange(p.value) %>% select(-statistic) %>%
  mutate_if(is.numeric,round,3)


###

results %>%
  select(w,gamma,r,N,y,abs_gold) %>%
  names() %>%
  paste('scale(abs_p)~',.)%>%
  map_df(~tidy(lm(as.formula(.x), 
                  data= results %>% select(-gamma_l) %>% mutate_all(scale)
  ))) %>%
  filter(!term %>% str_detect("Intercept")) %>%
  arrange(p.value) %>% select(term, estimate) %>%
  mutate_if(is.numeric,round,3) %>%
  left_join(
    results %>%
      select(w,gamma,r,N,y,abs_gold) %>%
      names() %>%
      paste('scale(abs_yu)~',.)%>%
      map_df(~tidy(lm(as.formula(.x), 
                      data= results %>% select(-gamma_l) %>% mutate_all(scale)
      ))) %>%
      filter(!term %>% str_detect("Intercept")) %>%
      arrange(p.value) %>% select(term, estimate) %>%
      mutate_if(is.numeric,round,3), by = "term"
  ) %>%
  left_join(
    results %>%
      select(w,gamma,r,N,y,abs_gold) %>%
      names() %>%
      paste('scale(abs_hp)~',.)%>%
      map_df(~tidy(lm(as.formula(.x), 
                      data= results %>% select(-gamma_l) %>% mutate_all(scale)
      ))) %>%
      filter(!term %>% str_detect("Intercept")) %>%
      arrange(p.value) %>% select(term, estimate) %>%
      mutate_if(is.numeric,round,3), by = "term"
  ) %>%
  left_join(
    results %>%
      select(w,gamma,r,N,y,abs_gold) %>%
      names() %>%
      paste('scale(abs_hyu)~',.)%>%
      map_df(~tidy(lm(as.formula(.x), 
                      data= results %>% select(-gamma_l) %>% mutate_all(scale)
      ))) %>%
      filter(!term %>% str_detect("Intercept")) %>%
      arrange(p.value) %>% select(term, estimate) %>%
      mutate_if(is.numeric,round,3), by = "term"
  ) %>% add_column(.after = "term", Concept = c(
    "Stage 0", "Homophily", "Familism",
    "Attrition", "Target Quota", "Pop. Size"
  )) %>%
  rename(I = estimate.x,
         II = estimate.y,
         III = estimate.x.x,
         IV = estimate.y.y) %>%
  kableExtra::kable("latex")

###

results %>% summarise(
  cor = cor(gamma,phi_k),
  cor2 = cor(gamma,phi_y),
  cor3 = cor(phi_k,phi_y)
)

results %>% ggplot() +
  geom_density(aes(err_gold), color = "goldenrod", size = 2) +
  #geom_density(aes(err_0), size = 2) +
  geom_density(aes(err_p), color = "purple", size = 2) +
  xlab("Error") +
  ylab("") +
  facet_wrap(~gamma_l)+
  theme_classic(base_size = 20)

###

results %>%
  mutate(id = cur_group_rows(), .before = "n.nuclei") %>%
  filter(abs(err_gold) > .055) %>% pull(id) -> redo

results$err_gold %>% max()


results %>% select(c(starts_with("err"),gamma_l)) %>%
  pivot_longer(cols = starts_with("err"),
               names_to = "Sampling") %>%
  ggplot() +
  geom_boxplot(aes(x = Sampling,
                   y = value))+
  facet_wrap(~gamma_l)
  

results$err_gold