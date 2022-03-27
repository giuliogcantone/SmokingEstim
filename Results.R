library(broom)

results =
  tibble(
    n.nuclei = params$n.nuclei,
    E_peers = params$E_peers,
    p_D = params$p_D,
    w = params$w,
    gamma = params$gamma,
    phi_y = stats$phi,
    phi_k = stats$phi_m,
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
    n_y = stats$n.y,
    y_y = stats$f.smokers_y,
    n_hy = stats$n.hy,
    y_hy = stats$f.smokers_hy,
  ) %>%
  mutate(
    err_gold = y_gold - y,
    err_p = y_p - y,
    err_hp = y_hp - y,
    err_y = y_y - y,
    err_hy = y_hy - y,
    abs_gold = abs(err_gold),
    abs_p = abs(err_p),
    abs_hp = abs(err_hp),
    abs_y = abs(err_y),
    abs_hy = abs(err_hy)
  ) %>% mutate(
    gamma_l = cut_number(gamma,4)
  )

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

results$phi %>% boxplot()
results$ERR_gold %>% boxplot()
results$ERR_0 %>% boxplot()
results$ERR_hp %>% boxplot()

lm(
  data = results,
  scale(abs(ERR_hp)) ~ 0 + scale(r) + scale(phi) +
    scale(N) + scale(f.smokers) + scale(n.hp) + scale(n.gold/n.hp) +
    scale(w) + scale(E_peers)
) %>% tidy()

lm(
  data = results,
  scale(abs(ERR_hp)) ~ 0 + scale(r) + scale(phi) +
    scale(N) + scale(f.smokers) + scale(n.hp) + scale(n.gold/n.hp)
) %>% tidy()

lm(
  data = results,
  scale(abs(ERR_hy)) ~ 0 + scale(r) + scale(phi) +
    scale(n.hp) + scale(n.hp - n.gold) + scale(w) + scale(E_peers)
) %>% tidy()

results %>% group_by(w_level) %>%
  summarise(
    n =n(),
    ERR_gold = median(abs(ERR_gold)),
    ERR_0 = median(abs(ERR_0)),
    ERR_hp = median(abs(ERR_hp)),
    ERR_hy = median(abs(ERR_hy)),
    n_hp = median(n.hp),
    n_hy = median(n.hy)
  )


results %>%
  mutate(id = cur_group_rows(), .before = "n.nuclei") %>%
  filter(abs(err_gold) > .055) %>% pull(id) -> redo

results$err_gold %>% min()


results$ERR_gold %>% min() 


stats$g %>% as_tibble %>%
  ggplot() +
  stat_bin(aes(e), bins = 10)

stats$g %>% as_tibble %>% pull(e) %>% mean()
