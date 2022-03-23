library(broom)

results =
  tibble(
    n.nuclei = params$n.nuclei,
    E_peers = params$E_peers,
    p_D = params$p_D,
    w = params$w,
    gamma = params$gamma,
    phi = stats$phi,
    r = params$r,
    N = stats$N,
    f.smokers = stats$f.smokers,
    n.gold = stats$n.gold,
    f.smokers_gold = stats$f.smokers_gold,
    n.0 = stats$n.0,
    f.smokers_0 = stats$f.smokers_0,
    n.hp = stats$n.hp,
    f.smokers_hp = stats$f.smokers_hp,
    n.hy = stats$n.hy,
    f.smokers_hy = stats$f.smokers_hy
  ) %>%
  mutate(
    ERR_gold = f.smokers_gold - f.smokers,
    ERR_0 = f.smokers_0 - f.smokers,
    ERR_hp = f.smokers_hp - f.smokers,
    ERR_hy =  f.smokers_hy - f.smokers
  ) %>%
  mutate(
    gamma_level = cut_number(gamma,4),
    phi_level = cut_number(phi,4),
    w_level = cut_number(w,4),
    n_level = cut_number(n.gold,4)
  )


results %>% summarise(
  cor = cor(gamma,phi),
  cor2 = cor(w,phi)
)

results %>% ggplot() +
  geom_density(aes(ERR_hp), color = "green") +
  geom_density(aes(ERR_hy), color = "red") +
  geom_density(aes(ERR_0)) +
  geom_density(aes(ERR_gold), color = "blue") +
  facet_wrap(~gamma_level)+
  theme_classic()

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

results %>% summarise(
  var_gold = var(ERR_gold),
  var_0 = var(ERR_0),
  var_hp = var(ERR_hp),
  var_hy = var(ERR_hy)
)

results %>%
  mutate(id = cur_group_rows(), .before = "n.nuclei") %>%
  filter(abs(ERR_gold) > .06) %>% pull(id) -> redo

results$ERR_gold %>% min() 
