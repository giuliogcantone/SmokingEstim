library(broom)

results =
  tibble(
    n.nuclei = params$n.nuclei,
    E_peers = params$E_peers,
    p_D = params$p_D,
    w = params$w,
    gamma = params$gamma,
    r = params$r,
    N = stats$N,
    f.smokers = stats$f.smokers,
    n.gold = stats$n.gold,
    f.smokers_gold = stats$f.smokers_gold,
    n.0 = stats$n.0,
    f.smokers_0 = stats$f.smokers_0,
    n.hp = stats$n.hp,
    f.smokers_hp = stats$f.smokers_hp,
    f.smokers_hp_m = stats$f.smokers_hp_m,
    n.hy = stats$n.hy,
    f.smokers_hy = stats$f.smokers_hy,
    f.smokers_hy_m = stats$f.smokers_hy_m,
  ) %>%
  mutate(
    ERR_gold = f.smokers - f.smokers_gold,
    ERR_0 = f.smokers - f.smokers_0,
    ERR_hp = f.smokers - f.smokers_hp,
    ERR_hp_m = f.smokers - f.smokers_hp_m,
    ERR_hy = f.smokers - f.smokers_hy,
    ERR_hy_m = f.smokers - f.smokers_hy_m,
  ) %>%
  mutate(
    gamma_level = cut_number(gamma,4)
  )

results$ERR_gold %>% boxplot()
results$ERR_0 %>% boxplot()
results$ERR_hp %>% boxplot()
results$ERR_hp_m %>% boxplot()

lm(
  data = results,
  scale(abs(ERR_hp)) ~ 0 + scale(r) + scale(gamma) +
    scale(N) + scale(f.smokers) + scale(n.hp) + scale(n.gold/n.hp) +
    scale(w) + scale(E_peers)
) %>% tidy()

lm(
  data = results,
  scale(abs(ERR_hy)) ~ 0 + scale(r) + scale(gamma) +
    scale(n.hp) + scale(n.hp - n.gold) + scale(w) + scale(E_peers)
) %>% tidy()

results %>% group_by(gamma_level) %>%
  summarise(
    ERR_gold = median(abs(ERR_gold)),
    ERR_0 = median(abs(ERR_0)),
    ERR_hp = median(abs(ERR_hp)),
    ERR_hy = median(abs(ERR_hy))
  )
