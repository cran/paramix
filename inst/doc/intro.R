## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
suppressPackageStartupMessages({
  library(paramix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(wpp2019)
})

pretty_kable <- function(dt, lims, var) {
  tbl2 <- as.data.frame(
    dcast.data.table(dt, model_partition ~ name, value.var = var)
  )
  tbl2$model_partition <- NULL
  row.names(tbl2) <- paste(
    head(model_agelimits, -1), tail(model_agelimits, -1) - 1L, sep = " - "
  )
  knitr::kable(tbl2)
}

## ----ifrfdef------------------------------------------------------------------
ifr_levin <- function(age_in_years) {
  (10^(-3.27 + 0.0524 * age_in_years))/100
}

## ----echo=FALSE, out.width = "100%", dpi = 600, fig.width = 8, fig.height = 3----
ggplot(
  data.frame(x = 0:100, y = ifr_levin(0:100))
) + aes(x, y) +
  geom_line() + theme_minimal() +
  scale_y_log10("Infection-fatality ratio") +
  labs(x = "Age")

## ----create_pops, cache = TRUE------------------------------------------------
# our model age group cut points
model_agelimits <- c(0, 5, 20, 65, 101)
# get select data from World Population Prospects estimates
data("popF", package = "wpp2019")
data("popM", package = "wpp2019")

pop_dt <- as.data.table(popF)[,
  .(name, age, popF = `2020`)
][
  as.data.table(popM), on = c("name", "age"),
  .(name, age, popF, popM = `2020`)
][, age := as.integer(gsub("^(\\d+)[-+].*$","\\1", age)) ][
  name %like% "Afghanistan|United Kingdom"
]

density_dt <- pop_dt[,
  .(
    from = c(age, max(age)+1),
    weight = c(popF + popM, 0)
  ), by = name
]

rm(popF)
rm(popM)

## ----ifrplotsetup, cache=TRUE-------------------------------------------------
plot_dt <- density_dt[, { # compute parameters for each country of interest
  paramix::parameter_summary(
    f_param = ifr_levin, f_pop = .SD, model_agelimits
  )
}, by = name]

## ----ifrfig, echo = FALSE, cache = TRUE, out.width = "100%", dpi = 600, fig.width = 8, fig.height = 5----
sub_plot_dt <- plot_dt[
  name == "Afghanistan" | method == "wm_f"
][,
  method := factor(fifelse(
    method != "wm_f", as.character(method),
    fifelse(name == "Afghanistan", "af_wm_f", "uk_wm_f")
  ))
][, .SD, .SDcols = -c("name")]
# parameter_summary yields method in (f_val, f_mean, mean_f, wm_f)
# with the exception of wm_f, these are independent of density, so we can drop
# duplicates

pop_p <- ggplot(density_dt[from != max(from)]) + aes(from, weight) +
  geom_col(width = 5, just = 1) +
  facet_grid(cols = vars(name)) +
  scale_x_continuous("Age", expand = expansion()) +
  scale_y_continuous("Population (thousands)") +
  theme_minimal() + theme(panel.spacing.x = unit(1.5, "line"))

ifr_p <- ggplot(sub_plot_dt[x <= density_dt[, max(from)-1]]) + aes(x = x, color = method, y = value) +
  geom_line(data = function(dt) dt[method == "f_val"], lwd = 0.8, lty = "dotted") +
  geom_step(data = function(dt) dt[method != "f_val"]) +
  scale_color_discrete(
    NULL,
    breaks = sub_plot_dt[, value[.N], by = method][order(-V1), method],
    labels = c(
      f_val = "f(age)", f_mid = "f(mid(Age))",
      mean_f = "E[f(age)]", f_mean = "f(E[age])",
      af_wm_f = "AF-weighted E[f(age)]", uk_wm_f = "UK-weighted E[f(age)]"
    )
  ) +
  scale_y_log10("Infection-fatality ratio", expand = expansion()) +
  scale_x_continuous("Age", breaks = 10 * (0:10), expand = expansion()) +
  theme_minimal() +
  theme(
    legend.position = "inside", legend.position.inside = c(0.025, 1),
    legend.justification.inside = c(0, 1)
  )

(pop_p / ifr_p) + plot_layout(heights = c(1, 2)) & theme(text = element_text(size = 10))

## ----params-------------------------------------------------------------------
# setup the model to outcome mapping using `alembic`s
mapping_dt <- density_dt[,
  paramix::alembic(
    f_param = ifr_levin, f_pop = .SD,
    model_partition = model_agelimits,
    output_partition = { res <- seq(min(from), max(from), by = 5L); res[length(res)] <- tail(model_agelimits, 1); res } 
  ),
  by = name
]

params <- mapping_dt[, paramix::blend(.SD), by = name]

## ----paramstbl, echo = FALSE--------------------------------------------------
pretty_kable(params, model_agelimits, "value")

## ----modeldeaths--------------------------------------------------------------
model_density_dt <- density_dt[, .(
  model_partition = model_agelimits[findInterval(from, model_agelimits)],
  weight
), by = name][, .(weight = sum(weight)), by = .(name, model_partition)][,
  weight := weight / sum(weight), by = name
]
model_deaths_dt <- model_density_dt[
  params, on = .(name, model_partition)
][,
  .(name, model_partition, deaths = weight * 1e6 * value)
]

## ----deathtbl, echo = FALSE---------------------------------------------------
pretty_kable(model_deaths_dt, model_agelimits, "deaths")

## ----alembicsetup-------------------------------------------------------------
distill_methods_dt <- model_deaths_dt[,
  paramix::distill_summary(
    mapping_dt[name == .BY],
    .SD[, .(model_partition, value = deaths)]
  ),
  by = name
]

## ----alembicplot, echo = FALSE, warning = FALSE, cache = TRUE, out.width = "100%", fig.width = 8, fig.height = 5, , dpi = 600----
distill_labels <- c(
  f_mean = "Deaths at mean age",
  f_mid = "Deaths uniform across age group",
  wm_f = "paramix approach",
  mean_f = "Deaths proportional to age distribution"
)

ggplot(distill_methods_dt) + aes(x = partition, y = value, color = method) +
  facet_grid(. ~ name) +
  geom_point() +
  theme_minimal() + theme(
    legend.position = "inside", legend.position.inside = c(1, 0),
    legend.justification.inside = c(1, 0), text = element_text(size = 10)
  ) +
  scale_y_log10("Deaths (log 10 scale)") +
  scale_x_continuous("Age") +
  scale_color_discrete(
    "Method",
    breaks = c(
      "f_mid", "f_mean",
      "mean_f", "wm_f"
    ),
    labels = function(b) { distill_labels[as.character(b)] }
  )

## ----figlex, echo = FALSE, cache = TRUE, out.width = "100%", fig.width = 8, fig.height = 3, dpi = 600----
data("mxF", package = "wpp2019")
data("mxM", package = "wpp2019")
life_expectancy_dt <- as.data.table(mxF)[,
  .(name, age, mxF = `2020-2025`)
][
  as.data.table(mxM), on = c("name", "age"),
  .(name, age, mxF, mxM = `2020-2025`)
][
  name %in% c("Afghanistan", "United Kingdom")
][pop_dt, on = .(name, age), .(name, age, mx = (mxF*popF + mxM*popM)/(popF + popM))]

life_expectancy_dt[,
  ax := 0.5
][,
  qx := fifelse(age == max(age), 1, c(diff(age), 0)*mx / (1 + c(diff(age), 0)*mx * (1 - ax)))
]
life_expectancy_dt[age == 0, lx := 1000]
life_expectancy_dt[, lx := {
  tmp <- lx
  for (i in 2:.N) {
    tmp[i] <- (1 - qx[i - 1]) * tmp[i - 1]
  }
  pmax(tmp, 0)
}, by = name]
life_expectancy_dt[,
  Lx := c(
    tail(lx, -1) + head(ax, -1) * (head(lx, -1) - tail(lx, -1)),
    tail(lx, 1)
  ), by = name
]
life_expectancy_dt[, ex := rev(cumsum(rev(Lx) / 1000)), by = name]

life_ex_interpolated_dt <- life_expectancy_dt[,
  .(age = 0:100, ex = stats::approxfun(age, ex)(0:100))
, by=name]

# because of the selected breakpoints, there are some half-ages, so will
# linearly interpolate between the relevant ages
expander <- distill_methods_dt[, .(
  age = unique(partition[partition != round(partition)])
), by = name][, .(
  lower = floor(age), upper = ceiling(age), tar = age, index = seq_len(.N)
), by = name]

lower_dt <- life_ex_interpolated_dt[expander, on = .(name, age = lower), .(name, index, age, ex)]
upper_dt <- life_ex_interpolated_dt[expander, on = .(name, age = upper), .(name, index, age, ex)]

expand_dt <- lower_dt[
  upper_dt, on = .(name, index)
][expander, on = .(name, index)][, .(
  age = tar,
  ex = approxfun(c(age, i.age), c(ex, i.ex))(tar)
), by = .(name, index)
][, .(name, age, ex)]

lex_dt <- rbind(
  life_ex_interpolated_dt[, .(name, age, ex)],
  expand_dt
)

yll_dt <- distill_methods_dt[
  lex_dt, on = .(name, partition = age)
][, .(YLL = sum(value * ex)), by = .(name, method)]

yllord <- yll_dt[name == "Afghanistan", method[order(YLL)]]

ggplot(yll_dt) + aes(x = name, y = YLL / 1e4, fill = method) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_minimal() +
  theme(
    legend.position = "inside", legend.position.inside = c(0.1, 0.9),
    legend.justification.inside = c(0, 1), text = element_text(size = 8)
  ) +
  scale_fill_brewer(
    "Method", palette = "Set1",
    labels = function(b) { distill_labels[as.character(b)] }
  ) + labs(
    x = "Country",
    y = "YLLs per one million infections\n(10,000 person-years)"
  )

