library(tidyverse) |> suppressMessages()
library(gt) |> suppressMessages()


######## Figures ########
dir.create(here::here("outputs", "figures"),
           recursive = TRUE, showWarnings = FALSE)

#### figure 2 ####
# vector con nombres de archivos
files <- list.files(here::here("data/performance"))

# vectores auxiliares de horizonte y frecuencia
file_names_segments <- map_vec(files, ~ str_split(.x, "_"))
fcst_hrz <- NULL
var_ind <- NULL

for(i in 1:length(files)) {
  fcst_hrz <- c(fcst_hrz, file_names_segments[[i]][2])
  var_ind <- c(var_ind, file_names_segments[[i]][1])
}

fcst_hrz <- str_remove(fcst_hrz, ".csv")

# medias de AE
ae_means <- files |>
  map_dbl(~ suppressMessages( read_csv(paste0(here::here("data/performance/",.x)),
                                      show_col_types = FALSE) ) |>
            mutate(AE = abs(IGAE - IOAE)) |>
            pull(AE) |>
            mean(na.rm = TRUE))

# Función que agrega elementos al conjunto de datos de entrada
gen_df <- function(file, fcst_hrz, var_ind) {
  
  df <- suppressMessages( read_csv(here::here(file), show_col_types = FALSE) ) |>
    rename(periodo = 1) |>
    mutate(fcst_hrz = fcst_hrz,
           var_ind = var_ind) 
  
  return(df)
} 

# Lectura de conjunto de datos de entrada y transformaciones requeridas
for(i in 1:length(files)) {
  gen_df(paste0("data/performance/", files[i]), fcst_hrz[i], var_ind[i]) |>
    mutate(lab_x = "2024/03",
           lab_y = if_else(periodo == lab_x, 5, NA_real_),
           txt = paste0("MAE: ", round(ae_means[i], 2)),
           txt = if_else(!is.na(lab_y), txt, NA_character_)) |>
    assign(paste("df", var_ind[i], fcst_hrz[i], sep = "_"), value = _)
}

# Consolidación de conjunto de datos
df <- `df_LEVELS_T+1` |>
  bind_rows(
    `df_LEVELS_T+2`,
    `df_MV_T+1`,
    `df_MV_T+2`,
    `df_AV_T+1`,
    `df_AV_T+2`
  ) 

df_plot <- df |> 
  mutate(fcst_hrz = factor(fcst_hrz, levels = c("T+1", "T+2")),
         var_ind = factor(var_ind, levels = c("LEVELS", "MV", "AV")))


# Graficación de datos
df |> mutate(fcst_hrz = factor(fcst_hrz, levels = c("T+1", "T+2")),
         var_ind = factor(var_ind, levels = c("LEVELS", "MV", "AV"))) |>
  ggplot(aes(periodo)) +
  geom_errorbar(aes(ymin = Inferior, ymax = Superior),
                linewidth = .5) +
  geom_point(aes(y = IOAE, col = "IOAE"), size = 1, alpha = 0.75) +
  geom_point(aes(y = IGAE, col = "IGAE"), size = 1, alpha = 0.75) +
  geom_text(aes(x = "2024/03", y = Inf, label = txt, vjust = 1.75),
            size = 4) +
  scale_color_manual(values = c("IOAE" = "black",
                                "IGAE" = "red")) +
  facet_grid(rows = vars(var_ind), cols = vars(fcst_hrz), scales = "free_y") +
  theme_bw() +
  guides(
    color = guide_legend(
      title = NULL,
      size = 3
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(paste0("%", 
              glue::glue_collapse(rep(" ", 25)),
              "%",
              glue::glue_collapse(rep(" ", 25)),
              "Index")) -> fig_2

ggsave(fig_2,
       filename = "figure_2.png",
       path = here::here("outputs/figures/"),
       device = "png",
       width = 1046,
       height = 706,
       units = "px",
       dpi = 96)


# Graficación de datos BW
df |> mutate(fcst_hrz = factor(fcst_hrz, levels = c("T+1", "T+2")),
             var_ind = factor(var_ind, levels = c("LEVELS", "MV", "AV"))) |>
  ggplot(aes(periodo)) +
  geom_errorbar(aes(ymin = Inferior, ymax = Superior),
                linewidth = .5) +
  geom_point(aes(y = IOAE, shape = "IOAE", fill = "IOAE"), size = 1, alpha = 0.75) +
  geom_point(aes(y = IGAE, shape = "IGAE", fill = "IGAE"), size = 1, alpha = 0.75) +
  geom_text(aes(x = "2024/03", y = Inf, label = txt, vjust = 1.75),
            size = 4) +
  scale_shape_manual(values = c("IOAE" = 21, "IGAE" = 22),
                     breaks = c("IOAE", "IGAE"),
                     name = NULL) +
  scale_fill_manual(values = c("IOAE" = "black", "IGAE" = "grey90"),
                    breaks = c("IOAE", "IAGE"),
                    name = NULL) +
  facet_grid(rows = vars(var_ind), cols = vars(fcst_hrz), scales = "free_y") +
  guides(
    shape = guide_legend(
      override.aes = list(
        shape = c(21, 22),
        fill  = c("black", "grey90"),
        color = "black",
        size  = 3
      )
    ),
    fill = "none"  
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(paste0("%", 
              glue::glue_collapse(rep(" ", 25)),
              "%",
              glue::glue_collapse(rep(" ", 25)),
              "Index")) -> fig_2_bw

ggsave(fig_2_bw,
       filename = "figure_2_bw.png",
       path = here::here("outputs/figures/"),
       device = "png",
       width = 1046,
       height = 706,
       units = "px",
       dpi = 96)




#### figure 3 ####

# vector con nombres de archivos
files <- list.files(here::here("data/ioae_pibo_gdp"))
files <- files[1:6]

# vectores auxiliares de horizonte y frecuencia
file_names_segments <- map_vec(files, ~ str_split(.x, "_"))
fcst_hrz <- NULL
var_freq <- NULL

for(i in 1:length(files)) {
  fcst_hrz <- c(fcst_hrz, file_names_segments[[i]][3])
  var_freq <- c(var_freq, file_names_segments[[i]][2])
}

fcst_hrz <- str_remove(fcst_hrz, ".csv")

# medias de AE
files_ <- paste(here::here("data/ioae_pibo_gdp/"), files, sep = "/")
ae_means <- files_ |>
  map_dbl(~ suppressMessages( read_csv(.x, show_col_types = FALSE) ) |>
            pull(AE) |>
            mean(na.rm = TRUE))

# Función que agrega elementos al conjunto de datos de entrada
gen_df <- function(file, fcst_hrz, var_freq) {
  
  df <- suppressMessages( read_csv(file, show_col_types = FALSE) ) |>
    rename(periodo = 1) |>
    mutate(mean_indicator = if_else(fcst_hrz == "T+0","IGAE","IOAE"),
           fcst_hrz = fcst_hrz,
           var_freq = if_else(var_freq == "Trimestral","Quarterly", "Annual")) 
  
  return(df)
} 

# Lectura de conjunto de datos de entrada y transformaciones requeridas
for(i in 1:length(files_)) {
  gen_df(files_[i], fcst_hrz[i], var_freq[i]) |>
    mutate(lab_x = "2024/03",
           lab_y = if_else(periodo == lab_x, 5, NA_real_),
           txt = paste0("MAE: ", round(ae_means[i], 2)),
           txt = if_else(!is.na(lab_y), txt, NA_character_)) |>
    assign(paste("df", var_freq[i], fcst_hrz[i], sep = "_"), value = _)
}

# Consolidación de conjunto de datos
df <- `df_Trimestral_T+0` |>
  bind_rows(
    `df_Trimestral_T+1`,
    `df_Trimestral_T+2`,
    `df_Anual_T+0`,
    `df_Anual_T+1`,
    `df_Anual_T+2`
  ) |> #-> df_
  mutate(across(2:6, ~ if_else(periodo == "2021/06" & var_freq == "Annual", NA_real_, .x)))


# Graficación de datos
df |> 
  filter(periodo != "2020/09") |>
  mutate(fcst_hrz = if_else(fcst_hrz == "T+0", "T", fcst_hrz),
         lower = if_else(fcst_hrz == "T", NA_real_, lower),
         upper = if_else(fcst_hrz == "T", NA_real_, upper),
         var_freq = factor(var_freq, levels = c("Quarterly", "Annual"))) |>
  ggplot(aes(periodo)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = mean, col = "IOAE"), size = 1, alpha = 0.75) +
  geom_point(aes(y = PIBO, col = "PIBO"), size = 1, alpha = 0.75) +
  geom_point(aes(y = PIB, col = "GDP"), size = 1, alpha = 0.75) +
  geom_text(aes(x = lab_x, y = lab_y, label = txt),
            size = 4, vjust = "left") +
  scale_color_manual(values = c("IOAE" = "black",
                                "PIBO" = "blue",
                                "GDP" = "red")) +
  facet_grid(rows = vars(fcst_hrz), cols = vars(var_freq)) +
  theme_bw() +
  guides(
    color = guide_legend(
      title = NULL,
      size = 3
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(label = "%") -> fig_3

ggsave(fig_3,
       filename = "figure_3.png",
       path = here::here("outputs/figures/"),
       device = "png",
       width = 791,
       height = 622,
       units = "px",
       dpi = 96)



# Graficación de datos BW
df |>
  filter(periodo != "2020/09") |>
  mutate(fcst_hrz = if_else(fcst_hrz == "T+0", "T", fcst_hrz),
         lower = if_else(fcst_hrz == "T", NA_real_, lower),
         upper = if_else(fcst_hrz == "T", NA_real_, upper),
         var_freq = factor(var_freq, levels = c("Quarterly", "Annual"))) |>
  ggplot(aes(periodo)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = mean, shape = "IOAE", fill = "IOAE"), size = 1.5, alpha = 0.75) +
  geom_point(aes(y = PIBO, shape = "PIBO", fill = "PIBO"), size = 1.5, alpha = 0.75) +
  geom_point(aes(y = PIB, shape = "GDP", fill = "GDP"), size = 1.5, alpha = 0.75) +
  geom_text(aes(x = lab_x, y = lab_y, label = txt),
            size = 4, vjust = "left") +
  scale_shape_manual(values = c("IOAE" = 21, "PIBO" = 22, "GDP" = 23),
                     breaks = c("IOAE", "PIBO", "GDP"),
                     name = NULL) +
  scale_fill_manual(values = c("IOAE" = "black", "PIBO" = "grey50", "GDP" = "grey90"),
                    breaks = c("IOAE", "PIBO", "GDP"),
                    name = NULL) +
  facet_grid(rows = vars(fcst_hrz), cols = vars(var_freq)) +
  guides(
    shape = guide_legend(
      override.aes = list(
        shape = c(21, 22, 23),
        fill  = c("black", "grey50", "grey90"),
        color = "black",
        size  = 3
      )
    ),
    fill = "none"  # evita una segunda leyenda separada
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(label = "%") -> fig_3_bw

ggsave(fig_3_bw,
       filename = "figure_3_bw.png",
       path = here::here("outputs/figures/"),
       device = "png",
       width = 791,
       height = 622,
       units = "px",
       dpi = 96)


#### figure 5 ####
df <- readRDS(here::here("outputs/rds/figures_inputs/f5.rds"))

df <- df |>
  mutate(variable = if_else(variable == "Conf Comer", "Conf Trade", variable))

TS <- readRDS(here::here("outputs/rds/figures_inputs/TS.rds"))

df |>
  mutate(variable = factor(variable,
                           levels = TS |> 
                             mutate(
                               Names = if_else(Names == "Conf Comer", "Conf Trade", Names)
                             )|> pull(Names))) |>
  ggplot(aes(x = variable)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = Phat), size = 1, alpha = 0.75) +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  labs(
    y = "Load",
    x = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1),
    panel.grid = element_line(color = "grey95")
  ) -> fig_5

ggsave(fig_5,
       filename = "figure_5.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 659,
       height = 540,
       units = "px",
       dpi = 96)


#  BW version
df |>
  mutate(variable = factor(variable,
                           levels = TS |> 
                             mutate(
                               Names = if_else(Names == "Conf Comer", "Conf Trade", Names)
                             )|> pull(Names))) |>
  ggplot(aes(x = variable)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = Phat), size = 1, alpha = 0.75) +
  geom_hline(yintercept = 0, color = "grey40", linetype = 2) +
  labs(
    y = "Load",
    x = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1),
    panel.grid = element_line(color = "grey95")
  ) -> fig_5_bw

ggsave(fig_5_bw,
       filename = "figure_5_bw.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 659,
       height = 540,
       units = "px",
       dpi = 96)

#### figure 6 ####
df <- readRDS(here::here("outputs/rds/figures_inputs/f6.rds"))

df <- df |>
  pivot_wider(names_from = name, values_from = value)


df |>
  ggplot(aes(period)) +
  geom_line(aes(y = igae, color = "igae"), linewidth = 1) +
  geom_line(aes(y = fhat, color = "fhat")) +
  scale_x_date(date_labels = "%Y") +
  scale_color_manual(labels = c(expression(hat(F)), "IGAE"),
    values = c("igae" = "black", "fhat" = "red")) +
  labs(
    col = ""
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_line(color = "grey95"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) -> fig_6

ggsave(fig_6,
       filename = "figure_6.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 659,
       height = 540,
       units = "px",
       dpi = 96)


#BW version
df |>
  ggplot(aes(period)) +
  geom_line(aes(y = igae, linetype = "igae", color = "IGAE"), linewidth = 1) +
  geom_line(aes(y = fhat, linetype = "fhat", color = "fhat")) +
  scale_x_date(date_labels = "%Y") +
  scale_linetype_manual(labels = c(expression(hat(F)), "IGAE"),
                        values = c("fhat" = "solid", "igae" = "solid")) +
  scale_color_manual(labels = c(expression(hat(F)), "IGAE"),
                     values = c("fhat" = "grey30", "IGAE" = "black")) +
  labs(
    col = "", linetype = ""
  ) +
  theme_bw() +
  guides(
    linetype = guide_legend(
      override.aes = list(
        linetype = c("solid", "solid"),
        color  = c("black", "grey30")
      )
    ),
    color = "none"
  ) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_line(color = "grey95"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) -> fig_6_bw

ggsave(fig_6_bw,
       filename = "figure_6_bw.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 659,
       height = 540,
       units = "px",
       dpi = 96)


#### figure 7 ####
df <- readRDS(here::here("outputs/rds/figures_inputs/f7.rds"))

maes <- df |>
  group_by(type, horizon) |>
  summarise(mae = mean(abs(mean - observed)), .groups = "drop_last") |>
  ungroup() |>
  mutate(mae = round(mae, 2)) 

df |>
  mutate(id = paste(type, horizon, sep = "_")) |>
  mutate(type = case_when(
    type == "annual" ~ "AV",
    type == "monthly" ~ "MV",
    TRUE ~ "Levels"
  ),
  type = factor(type, levels = c("Levels", "MV", "AV"))) |>
  left_join(maes |>
              mutate(id = paste(type, horizon, sep = "_")) |>
              select(id, mae),
            by = join_by(id)) |>
  mutate(horizon = if_else(horizon == 1, "T+1", "T+2"),
         horizon = factor(horizon, levels = c("T+1", "T+2")),
         mae_txt = if_else(period == "2023/11",
                           paste0("MAE: ",mae), NA_character_),
         lab_x = "2023/11") |>
  ggplot(aes(period)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = mean, col = "IOAE"), size = 0.8, alpha = 0.75) +
  geom_point(aes(y = observed, col = "IGAE"), size = 0.8, alpha = 0.75) +
  geom_text(aes(x = lab_x, y = -Inf, label = mae_txt),
            size = 4, vjust = -1.5) +
  scale_color_manual(values = c("IOAE" = "black",
                                "IGAE" = "red")) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) +
  facet_grid(rows = vars(type), cols = vars(horizon),
             scales = "free_y") +
  theme_bw() +
  guides(
    color = guide_legend(
      color = c("black", "red")
    )
  ) +
  labs(color = "") +
  theme(
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(label = paste0("%",
                      glue::glue_collapse(rep(" ", 40)),
                      "%",
                      glue::glue_collapse(rep(" ", 40)),
                      "Index")) -> fig_7

ggsave(fig_7,
       filename = "figure_7.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 866,
       height = 640,
       units = "px",
       dpi = 90)


# BW version
df |>
  mutate(id = paste(type, horizon, sep = "_")) |>
  mutate(type = case_when(
    type == "annual" ~ "AV",
    type == "monthly" ~ "MV",
    TRUE ~ "Levels"
  ),
  type = factor(type, levels = c("Levels", "MV", "AV"))) |>
  left_join(maes |>
              mutate(id = paste(type, horizon, sep = "_")) |>
              select(id, mae),
            by = join_by(id)) |>
  mutate(horizon = if_else(horizon == 1, "T+1", "T+2"),
         horizon = factor(horizon, levels = c("T+1", "T+2")),
         mae_txt = if_else(period == "2023/11",
                           paste0("MAE: ",mae), NA_character_),
         lab_x = "2023/11") |>
  ggplot(aes(period)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                linewidth = .5) +
  geom_point(aes(y = mean, shape = "IOAE", fill = "IOAE"), size = 0.8, alpha = 0.75) +
  geom_point(aes(y = observed, shape = "IGAE", fill = "IGAE"), size = 0.8, alpha = 0.75) +
  geom_text(aes(x = lab_x, y = -Inf, label = mae_txt),
            size = 4, vjust = -1.5) +
  scale_shape_manual(values = c("IOAE" = 21, "IGAE" = 23),
                     breaks = c("IOAE", "IGAE"),
                     name = NULL) +
  scale_fill_manual(values = c("IOAE" = "black", "IGAE" = "grey90"),
                    breaks = c("IOAE", "IGAE"),
                    name = NULL) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 2)]) +
  facet_grid(rows = vars(type), cols = vars(horizon),
             scales = "free_y") +
  theme_bw() +
  guides(
    shape = guide_legend(
      override.aes = list(
        shape = c(21, 23),
        fill  = c("black", "grey90"),
        color = "black",
        size  = 3
      )
    ),
    fill = "none"  
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid = element_line(color = "grey95")
  ) +
  ylab(label = paste0("%",
                      glue::glue_collapse(rep(" ", 40)),
                      "%",
                      glue::glue_collapse(rep(" ", 40)),
                      "Index")) -> fig_7_bw

ggsave(fig_7_bw,
       filename = "figure_7_bw.png",
       path = here::here("outputs/figures"),
       device = "png",
       width = 866,
       height = 640,
       units = "px",
       dpi = 90)


#### figures appendix ####
## cumulative MAE
archivos <- list.files(here::here("outputs/rds/figures_inputs"))
archivos <- archivos[str_detect(archivos, "mat_mae")]
archivos_nom <- str_extract(archivos, "(?<=mat_mae_).+(?=\\.rds)")

l_arch_nom <- list(
  arch_nom_1 <- archivos[str_detect(archivos_nom, "(?<=_)1")],
  arch_nom_2 <- archivos[str_detect(archivos_nom, "(?<=_)2")]
)


modelos_nom <- str_extract(archivos_nom, ".+(?=_)") |> unique()


# Preparing data with tibble format
f_df_tm <- function(mts, ind_mod, ind_id) {
  
  periodo <- mts |>
    time() |>
    zoo::as.Date.yearmon()
  
  modelo_nom <- modelos_nom[ind_mod]
  modelo_id <- str_extract(v_arch_target[ind_id], "\\d{2}(?=\\.)")
  
  df <- mts |>
    as_tibble() |>
    mutate(
      periodo = periodo
    ) |>
    pivot_longer(1:2) |>
    mutate(
      modelo = if_else(name == "Series 1", "dfm", modelo_nom)
    ) |>
    select(-name) |>
    pivot_wider(
      names_from = modelo,
      values_from = value
    ) |>
    mutate(
      id = str_sub(modelo_id, 1, 1) |> as.integer(),
      h =str_sub(modelo_id, 2, 2) |> as.integer()
    )
  
  return(df)
}


l_df_mod <- list()
l_df_tipo <- list()


for(k in 1:2) {
  
  # Generatin models list according the type of series (full vs selected)
  for(j in 1:length(modelos_nom)) {
    
    v_arch_target <- l_arch_nom[[k]][str_detect(l_arch_nom[[k]], modelos_nom[j])]
    
    # Integrating by horizon
    for(i in 1:length(v_arch_target)) {
      if(i == 1) {
        df <- f_df_tm(readRDS(here::here(paste0("outputs/rds/figures_inputs/", v_arch_target[1]))), j, i)
      } else {
        df <- df |> bind_rows(
          f_df_tm(readRDS(here::here(paste0("outputs/rds/figures_inputs/", v_arch_target[i]))), j, i)
        )
      }
    }
    
    l_df_mod[[j]] <- df
    
  }
  
  l_df_tipo[[k]] <- l_df_mod
  
}

for(i in 1:4) {
  if(i == 1) {
    df_t1 <- l_df_tipo[[1]][[i]]
  } else {
    df_t1 <- df_t1 |>
      bind_cols(
        l_df_tipo[[1]][[i]] |>
          select(-c(periodo, dfm, id, h))
      )
  }
}

for(i in 1:4) {
  if(i == 1) {
    df_t2 <- l_df_tipo[[2]][[i]]
  } else {
    df_t2 <- df_t2 |>
      bind_cols(
        l_df_tipo[[2]][[i]] |>
          select(-c(periodo, dfm, id, h))
      )
  }
}

df <- bind_rows(df_t1, df_t2)

# Plotting both kinds of series: full and selected
for(i in 1:2) {
  # selecting series
  filtro <- i # 1 o 2
  # plot
  df |>
    pivot_longer(
      -c(periodo, id, h, dfm)
    ) |>
    filter(id == filtro) |>
    mutate(
      h = if_else(h == 1, "T+1", "T+2"),
      name = str_to_upper(name),
      name = factor(name, levels = c("LASSO", "FAVAR", "MLP", "ARIMA"))
    ) |>
    ggplot(aes(periodo, group = name)) +
    geom_rect(aes(xmin = as.Date("2020/03/01"),
                  xmax = as.Date("2020/06/01"),
                  ymin = -Inf,
                  ymax = Inf),
              fill = "grey90",
              alpha = 0.5) +
    geom_line(aes(y = dfm, color = "dfm"), linewidth = 1, alpha = 0.8) +
    geom_line(aes(y = value, color = "value"), linewidth = 1, alpha = 0.8) +
    scale_x_date(date_labels = "%Y", breaks = as.Date(c("2019/01/01",
                                                        "2020/01/01",
                                                        "2021/01/01",
                                                        "2022/01/01",
                                                        "2023/01/01",
                                                        "2024/01/01"))) +
    scale_color_manual(
      values = c("dfm" = "black", "value" = "red"),
      labels = c("DFM", "Competitor")
    ) +
    facet_grid(
      rows = vars(name),
      cols = vars(h)
    ) +
    labs(x = "", y = "Cumulative MAE", color = "") +
    guides(
      color = guide_legend(
        override.aes = list(
          color = c("black","red"),
          labels = c("DFM", "Competitor")
        )
      )
    ) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey95"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) 
  
  ggsave(filename = paste0("figure_app_", filtro, ".png"),
         path = here::here("outputs/figures"),
         device = "png",
         width = 890,
         height = 940,
         units = "px",
         dpi = 90)
  
  
  # BW plot
  df |>
    pivot_longer(
      -c(periodo, id, h, dfm)
    ) |>
    filter(id == filtro) |>
    mutate(
      h = if_else(h == 1, "T+1", "T+2"),
      name = str_to_upper(name),
      name = factor(name, levels = c("LASSO", "FAVAR", "MLP", "ARIMA"))
    ) |>
    ggplot(aes(periodo, group = name)) +
    geom_rect(aes(xmin = as.Date("2020/03/01"),
                  xmax = as.Date("2020/06/01"),
                  ymin = -Inf,
                  ymax = Inf),
              fill = "grey90",
              alpha = 0.5) +
    geom_line(aes(y = dfm, linetype = "dfm", color = "dfm"), linewidth = 1, alpha = 0.8) +
    geom_line(aes(y = value, linetype = "value", color = "value"), alpha = 0.8) +
    scale_x_date(date_labels = "%Y", breaks = as.Date(c("2019/01/01",
                                                        "2020/01/01",
                                                        "2021/01/01",
                                                        "2022/01/01",
                                                        "2023/01/01",
                                                        "2024/01/01"))) +
    scale_color_manual(
      values = c("dfm" = "black", "value" = "grey30"),
      labels = c("DFM", "Competitor")
    ) +
    scale_linetype_manual(
      values = c("dfm" = "solid", "value" = "solid")
    ) +
    facet_grid(
      rows = vars(name),
      cols = vars(h)
    ) +
    labs(x = "", y = "Cumulative MAE", color = "", linetype = "") +
    guides(
      color = guide_legend(
        override.aes = list(
          color = c("black","grey30"),
          linetype = c("solid", "solid")
        )
      ),
      linetype = "none"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey95"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) 
  
  ggsave(filename = paste0("figure_app_", filtro, "_bw.png"),
         path = here::here("outputs/figures"),
         device = "png",
         width = 890,
         height = 940,
         units = "px",
         dpi = 90)
  
}


## cumulative RMSE
archivos <- list.files(here::here("outputs/rds/figures_inputs"))
archivos <- archivos[str_detect(archivos, "mat_rmse")]
archivos_nom <- str_extract(archivos, "(?<=mat_rmse_).+(?=\\.rds)")

l_arch_nom <- list(
  arch_nom_1 <- archivos[str_detect(archivos_nom, "(?<=_)1")],
  arch_nom_2 <- archivos[str_detect(archivos_nom, "(?<=_)2")]
)

modelos_nom <- str_extract(archivos_nom, ".+(?=_)") |> unique()

# Preparing data with tibble format
f_df_tm <- function(mts, ind_mod, ind_id) {
  
  periodo <- mts |>
    time() |>
    zoo::as.Date.yearmon()
  
  modelo_nom <- modelos_nom[ind_mod]
  modelo_id <- str_extract(v_arch_target[ind_id], "\\d{2}(?=\\.)")
  
  df <- mts |>
    as_tibble() |>
    mutate(
      periodo = periodo
    ) |>
    pivot_longer(1:2) |>
    mutate(
      modelo = if_else(name == "Series 1", "dfm", modelo_nom)
    ) |>
    select(-name) |>
    pivot_wider(
      names_from = modelo,
      values_from = value
    ) |>
    mutate(
      id = str_sub(modelo_id, 1, 1) |> as.integer(),
      h =str_sub(modelo_id, 2, 2) |> as.integer()
    )
  
  return(df)
}


l_df_mod <- list()
l_df_tipo <- list()


for(k in 1:2) {
  
  # Generatin models list according the type of series (full vs selected)
  for(j in 1:length(modelos_nom)) {
    
    v_arch_target <- l_arch_nom[[k]][str_detect(l_arch_nom[[k]], modelos_nom[j])]
    
    # Integrating by horizon
    for(i in 1:length(v_arch_target)) {
      if(i == 1) {
        df <- f_df_tm(readRDS(here::here(paste0("outputs/rds/figures_inputs/", v_arch_target[1]))), j, i)
      } else {
        df <- df |> bind_rows(
          f_df_tm(readRDS(here::here(paste0("outputs/rds/figures_inputs/", v_arch_target[i]))), j, i)
        )
      }
    }
    
    l_df_mod[[j]] <- df
    
  }
  
  l_df_tipo[[k]] <- l_df_mod
  
}

for(i in 1:4) {
  if(i == 1) {
    df_t1 <- l_df_tipo[[1]][[i]]
  } else {
    df_t1 <- df_t1 |>
      bind_cols(
        l_df_tipo[[1]][[i]] |>
          select(-c(periodo, dfm, id, h))
      )
  }
}

for(i in 1:4) {
  if(i == 1) {
    df_t2 <- l_df_tipo[[2]][[i]]
  } else {
    df_t2 <- df_t2 |>
      bind_cols(
        l_df_tipo[[2]][[i]] |>
          select(-c(periodo, dfm, id, h))
      )
  }
}

df <- bind_rows(df_t1, df_t2)

# Plotting both kinds of series: full and selected
for(i in 1:2) {
  # selecting series
  filtro <- i # 1 o 2
  # plot
  df |>
    pivot_longer(
      -c(periodo, id, h, dfm)
    ) |>
    filter(id == filtro) |>
    mutate(
      h = if_else(h == 1, "T+1", "T+2"),
      name = str_to_upper(name),
      name = factor(name, levels = c("LASSO", "FAVAR", "MLP", "ARIMA"))
    ) |>
    ggplot(aes(periodo, group = name)) +
    geom_rect(aes(xmin = as.Date("2020/03/01"),
                  xmax = as.Date("2020/06/01"),
                  ymin = -Inf,
                  ymax = Inf),
              fill = "grey90",
              alpha = 0.5) +
    geom_line(aes(y = dfm, color = "dfm"), linewidth = 1, alpha = 0.8) +
    geom_line(aes(y = value, color = "value"), linewidth = 1, alpha = 0.8) +
    scale_x_date(date_labels = "%Y", breaks = as.Date(c("2019/01/01",
                                                        "2020/01/01",
                                                        "2021/01/01",
                                                        "2022/01/01",
                                                        "2023/01/01",
                                                        "2024/01/01"))) +
    scale_color_manual(
      values = c("dfm" = "black", "value" = "red"),
      labels = c("DFM", "Competitor")
    ) +
    facet_grid(
      rows = vars(name),
      cols = vars(h)
    ) +
    labs(x = "", y = "Cumulative RMSE", color = "") +
    guides(
      color = guide_legend(
        override.aes = list(
          color = c("black","red"),
          labels = c("DFM", "Competitor")
        )
      )
    ) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey95"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) 
  
  ggsave(filename = paste0("figure_app_", filtro + 2, ".png"),
         path = here::here("outputs/figures"),
         device = "png",
         width = 890,
         height = 940,
         units = "px",
         dpi = 90)
  
  
  # BW plot
  df |>
    pivot_longer(
      -c(periodo, id, h, dfm)
    ) |>
    filter(id == filtro) |>
    mutate(
      h = if_else(h == 1, "T+1", "T+2"),
      name = str_to_upper(name),
      name = factor(name, levels = c("LASSO", "FAVAR", "MLP", "ARIMA"))
    ) |>
    ggplot(aes(periodo, group = name)) +
    geom_rect(aes(xmin = as.Date("2020/03/01"),
                  xmax = as.Date("2020/06/01"),
                  ymin = -Inf,
                  ymax = Inf),
              fill = "grey90",
              alpha = 0.5) +
    geom_line(aes(y = dfm, linetype = "dfm", color = "dfm"), linewidth = 1, alpha = 0.8) +
    geom_line(aes(y = value, linetype = "value", color = "value"), alpha = 0.8) +
    scale_x_date(date_labels = "%Y", breaks = as.Date(c("2019/01/01",
                                                        "2020/01/01",
                                                        "2021/01/01",
                                                        "2022/01/01",
                                                        "2023/01/01",
                                                        "2024/01/01"))) +
    scale_color_manual(
      values = c("dfm" = "black", "value" = "grey30"),
      labels = c("DFM", "Competitor")
    ) +
    scale_linetype_manual(
      values = c("dfm" = "solid", "value" = "solid")
    ) +
    facet_grid(
      rows = vars(name),
      cols = vars(h)
    ) +
    labs(x = "", y = "Cumulative RMSE", color = "", linetype = "") +
    guides(
      color = guide_legend(
        override.aes = list(
          color = c("black","grey30"),
          linetype = c("solid", "solid")
        )
      ),
      linetype = "none"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_line(color = "grey95"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ) 
  
  ggsave(filename = paste0("figure_app_", filtro + 2, "_bw.png"),
         path = here::here("outputs/figures"),
         device = "png",
         width = 890,
         height = 940,
         units = "px",
         dpi = 90)
  
}

######## Tables ########
dir.create(here::here("outputs", "rds", "prepared_tables"),
           recursive = TRUE, showWarnings = FALSE)

#### table 1 ####
readRDS(here::here("outputs/rds/models_objects/Trans.rds")) |>
  gt() |>
  fmt_number(
    columns = 3,
    decimals = 4
  ) |>
  cols_width(
    1 ~ px(160),
    everything() ~ px(110)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t1.rds"))


#### table n_factors ####
readRDS(here::here("outputs/rds/models_objects/df_n_factors.rds")) |>
  rename(Period = Date) |>
  mutate(Period = format(Period, "%Y-%m")) |>
  gt() |>
  cols_label(
    Bai_Ng_IC1 = html("IC<sub>1</sub>"),
    Bai_Ng_IC2 = html("IC<sub>2</sub>"),
    Bai_Ng_IC3 = html("IC<sub>3</sub>"),
    Onatski_ed = "ed",
    Ahn_Horenstein_ker = "ker",
    Ahn_Horenstein_ker = "kgr"
  ) |>
  tab_spanner(
    label = "Bai & Ng",
    columns = 2:4
  ) |>
  tab_spanner(
    label = "Onatski",
    columns = 5
  ) |>
  tab_spanner(
    label = "Ahn & Horenstein",
    columns = 6:7
  ) |>
  cols_width(
    everything() ~ px(80)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t_n_factors.rds"))





#### table common factor ADF test ####
readRDS(here::here("outputs/rds/models_objects/df_fct_adf_test_pv.rds")) |>
  gt() |>
  cols_label(
    DF_statistic = "DF statistic",
    p_value = "p-value"
  ) |>
  fmt_number(
    decimals = 3
  ) |>
  cols_width(
    1:2 ~ px(110)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/fct_adf_test.rds"))


#### table pooled test ####
readRDS(here::here("outputs/rds/models_objects/df_pooled_tests_pv.rds")) |>
  gt() |>
  cols_label(
    P_statistic = "P statistic",
    p_value = "p-value"
  ) |>
  fmt_number(
    decimals = 3
  ) |>
  cols_width(
    1:2 ~ px(110)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/pooled_test.rds"))



#### table 2 ####
models <- c("LASSO", "FAVAR", "MLP", "ARIMA")
readRDS(here::here("outputs/rds/models_objects/dm_mae.rds")) |>
  mutate(
    groups = if_else(`Variable / Model` %in% models,
                     "Econometric and ML models",
                     "Individual indactors")
  ) |>
  relocate(groups) |>
  gt(groupname_col = "groups") |>
  cols_label(
    `DM h = 1` = html("DM<sub>HAC</sub>"),
    `p-value h = 1` = "p-value",
    `DM h = 2` = html("DM<sub>HAC</sub>"),
    `p-value h = 2` = "p-value"
  ) |>
  tab_spanner(
    label = "h = 1",
    columns = 3:4
  ) |>
  tab_spanner(
    label = "h = 2",
    columns = 5:6
  ) |>
  cols_width(
    2 ~ px(160),
    everything() ~ px(120)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t2.rds"))


#### SPA tests ####
readRDS(here::here("outputs/rds/models_objects/spa_test.rds")) |>
  pivot_wider(
    names_from = h,
    values_from = c(Statistic, p_value)
  ) |>
  relocate(
    p_value_1,
    .after = Statistic_1
  ) |>
  gt() |>
  tab_spanner(
    label = "h = 1",
    columns = 2:3
  ) |>
  tab_spanner(
    label = "h = 2",
    columns = 4:5
  ) |>
  fmt_number(
    columns = 2:5,
    decimals = 4
  ) |>
  cols_label(
    Metric ~ "",
    starts_with("S") ~ "Statistic",
    starts_with("p") ~ "p-value"
  ) |>
  cols_width(
    everything() ~ px(110)
  ) |>
  saveRDS("outputs/rds/prepared_tables/t_spa_test.rds")


#### table 3 ####
v_pvalue <- readRDS(here::here("outputs/rds/models_objects/mcs_mae_pvalue.rds")) |>
  mutate(p_value = as.character(p_value)) |>
  pull(p_value)
v_pvalue <- c("p-value", v_pvalue[1], NA, v_pvalue[2], NA)

df_aux <- readRDS(here::here("outputs/rds/models_objects/mcs_mae.rds"))

lbl <- names(df_aux)

df_aux |>
  mutate(
    `Rank h = 1` = if_else(is.na(`TR h = 1`),
                           NA_character_,
                           paste0("(", `Rank h = 1`, ")")),
    `Rank h = 2` = if_else(is.na(`TR h = 2`),
                           NA_character_,
                           paste0("(", `Rank h = 2`, ")"))
  ) |>
  bind_rows(
    {
      tibble(
        pivot = rep(1, length(v_pvalue)),
        lbl = lbl,
        p_value = v_pvalue
      ) |>
        pivot_wider(
          names_from = lbl,
          values_from = p_value
        ) |>
        select(-pivot) |>
        mutate(
          across(c(3,5), as.numeric)
        )
    }
  ) |>
  relocate(1, 3, 2, 5, 4) |>
  mutate(
    groups = case_when(
      `Variable / Model` == "p-value" ~ "Global",
      `Variable / Model` %in% models ~ "Econometric and ML models",
      TRUE ~ "Individual indicators"
    )
  ) |>
  mutate(
    `Variable / Model` = if_else(`Variable / Model` == "DFM", "IOAE", `Variable / Model`)
  ) |>
  relocate(groups) |>
  gt(groupname_col = "groups") |>
  tab_spanner(
    label = "h = 1",
    columns = 3:4
  ) |>
  tab_spanner(
    label = "h = 2",
    columns = 5:6
  ) |>
  cols_label(
    `TR h = 1` = html("T<sub>R,i</sub>"),
    `Rank h = 1` = html("(Rank<sub>R</sub>)"),
    `TR h = 2` = html("T<sub>R,i</sub>"),
    `Rank h = 2` = html("(Rank<sub>R</sub>)")
  ) |>
  sub_missing(
    columns = c(3, 5),
    rows = `Variable / Model` == "p-value",
    missing_text = ""
  ) |>
  cols_width(
    2 ~ px(160),
    everything() ~ px(120)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t3.rds"))


#### DM-HAC with RMSE ####
readRDS(here::here("outputs/rds/models_objects/dm_rmse.rds")) |>
  mutate(
    groups = if_else(`Variable / Model` %in% models,
                     "Econometric and ML models",
                     "Individual indactors")
  ) |>
  relocate(groups) |>
  gt(groupname_col = "groups") |>
  cols_label(
    `DM h = 1` = html("DM<sub>HAC</sub>"),
    `p-value h = 1` = "p-value",
    `DM h = 2` = html("DM<sub>HAC</sub>"),
    `p-value h = 2` = "p-value"
  ) |>
  tab_spanner(
    label = "h = 1",
    columns = 3:4
  ) |>
  tab_spanner(
    label = "h = 2",
    columns = 5:6
  ) |>
  cols_width(
    2 ~ px(160),
    everything() ~ px(120)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t_dm_rmse.rds"))


#### MCS with RMSE ####
v_pvalue <- readRDS(here::here("outputs/rds/models_objects/mcs_rmse_pvalue.rds")) |>
  mutate(p_value = as.character(p_value)) |>
  pull(p_value)
v_pvalue <- c("p-value", v_pvalue[1], NA, v_pvalue[2], NA)

df_aux <- readRDS(here::here("outputs/rds/models_objects/mcs_rmse.rds"))

lbl <- names(df_aux)

df_aux |>
  mutate(
    `Rank h = 1` = if_else(is.na(`TR h = 1`),
                           NA_character_,
                           paste0("(", `Rank h = 1`, ")")),
    `Rank h = 2` = if_else(is.na(`TR h = 2`),
                           NA_character_,
                           paste0("(", `Rank h = 2`, ")"))
  ) |>
  bind_rows(
    {
      tibble(
        pivot = rep(1, length(v_pvalue)),
        lbl = lbl,
        p_value = v_pvalue
      ) |>
        pivot_wider(
          names_from = lbl,
          values_from = p_value
        ) |>
        select(-pivot) |>
        mutate(
          across(c(3,5), as.numeric)
        )
    }
  ) |>
  relocate(1, 3, 2, 5, 4) |>
  mutate(
    groups = case_when(
      `Variable / Model` == "p-value" ~ "Global",
      `Variable / Model` %in% models ~ "Econometric and ML models",
      TRUE ~ "Individual indicators"
    )
  ) |>
  mutate(
    `Variable / Model` = if_else(`Variable / Model` == "DFM", "IOAE", `Variable / Model`)
  ) |>
  relocate(groups) |>
  gt(groupname_col = "groups") |>
  tab_spanner(
    label = "h = 1",
    columns = 3:4
  ) |>
  tab_spanner(
    label = "h = 2",
    columns = 5:6
  ) |>
  cols_label(
    `TR h = 1` = html("T<sub>R,i</sub>"),
    `Rank h = 1` = html("(Rank<sub>R</sub>)"),
    `TR h = 2` = html("T<sub>R,i</sub>"),
    `Rank h = 2` = html("(Rank<sub>R</sub>)")
  ) |>
  sub_missing(
    columns = c(3, 5),
    rows = `Variable / Model` == "p-value",
    missing_text = ""
  ) |>
  cols_width(
    2 ~ px(160),
    everything() ~ px(120)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t_mcs_rmse.rds"))


#### table 4 ####
readRDS(here::here("outputs/rds/models_objects/df_dfm_err.rds")) |>
  arrange(window, r) |>
  slice(1:48) |>
  mutate(
    sample = rep(c(rep("Full sample", 4), rep("Excluding COVID-19", 4)), 6)
  ) |>
  relocate(sample) |>
  pivot_wider(
    #names_from = lbl,
    names_from = c(window, r, h),
    values_from = value
  ) |>
  gt(groupname_col = "sample") |>
  fmt_number(
    columns = 3:14,
    decimals = 3
  ) |>
  tab_spanner(
    label = "r = 1",
    columns = 3:4,
    id = "exp_r1"
  ) |>
  tab_spanner(
    label = "r = 2",
    columns = 5:6,
    id = "exp_r2"
  ) |>
  tab_spanner(
    label = "r = 3",
    columns = 7:8,
    id = "exp_r3"
  ) |>
  tab_spanner(
    label = "r = 1",
    columns = 9:10,
    id = "rol_r1"
  ) |>
  tab_spanner(
    label = "r = 2",
    columns = 11:12,
    id = "rol_r2"
  ) |>
  tab_spanner(
    label = "r = 3",
    columns = 13:14,
    id = "rol_r3"
  ) |>
  tab_spanner(
    label = "Expanding window",
    spanners = c("exp_r1", "exp_r2", "exp_r3")
  ) |>
  tab_spanner(
    label = "Rolling window",
    spanners =  c("rol_r1", "rol_r2", "rol_r3")
  ) |>
  cols_label(
    ends_with("1") ~ "h = 1",
    ends_with("2") ~ "h = 2",
    metric ~ ""
  ) |>
  cols_width(
    everything() ~ px(60)
  ) |>
  saveRDS(here::here("outputs/rds/prepared_tables/t4.rds"))


cat("\n\nReplication concluded!!\n")
cat("\n")
