#' Trace plots
#'
#' \code{trace_plots} implements function to plot Markov trace.
#'
#' @param outcomes List with model outcomes
#' @return 
#' main_states_trace_plot
#' sero_states_trace_plot
#' main_states_time
#' plots
#' layout
#' @export

trace_plots <- function(outcomes){
# Prepare data
df_M_agg_trace <- as.data.frame(outcomes$m_M_agg_trace)
df_M_agg_trace$month <- as.numeric(rownames(df_M_agg_trace))
df_M_agg_trace_plot <- df_M_agg_trace %>% gather(state, proportion, "Death", "ODF", "ODN", "REL", "BUP", "BUPC", "MET", "METC", "ABS") # health states to plot

df_M_agg_trace_sero <- as.data.frame(outcomes$m_M_agg_trace_sero)
df_M_agg_trace_sero$month <- as.numeric(rownames(df_M_agg_trace_sero))
df_M_agg_trace_sero_plot <- df_M_agg_trace_sero %>% gather(state, proportion, "NEG-Dead", "HIV-Dead", "HCV-Dead", "COI-Dead", "NEG-Alive", "HIV-Alive", "HCV-Alive", "COI-Alive") # health states to plot

df_M_agg_state_time <- df_M_agg_trace %>% gather(state, proportion, "ODF", "ODN", "REL", "BUP", "BUPC", "MET", "METC", "ABS") # alive health states to plot
df_M_agg_state_time <- df_M_agg_state_time %>% 
  group_by(state) %>% 
  summarise_each(funs(sum), proportion) %>%
  mutate(percentage = round((proportion / sum(proportion)) * 100,1))

# Preserve order for plotting
state_order_trace <- factor(df_M_agg_trace_plot$state, levels = c("Death", "ODF", "ODN", "REL", "BUP", "BUPC", "MET", "METC", "ABS"))
state_order_trace_sero <- factor(df_M_agg_trace_sero_plot$state, levels = c("NEG-Dead", "HIV-Dead", "HCV-Dead", "COI-Dead", "NEG-Alive", "HIV-Alive", "HCV-Alive", "COI-Alive"))
state_order_time  <- factor(df_M_agg_state_time$state, levels = c("ODF", "ODN", "REL", "BUP", "BUPC", "MET", "METC", "ABS"))
#state_colours_trace <- c("#d9d9d9", "#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598", "#abdda4", "#9ecae1", "#bcbddc") # colour pallette 1
state_colours_trace2 <- c("#d9d9d9", "#252525", "#b2182b", "#d6604d", "#d9f0d3", "#1b7837", "#d1e5f0", "#9ecae1", "#bcbddc") # colour pallette 2
state_colours_trace_sero <- c("#252525", "#cb181d", "#2171b5", "#6a51a3", "#969696", "#fc9272", "#9ecae1", "#bcbddc")

#state_colours_time <- c("#d53e4f", "#f46d43", "#fdae61", "#ffffbf", "#e6f598")
state_colours_time2 <- c("#252525", "#b2182b", "#fddbc7", "#d6604d", "#d9f0d3", "#1b7837", "#9ecae1", "#bcbddc") # colour pallette 2

### Markov trace plots ###
# Model 1 (Primary health state definition)
# Base model states
main_states_trace_plot <- ggplot(df_M_agg_trace_plot, aes(x = month, y = proportion, fill = state_order_trace)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Time (months)") + ylab("Proportion in state") +
  geom_area() +
  scale_fill_manual(name = "Health States", values = state_colours_trace2)

# Serostatus
sero_states_trace_plot <- ggplot(df_M_agg_trace_sero_plot, aes(x = month, y = proportion, fill = state_order_trace_sero)) + 
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Time (months)") + ylab("Proportion in state") +
  geom_area() +
  scale_fill_manual(name = "Serostatus", values = state_colours_trace_sero)

### Time spent in health states ###
main_states_time <- ggplot(df_M_agg_state_time, aes(x = state_order_time, y = proportion, fill = state_order_time)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Health State") + ylab("Time") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = state_colours_time2) +
  geom_text(aes(label = paste0(round(proportion,1)," (",percentage,"%)")), hjust = -0.25, size = 3.5) +
  coord_flip(ylim = c(0, 780))

### Combined plot ###
plots <- list()
plots[[1]] <- main_states_trace_plot
plots[[2]] <- main_states_time
layout <- matrix(c(1, 1, 2), nrow = 3, byrow = TRUE)

return(list(main_states_trace_plot,
            sero_states_trace_plot,
            main_states_time,
            plots,
            layout))
}