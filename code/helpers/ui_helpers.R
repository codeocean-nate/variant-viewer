# Shared UI helpers for empty-state rendering

EMPTY_STATE_MESSAGE <- "No variant data loaded. Use the Data panel in the sidebar to select a VCF file from /data/, then click Process VCF."

empty_state_plotly <- function(msg = EMPTY_STATE_MESSAGE) {
  plotly::plot_ly() %>%
    plotly::layout(
      title = list(text = msg, x = 0.5, xanchor = "center", y = 0.5, yanchor = "middle"),
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE)
    )
}

empty_state_baseplot <- function(msg = EMPTY_STATE_MESSAGE) {
  plot.new()
  text(0.5, 0.5, msg, cex = 1.1, col = "gray40")
}

empty_state_html <- function(msg = EMPTY_STATE_MESSAGE) {
  htmltools::p(msg, style = "color: gray; font-size: 14px; padding: 20px;")
}
