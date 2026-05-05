source("riley_continuous_sample_size.R")

library(shiny)
library(bslib)

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

theme <- bs_theme(
  version = 5,
  bootswatch = "morph",
  primary = "#3563E9",
  secondary = "#E56B6F",
  base_font = "system-ui",
  heading_font = "system-ui"
)

ui <- page_sidebar(
  theme = theme,
  title = "Riley Continuous Sample Size",
  sidebar = sidebar(
    width = 340,
    class = "sidebar-glass",
    h4("Inputs"),
    sliderInput("p", "Predictor parameters (p)", min = 2, max = 100, value = 25, step = 1),
    sliderInput("r2_adj", "Anticipated adjusted R²", min = 0.01, max = 0.95, value = 0.20, step = 0.01),
    sliderInput("target_sc", "Target shrinkage", min = 0.80, max = 0.99, value = 0.90, step = 0.01),
    sliderInput("delta_r2", "Max. R² optimism", min = 0.01, max = 0.20, value = 0.05, step = 0.01),
    sliderInput("target_mmoe", "Target MMOE for sigma", min = 1.01, max = 1.30, value = 1.10, step = 0.01),
    sliderInput("alpha", "Confidence level alpha", min = 0.01, max = 0.20, value = 0.05, step = 0.01),
    hr(),
    actionButton("preset", "Load paper example", class = "btn-primary w-100")
  ),
  layout_columns(
    col_widths = c(4, 4, 4),
    card(
      card_header("Criterion I"),
      card_body(
        div(class = "value-box", textOutput("n_i")),
        p("Sample size for expected shrinkage.")
      )
    ),
    card(
      card_header("Criterion II"),
      card_body(
        div(class = "value-box", textOutput("n_ii")),
        p("Sample size for R² optimism threshold.")
      )
    ),
    card(
      card_header("Criterion III"),
      card_body(
        div(class = "value-box", textOutput("n_iii")),
        p("Sample size for sigma precision.")
      )
    )
  ),
  layout_columns(
    col_widths = c(5, 7),
    card(
      card_header("Recommended minimum"),
      card_body(
        div(class = "hero-n", textOutput("n_min")),
        p("This is the largest sample size across the criteria.")
      )
    ),
    card(
      card_header("Quick plot"),
      card_body(plotOutput("curve", height = 260))
    )
  ),
  card(
    card_header("Optional intercept check"),
    card_body(
      layout_columns(
        col_widths = c(4, 4, 4),
        numericInput("alpha_hat", "Anticipated mean outcome", value = 1.90, min = -1e6),
        numericInput("sigma2_null", "Null-model variance", value = 0.36, min = 1e-12),
        uiOutput("ci_box")
      )
    )
  ),
  card(
    card_header("Notes"),
    card_body(
      tags$ul(
        tags$li("Criterion I solves the paper's shrinkage equation numerically."),
        tags$li("Criterion II and III use the paper's closed-form expressions."),
        tags$li("For publication, this app can be deployed as a standard Shiny app.")
      )
    )
  ),
  tags$style(HTML("
    body {
      background:
        radial-gradient(circle at top left, rgba(53, 99, 233, 0.16), transparent 35%),
        radial-gradient(circle at right, rgba(229, 107, 111, 0.14), transparent 32%),
        linear-gradient(180deg, #f8fbff 0%, #f5f7ff 100%);
    }
    .sidebar-glass {
      background: rgba(255, 255, 255, 0.72);
      backdrop-filter: blur(12px);
      border-right: 1px solid rgba(20, 30, 60, 0.08);
    }
    .value-box {
      font-size: 2.2rem;
      font-weight: 800;
      line-height: 1;
      color: #1f2b4d;
      letter-spacing: -0.03em;
    }
    .hero-n {
      font-size: 3.2rem;
      font-weight: 900;
      line-height: 1;
      color: #3563E9;
      letter-spacing: -0.05em;
    }
  "))
)

server <- function(input, output, session) {
  observeEvent(input$preset, {
    updateSliderInput(session, "p", value = 25)
    updateSliderInput(session, "r2_adj", value = 0.20)
    updateSliderInput(session, "target_sc", value = 0.90)
    updateSliderInput(session, "delta_r2", value = 0.05)
    updateSliderInput(session, "target_mmoe", value = 1.10)
    updateSliderInput(session, "alpha", value = 0.05)
    updateNumericInput(session, "alpha_hat", value = 1.90)
    updateNumericInput(session, "sigma2_null", value = 0.36)
  })

  results <- reactive({
    list(
      i = riley_continuous_n_shrinkage(input$p, input$r2_adj, input$target_sc),
      ii = riley_continuous_n_r2diff(input$p, input$r2_adj, input$delta_r2),
      iii = riley_continuous_n_sigma(input$p, input$target_mmoe, input$alpha),
      min = riley_continuous_min_sample_size(
        input$p, input$r2_adj, input$target_sc,
        input$delta_r2, input$target_mmoe, input$alpha
      )
    )
  })

  output$n_i <- renderText(fmt_n(results()$i))
  output$n_ii <- renderText(fmt_n(results()$ii))
  output$n_iii <- renderText(fmt_n(results()$iii))
  output$n_min <- renderText(fmt_n(results()$min$minimum))

  output$curve <- renderPlot({
    n_grid <- seq(max(input$p + 3, 10), max(results()$min$minimum * 1.5, input$p + 20))
    sc <- vapply(n_grid, function(n) riley_continuous_shrinkage(n, input$p, input$r2_adj), numeric(1))
    plot(
      n_grid, sc, type = "l", lwd = 3, col = "#3563E9",
      xlab = "Sample size (n)", ylab = "Expected shrinkage",
      ylim = c(min(sc, input$target_sc) - 0.02, 1.01),
      xaxs = "i", yaxs = "i", bty = "l"
    )
    abline(h = input$target_sc, col = "#E56B6F", lty = 2, lwd = 2)
    abline(v = results()$min$minimum, col = "#1f2b4d", lty = 3, lwd = 2)
    legend(
      "bottomright",
      legend = c("Expected shrinkage", "Target", "Recommended n"),
      col = c("#3563E9", "#E56B6F", "#1f2b4d"),
      lty = c(1, 2, 3), lwd = c(3, 2, 2), bty = "n", cex = 0.9
    )
  }, res = 96)

  output$ci_box <- renderUI({
    ci <- riley_continuous_intercept_ci(
      n = results()$min$minimum,
      p = input$p,
      alpha_hat = input$alpha_hat,
      sigma2_null = input$sigma2_null,
      r2_adj = input$r2_adj,
      alpha = input$alpha
    )
    card(
      card_header("95% CI for intercept"),
      card_body(
        div(class = "value-box", sprintf("%.3f to %.3f", ci[["lower"]], ci[["upper"]])),
        p("Using the current recommended minimum sample size.")
      )
    )
  })
}

shinyApp(ui, server)
