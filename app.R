# ==============================================================================
# R SHINY APP: Pulse & Pitch (Performer's Biofeedback Dashboard)
# ==============================================================================

# --- 1. Load All Necessary Packages ---
library(shiny)
library(shinyjs)
library(bslib)
library(ggplot2)
library(dplyr)
library(tuneR)
library(seewave)
library(av)
library(plotly)
library(htmlwidgets)
library(shinycssloaders)

# --- 2. Comprehensive Analysis Function ---
run_full_analysis <- function(hrv_path, audio_path, duration_s = 60) {
  tryCatch({
    rr_raw <- read.table(hrv_path, header = FALSE, col.names = "RR_ms")$RR_ms
    rr_corrected <- rr_raw
    percent_change <- abs(diff(rr_raw) / rr_raw[-length(rr_raw)])
    artifact_indices <- which(percent_change > 0.20) + 1
    if (length(artifact_indices) > 0) {
      for (idx in artifact_indices) {
        prev_val <- rr_corrected[idx - 1]; next_idx <- idx + 1
        while (next_idx %in% artifact_indices && next_idx <= length(rr_corrected)) next_idx <- next_idx + 1
        rr_corrected[idx] <- if (next_idx <= length(rr_corrected)) mean(c(prev_val, rr_corrected[next_idx])) else prev_val
      }
    }
    time_s <- cumsum(rr_corrected) / 1000; time_s <- time_s - time_s[1]; instant_hr <- 60000 / rr_corrected
    rolling_rmssd <- rep(NA_real_, length(rr_corrected)); win_max <- 5
    for (i in seq_along(rr_corrected)) {
      window_start_time <- max(0, time_s[i] - win_max)
      idx <- which(time_s <= time_s[i] & time_s > window_start_time)
      if (length(idx) >= 2) { d <- diff(rr_corrected[idx]); rolling_rmssd[i] <- sqrt(mean(d^2)) }
    }
    # Create the full, un-binned data first for accurate summaries
    hrv_data_full <- data.frame(Time_s = time_s, HeartRate_BPM = instant_hr, Rolling_RMSSD = rolling_rmssd) %>% 
      dplyr::filter(Time_s <= duration_s)
    
    wav_obj <- tuneR::readWave(audio_path); sr <- wav_obj@samp.rate; n_keep <- floor(duration_s * sr)
    if (length(wav_obj@left) > n_keep) wav_obj@left <- wav_obj@left[1:n_keep]
    win <- floor(1.0 * sr); n_win <- floor(length(wav_obj@left) / win)
    audio_features <- data.frame(Time_s = 0:(n_win - 1), Loudness_RMS = numeric(n_win), Brightness_Centroid = numeric(n_win), Mean_Pitch_Hz = numeric(n_win))
    for (k in seq_len(n_win)) {
      a <- (k - 1) * win + 1; b <- k * win
      chunk <- tuneR::Wave(left = wav_obj@left[a:b], samp.rate = sr, bit = wav_obj@bit)
      audio_features$Loudness_RMS[k] <- seewave::rms(chunk@left)
      pc <- seewave::fund(chunk, fmax = 2000, plot = FALSE)
      audio_features$Mean_Pitch_Hz[k] <- if (is.matrix(pc) && nrow(pc) > 0) mean(pc[,2], na.rm = TRUE) else NA_real_
      if (audio_features$Loudness_RMS[k] > 0.001) {
        sp <- seewave::spec(chunk, f = sr, wl = 2048, plot = FALSE); spp <- seewave::specprop(sp)
        audio_features$Brightness_Centroid[k] <- spp$cent
      } else audio_features$Brightness_Centroid[k] <- NA_real_
    }
    hrv_binned <- hrv_data_full %>%
      dplyr::mutate(Time_s = floor(Time_s)) %>%
      dplyr::group_by(Time_s) %>%
      dplyr::summarise(HR_mean = mean(HeartRate_BPM, na.rm = TRUE), RMSSD_mean = mean(Rolling_RMSSD, na.rm = TRUE), .groups = "drop")
    
    # Calculate summary stats from the un-binned data
    hr_summary <- list(
      min = round(min(hrv_data_full$HeartRate_BPM, na.rm = TRUE)),
      max = round(max(hrv_data_full$HeartRate_BPM, na.rm = TRUE)),
      avg = round(mean(hrv_data_full$HeartRate_BPM, na.rm = TRUE))
    )
    
    joint_data <- dplyr::full_join(audio_features, hrv_binned, by = "Time_s") %>%
      dplyr::arrange(Time_s) %>%
      dplyr::filter(Time_s >= 0, Time_s <= duration_s)
    
    # Return both the timeseries and the summary list
    return(list(timeseries = joint_data, summary = hr_summary))
    
  }, error = function(e) { print(e); return(NULL) })
}

compute_waveform_bars <- function(wav_path, duration_s, target_bins = 900L) {
  wav <- tuneR::readWave(wav_path); sr  <- wav@samp.rate; n   <- min(length(wav@left), floor(duration_s * sr))
  x   <- wav@left[1:n]; amp <- as.numeric(x) / (2^(wav@bit - 1)); bins <- target_bins
  idx  <- floor(seq(0, n, length.out = bins + 1)); idx[length(idx)] <- n
  peaks <- vapply(seq_len(bins), function(i) { if (idx[i+1] > idx[i]) max(abs(amp[(idx[i]+1):idx[i+1]])) else 0 }, numeric(1))
  t_mid <- ((idx[-1] + idx[-length(idx)]) / 2) / sr; bin_sec <- (n / sr) / bins; width_sec <- 0.9 * bin_sec
  rbind(data.frame(Time_s = t_mid, Amp =  peaks, width = width_sec), data.frame(Time_s = t_mid, Amp = -peaks, width = width_sec))
}

# --- 3. UI ---
ui <- fluidPage(
  useShinyjs(), 
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$script(src = "interactions.js"),
    tags$script(HTML("
      function enableScrolling() {
        document.body.classList.add('allow-scroll');
      }
    "))
  ),
  
  theme = bslib::bs_theme(bootswatch = "cyborg"),
  
  # Single centered container
  div(class = "main-container",
    # Title
    div(class = "title-container",
      h1("Pulse & Pitch"),
      div(class = "subtitle", "Performance Analysis Dashboard"),
      div(class = "project-description",
        p("This project was inspired by a deep-dive analysis of pianist Yuja Wang's marathon concert at Carnegie Hall, where I correlated her heart rate variability (HRV) with the performance. To make this type of analysis accessible to any performer or researcher, I developed this interactive R Shiny dashboard. The application allows users to upload their own HRV and audio data to visualize the real-time synchronization between physiological stress, heart rate, and key acoustic features like loudness and pitch."),
        br(),
        p("Learn more about the original Yuja Wang project and my other research at", 
          tags$a(href = "https://tristin.org", target = "_blank", "tristin.org", style = "color: #7DE0F6; text-decoration: none;"), ".")
      )
    ),
    
    # Upload section - centered
    div(class = "upload-section",
      h3("Upload & Analysis"),
      
      div(class = "file-input-container",
        fileInput("hrv_file", "1. Upload Performance HRV (.txt)")
      ),
      
      div(class = "file-input-container",
        fileInput("audio_file", "2. Upload Performance Audio (.m4a, etc.)")
      ),
      
      div(class = "numeric-input-container",
        numericInput("duration", "3. Analysis Duration (seconds):", value = 60, min = 10, max = 300)
      ),
      
      # Apply your custom class to the button
      actionButton("run_analysis", "4. Run Analysis", class = "btn-primary btn-lg btn-run-analysis")
    ),
    
    # Results section
    div(class = "results-section",
      uiOutput("main_content")
    )
  )
)

# --- 4. Server ---
server <- function(input, output, session) {
  output$main_content <- renderUI({
    
    # If the analysis_data() reactive is NULL (i.e., the app has just started or been reset),
    # show empty main panel
    if (is.null(analysis_data())) {
      
      div(
        style = "text-align: center; padding: 3rem 2rem; color: #e2e8f0; min-height: 400px; display: flex; align-items: center; justify-content: center;"
      )
      
    } else {
      
      # Otherwise, once analysis_data() has data, render the full results UI.
      tagList(
        uiOutput("audio_player_ui"),
        tabsetPanel(
          id = "results_tabs",
          
          tabPanel("Performance (Audio + HR)", 
                   shinycssloaders::withSpinner(
                     plotlyOutput("combo_wave_hr", height = "450px"), 
                     color = "#FF004D", type = 6
                   ),
                   fluidRow(
                     column(4, uiOutput("min_hr_box")),
                     column(4, uiOutput("max_hr_box")),
                     column(4, uiOutput("avg_hr_box"))
                   )
          ),
          
          tabPanel("Relative Stress Sync",
                   shinycssloaders::withSpinner(
                     plotlyOutput("relative_combo", height="680px"), 
                     color="#FFA500"
                   )
          ),
          
          tabPanel("Audio Features Sync",
                   shinycssloaders::withSpinner(plotlyOutput("loudness_plot", height = "250px"), color = "#32CD32"),
                   shinycssloaders::withSpinner(plotlyOutput("brightness_plot", height = "250px"), color = "#32CD32"),
                   shinycssloaders::withSpinner(plotlyOutput("pitch_plot", height = "250px"), color = "#32CD32")
          ),
          
          tabPanel("Joint Analysis & Correlation", 
                   # Example of adding the customizable lag slider here
                   sliderInput("corr_lag", "Correlation Lag (seconds):", 
                               min = -3, max = 3, value = -2, step = 1),
                   hr(), # Adds a horizontal line for separation
                   shinycssloaders::withSpinner(plotlyOutput("correlation_table", height="300px")),
                   shinycssloaders::withSpinner(plotlyOutput("joint_hrv_loudness_plot", height = "350px")),
                   shinycssloaders::withSpinner(plotlyOutput("joint_hrv_pitch_plot", height = "350px")),
                   shinycssloaders::withSpinner(plotlyOutput("joint_hrv_brightness_plot", height = "350px"))
          )
        )
      ) # End tagList
    } # End else
  }) # End renderUI
  
  output$audio_player_ui <- renderUI({
    req(audio_src())
    tags$audio(
      id = "audioplayer",
      src = audio_src(),
      controls = TRUE,
      preload = "auto",                    # <-- NEW
      style = "width: 100%; margin-bottom: 20px;"
    )
  })
  
  
  # Reactive values now store a list with timeseries and summary data
  analysis_data <- reactiveVal(NULL)
  audio_src <- reactiveVal(NULL)
  waveform_bars <- reactiveVal(NULL)
  
  observeEvent(input$run_analysis, {
    req(input$hrv_file, input$audio_file)
    showNotification("Processing files...", type = "message", duration = 5)
    temp_wav_path <- tempfile(fileext = ".wav")
    av::av_audio_convert(input$audio_file$datapath, output = temp_wav_path)
    # The analysis function now returns a list
    full_data <- run_full_analysis(input$hrv_file$datapath, temp_wav_path, duration_s = input$duration)
    if (is.null(full_data) || nrow(full_data$timeseries) == 0) { showNotification("Error: Analysis failed.", type = "error"); return() }
    analysis_data(full_data)
    wf_df <- compute_waveform_bars(temp_wav_path, duration_s = input$duration, target_bins = 900L)
    waveform_bars(wf_df %>% dplyr::filter(Time_s <= input$duration))
    if (!dir.exists("www")) dir.create("www")
    file.copy(temp_wav_path, "www/current_audio.wav", overwrite = TRUE)
    audio_src(paste0("current_audio.wav?", runif(1)))
    showNotification("Analysis complete.", type = "message")
    
    # Enable scrolling when analysis completes
    shinyjs::runjs("enableScrolling();")
  })
  

  dark_theme_visible_text <- reactive({
    theme_dark() + theme(
      plot.background  = element_rect(fill = "#1C1C27", color = NA),
      panel.background = element_rect(fill = "#1C1C27"),
      panel.grid.major = element_line(color = "#404040"),
      panel.grid.minor = element_blank(), text = element_text(color = "white"),
      plot.title = element_text(hjust = 0.5, size = 14, color = "white"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "white"),
      axis.title = element_text(color = "white", size = 12), axis.text = element_text(color = "white"),
      legend.background= element_rect(fill = "#1C1C27"), legend.key = element_rect(fill = "#1C1C27")
    )
  })
  
  # --- Tab 1: Performance (Audio + HR) - MODIFIED ---
  
  # Helper function to create the styled summary boxes
  create_summary_box <- function(value, label) {
    div(style = "text-align: center; color: white; padding: 10px;",
        div(style = "font-size: 44px; font-weight: bold; line-height: 1.1;", value),
        div(style = "font-size: 14px; text-transform: uppercase; color: #aaa;", label)
    )
  }
  
  # Render the summary boxes by accessing the 'summary' part of the results list
  output$min_hr_box <- renderUI({ req(analysis_data()); create_summary_box(analysis_data()$summary$min, "Min HR") })
  output$max_hr_box <- renderUI({ req(analysis_data()); create_summary_box(analysis_data()$summary$max, "Max HR") })
  output$avg_hr_box <- renderUI({ req(analysis_data()); create_summary_box(analysis_data()$summary$avg, "Average HR") })
  
  # The combo plot now accesses the 'timeseries' part of the results list
  output$combo_wave_hr <- renderPlotly({
    req(analysis_data(), waveform_bars())
    p_wave <- ggplot(waveform_bars(), aes(Time_s, Amp)) + geom_col(fill = "#7DE0F6", width = waveform_bars()$width[1], alpha = 0.95) + labs(x = NULL, y = NULL) + dark_theme_visible_text() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), plot.margin = margin(5, 10, 0, 5))
    p_hr <- ggplot(analysis_data()$timeseries, aes(x = Time_s, y = HR_mean)) + geom_line(color = "#FF004D", size = 1.0) + labs(x = "Time (seconds)", y = "BPM") + dark_theme_visible_text() + theme(plot.margin = margin(0, 10, 5, 5))
    gp_wave <- ggplotly(p_wave, tooltip = "x") %>%
      layout(
        yaxis = list(range = c(-1, 1), zeroline = TRUE, zerolinecolor = "#404040"),
        xaxis = list(range = c(0, input$duration))
      ) %>%
      config(displayModeBar = FALSE)
    
    gp_hr <- ggplotly(p_hr, tooltip = c("x","y")) %>%
      layout(xaxis = list(range = c(0, input$duration))) %>%
      config(displayModeBar = FALSE)
    
    combo <- subplot(gp_wave, gp_hr, nrows = 2, shareX = TRUE,
                     heights = c(0.38, 0.62), margin = 0.03) %>%
      layout(
        annotations = list(
          list(x = 0.5, y = 1.00, text = "Audio Waveform", showarrow = FALSE,
               xref = "paper", yref = "paper", yanchor = "bottom",
               font = list(size = 14, color = "white")),
          list(x = 0.5, y = 0.45, text = "Heart Rate", showarrow = FALSE,
               xref = "paper", yref = "paper", yanchor = "bottom",
               font = list(size = 14, color = "white"))
        ),
        clickmode = "event+select"  # reliable clicks on line traces
      )
    
    htmlwidgets::onRender(
      combo,
      "function(el, x){ if (window.syncPlotWithAudio) window.syncPlotWithAudio(el); }"
    )
  })
  
  relative_data <- reactive({ req(analysis_data()); base <- analysis_data()$timeseries %>% dplyr::select(Time_s, RMSSD_mean) %>% dplyr::arrange(Time_s); valid <- base %>% dplyr::filter(is.finite(RMSSD_mean)); if (nrow(valid) < 3) { return(list( data  = base %>% dplyr::mutate(Relative_Level = factor(NA_character_, levels = c("Relatively Relaxed","Near Session Average","Relatively High Strain"))), mean  = NA_real_, lower = NA_real_, upper = NA_real_, y_min = NA_real_, y_max = NA_real_ )) }; m <- mean(valid$RMSSD_mean, na.rm = TRUE); s <- sd(valid$RMSSD_mean, na.rm = TRUE); lower <- m - 0.5 * s; upper <- m + 0.5 * s; y_min <- min(valid$RMSSD_mean, na.rm = TRUE); y_max <- max(valid$RMSSD_mean, na.rm = TRUE); base <- base %>% dplyr::mutate( Relative_Level = dplyr::case_when( !is.finite(RMSSD_mean) ~ NA_character_, RMSSD_mean < lower ~ "Relatively High Strain", RMSSD_mean <= upper ~ "Near Session Average", TRUE ~ "Relatively Relaxed" ), Relative_Level = factor(Relative_Level, levels = c("Relatively Relaxed","Near Session Average","Relatively High Strain")) ); list(data = base, mean = m, lower = lower, upper = upper, y_min = y_min, y_max = y_max) })

  output$relative_combo <- renderPlotly({
    info <- relative_data()
    relative_colors <- c(
      "Relatively Relaxed"     = "#4575b4",
      "Near Session Average"   = "#abd9e9",
      "Relatively High Strain" = "#f46d43"
    )
    levs <- names(relative_colors)
    
    x_min <- 0
    x_max <- input$duration
    have_thr <- is.finite(info$mean) && is.finite(info$lower) &&
      is.finite(info$upper) && is.finite(info$y_min) && is.finite(info$y_max)
    
    # --------- build shading bands for geom_ribbon (needs full x sequence) ---------
    bands <- NULL
    if (have_thr) {
      xgrid <- seq(x_min, x_max, by = 1)
      n <- length(xgrid)
      bands <- dplyr::bind_rows(
        dplyr::tibble(Time_s = xgrid,
                      ymin   = rep(info$upper, n),
                      ymax   = rep(info$y_max, n),
                      Zone   = factor("Relatively Relaxed", levels = levs)),
        dplyr::tibble(Time_s = xgrid,
                      ymin   = rep(info$lower, n),
                      ymax   = rep(info$upper, n),
                      Zone   = factor("Near Session Average", levels = levs)),
        dplyr::tibble(Time_s = xgrid,
                      ymin   = rep(info$y_min, n),
                      ymax   = rep(info$lower, n),
                      Zone   = factor("Relatively High Strain", levels = levs))
      )
    }
    
    # --------- TOP: RMSSD with zones ---------
    # Midpoints of each zone
    y_relaxed <- (info$upper + info$y_max)/2
    y_avg     <- (info$lower + info$upper)/2
    y_strain  <- (info$y_min  + info$lower)/2
    
    p1 <- ggplot(info$data, aes(x = Time_s, y = RMSSD_mean)) +
      { if (!is.null(bands)) geom_ribbon(
        data = bands,
        aes(x = Time_s, ymin = ymin, ymax = ymax, fill = Zone),
        inherit.aes = FALSE, alpha = 0.5, color = NA
      ) } +
      geom_line(color = "white", linewidth = 1.2, na.rm = TRUE) +
      { if (have_thr) geom_hline(yintercept = info$mean,  color = "yellow", linetype = "dashed") } +
      { if (have_thr) geom_hline(yintercept = c(info$lower, info$upper),
                                 color = "yellow", linetype = "dotted") } +
      scale_fill_manual(name = "Relative Zone", values = relative_colors, drop = FALSE) +
      scale_y_continuous(
        name = "Relative Zone",
        breaks = c(y_strain, y_avg, y_relaxed),
        labels = c("High Strain RMSSD", "Average RMSSD", "Relaxed RMSSD")
      ) +
      dark_theme_visible_text() +
      theme(legend.position = "none")
    
    
    gp1 <- ggplotly(p1, tooltip = c("x","y")) %>%
      layout(
        margin = list(l = 70, r = 20, t = 60, b = 0),
        xaxis  = list(range = c(0, input$duration)),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.1)
      ) %>%
      config(displayModeBar = FALSE)
    
    # --------- BOTTOM: Relative state timeline (contiguous segments, no gaps) ---------
    df_runs <- info$data %>%
      dplyr::filter(Time_s >= 0, Time_s <= input$duration) %>%
      tidyr::complete(Time_s = 0:input$duration) %>%
      dplyr::arrange(Time_s) %>%
      tidyr::fill(Relative_Level, .direction = "down") %>%
      dplyr::filter(!is.na(Relative_Level)) %>%
      dplyr::mutate(Relative_Level = factor(Relative_Level, levels = levs)) %>%
      dplyr::mutate(run_id = cumsum(Relative_Level != dplyr::lag(Relative_Level,
                                                                 default = dplyr::first(Relative_Level)))) %>%
      dplyr::group_by(run_id) %>%
      dplyr::summarise(Relative_Level = dplyr::first(Relative_Level),
                       xstart = min(Time_s),
                       xend   = max(Time_s) + 1L,  # cover the last second fully
                       .groups = "drop")
    
    p2 <- ggplot(df_runs) +
      geom_segment(aes(x = xstart, xend = xend,
                       y = Relative_Level, yend = Relative_Level,
                       color = Relative_Level),
                   linewidth = 6, lineend = "butt") +
      scale_color_manual(values = relative_colors, name = "Relative Level") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
      labs(title = "RMSSD (Relative Stress) During Performance",
           subtitle = "Dashed = session mean; Dotted = ±0.5 SD",
           x = "Time (seconds)", y = "Relative State") +
      dark_theme_visible_text() +
      theme(legend.position = "none")
    
    gp2 <- ggplotly(p2, tooltip = c("x","colour")) %>%
      layout(
        margin = list(l = 70, r = 20, t = 50, b = 50),
        xaxis  = list(range = c(0, input$duration))
      ) %>%
      config(displayModeBar = FALSE)
    
    # --------- Combine (perfect alignment) ---------
    combo <- subplot(gp1, gp2, nrows = 2, shareX = TRUE,
                     heights = c(0.65, 0.35), margin = 0.05) %>%
      layout(plot_bgcolor = '#1C1C27', paper_bgcolor = '#1C1C27', clickmode = "event+select")
    
    htmlwidgets::onRender(combo, "function(el, x){ if (window.syncPlotWithAudio) window.syncPlotWithAudio(el); }")
  })
  
  
  output$loudness_plot <- renderPlotly({
    req(analysis_data())
    p <- ggplot(analysis_data()$timeseries, aes(x = Time_s, y = Loudness_RMS)) + geom_line(color = "#00CED1", size = 1.0) + labs(title = "Musical Loudness (Dynamics)", x = NULL, y = "RMS Energy") + dark_theme_visible_text()
    
    gp <- ggplotly(p) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(gp, "syncPlotWithAudio")
  })
  
  output$brightness_plot <- renderPlotly({
    req(analysis_data())
    p <- ggplot(na.omit(analysis_data()$timeseries), aes(x = Time_s, y = Brightness_Centroid)) + geom_line(color = "#FFA500", size = 1.0) + labs(title = "Sound Brightness (Timbre)", x = NULL, y = "Centroid (kHz)") + dark_theme_visible_text()
    
    gp <- ggplotly(p) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(gp, "syncPlotWithAudio")
  })
  
  output$pitch_plot <- renderPlotly({
    req(analysis_data())
    p <- ggplot(na.omit(analysis_data()$timeseries), aes(x = Time_s, y = Mean_Pitch_Hz)) + geom_line(color = "#32CD32", size = 1.0) + labs(title = "Pitch (Melody)", x = "Time (seconds)", y = "Frequency (Hz)") + dark_theme_visible_text()
    
    gp <- ggplotly(p) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(gp, "syncPlotWithAudio")

  })
  
  output$joint_hrv_loudness_plot <- renderPlotly({
    req(analysis_data())
    data <- analysis_data()$timeseries
    
    p <- plot_ly(data, x = ~Time_s) %>% 
      add_trace(y = ~HR_mean, name = 'Heart Rate', type = 'scatter', mode = 'lines', line = list(color = '#FF004D', width = 3)) %>% 
      add_trace(y = ~Loudness_RMS, name = 'Loudness (RMS)', type = 'scatter', mode = 'lines', line = list(color = '#00CED1', width = 3), yaxis = 'y2') %>% 
      layout(title = list(text = "Joint Analysis: Heart Rate vs. Loudness", font = list(color="white")), yaxis = list(title = list(text="Heart Rate", font=list(color="#00BFFF")), color = "white", gridcolor = "#404040"), yaxis2 = list(title = list(text="Loudness (RMS)", font=list(color="#FF004D")), overlaying = 'y', side = 'right', color = "white", showgrid=FALSE), xaxis = list(title = "Time (seconds)", color = "white", gridcolor = "#404040"), plot_bgcolor = '#1C1C27', paper_bgcolor = '#1C1C27', legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2, font = list(color='white')), margin = list(r = 80)) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(p, "syncPlotWithAudio")
  })
  
  output$joint_hrv_pitch_plot <- renderPlotly({
    req(analysis_data())
    data <- analysis_data()$timeseries
    
    p <- plot_ly(data, x = ~Time_s) %>% 
      add_trace(y = ~RMSSD_mean, name = 'RMSSD (ms)', type = 'scatter', mode = 'lines', line = list(color = '#00BFFF', width = 3)) %>% 
      add_trace(y = ~Mean_Pitch_Hz, name = 'Pitch (Hz)', type = 'scatter', mode = 'lines', line = list(color = '#32CD32', width = 3), yaxis = 'y2') %>% 
      layout(title = list(text = "Joint Analysis: RMSSD vs. Pitch", font = list(color="white")), yaxis = list(title = list(text="RMSSD (ms)", font=list(color="#00BFFF")), color = "white", gridcolor = "#404040"), yaxis2 = list(title = list(text="Pitch (Hz)", font=list(color="#32CD32")), overlaying = 'y', side = 'right', color = "white", showgrid=FALSE), xaxis = list(title = "Time (seconds)", color = "white", gridcolor = "#404040"), plot_bgcolor = '#1C1C27', paper_bgcolor = '#1C1C27', legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2, font = list(color='white')), margin = list(r = 80)) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(p, "syncPlotWithAudio")
  })
  
  output$joint_hrv_brightness_plot <- renderPlotly({
    req(analysis_data())
    data <- analysis_data()$timeseries
    
    p <- plot_ly(data, x = ~Time_s) %>% 
      add_trace(y = ~RMSSD_mean, name = 'RMSSD (ms)', type = 'scatter', mode = 'lines', line = list(color = '#00BFFF', width = 3)) %>% 
      add_trace(y = ~Brightness_Centroid, name = 'Brightness (kHz)', type = 'scatter', mode = 'lines', line = list(color = '#FFA500', width = 3), yaxis = 'y2') %>% 
      layout(title = list(text = "Joint Analysis: RMSSD vs. Brightness", font = list(color="white")), yaxis = list(title = list(text="RMSSD (ms)", font=list(color="#00BFFF")), color = "white", gridcolor = "#404040"), yaxis2 = list(title = list(text="Brightness (kHz)", font=list(color="#FFA500")), overlaying = 'y', side = 'right', color = "white", showgrid=FALSE), xaxis = list(title = "Time (seconds)", color = "white", gridcolor = "#404040"), plot_bgcolor = '#1C1C27', paper_bgcolor = '#1C1C27', legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2, font = list(color='white')), margin = list(r = 80)) %>% 
      config(displayModeBar = FALSE)
    
    htmlwidgets::onRender(p, "syncPlotWithAudio")
  })
  
  output$correlation_table <- renderPlotly({ req(analysis_data()); dat <- na.omit(analysis_data()$timeseries); if (nrow(dat) < 3) return(NULL); 
  current_lag <- input$corr_lag; dynamic_title <- sprintf("Cross-Correlation at %ds Lag (Physiology vs Audio)", current_lag);
  safe_lag_cor <- function(audio_vec, physio_vec, lag_seconds = -2) { if(sd(audio_vec, na.rm=TRUE) > 0 && sd(physio_vec, na.rm=TRUE) > 0) { df <- data.frame(audio = audio_vec, physio = physio_vec); df$physio_lagged <- dplyr::lead(df$physio, n = abs(lag_seconds)); cor(df$audio, df$physio_lagged, use = "pairwise.complete.obs") } else { NA } }; lag_val <- input$corr_lag; r_hrv_l <- safe_lag_cor(dat$Loudness_RMS, dat$RMSSD_mean, lag_val); r_hr_l  <- safe_lag_cor(dat$Loudness_RMS, dat$HR_mean,    lag_val); r_hrv_p <- safe_lag_cor(dat$Mean_Pitch_Hz, dat$RMSSD_mean, lag_val); r_hr_p  <- safe_lag_cor(dat$Mean_Pitch_Hz, dat$HR_mean,    lag_val); r_hrv_b <- safe_lag_cor(dat$Brightness_Centroid, dat$RMSSD_mean, lag_val); r_hr_b  <- safe_lag_cor(dat$Brightness_Centroid, dat$HR_mean,    lag_val); fmt <- function(x) ifelse(is.na(x), "N/A", sprintf("%.2f", x)); stronger_l <- if (is.na(r_hrv_l) || is.na(r_hr_l)) "N/A" else if (abs(r_hrv_l) >= abs(r_hr_l)) "<b>HRV</b>" else "<b>HR</b>"; stronger_p <- if (is.na(r_hrv_p) || is.na(r_hr_p)) "N/A" else if (abs(r_hrv_p) >= abs(r_hr_p)) "<b>HRV</b>" else "<b>HR</b>"; stronger_b <- if (is.na(r_hrv_b) || is.na(r_hr_b)) "N/A" else if (abs(r_hrv_b) >= abs(r_hr_b)) "<b>HRV</b>" else "<b>HR</b>"; plot_ly( type = "table", header = list( values = c("<b>Metric</b>", "<b>Heart Rate Variability (HRV, RMSSD)</b>", "<b>Heart Rate</b>", "<b>Stronger Correlation</b>"), fill   = list(color = "#2a3f5f"), font   = list(color = "white", size = 14), align  = c("left", rep("center", 3)), line = list(color = "#506784", width = 1) ), cells = list( values = list( c("Loudness", "Pitch", "Brightness"), c(fmt(r_hrv_l), fmt(r_hrv_p), fmt(r_hrv_b)), c(fmt(r_hr_l), fmt(r_hr_p), fmt(r_hr_b)), c(stronger_l, stronger_p, stronger_b) ), fill  = list(color = "#1C1C27"), font  = list(color = "white", size = 13), align = c("left", rep("center", 3)), line = list(color = "#506784", width = 1), height = 30 ) ) %>% layout( title = list(text = dynamic_title, font = list(color="white"), y=0.95), paper_bgcolor = "#1C1C27", plot_bgcolor  = "#1C1C27", margin = list(l = 20, r = 20, t = 50, b = 10) ) %>% config(displayModeBar = FALSE) })

}

# --- 5. Run the Application ---
shinyApp(ui = ui, server = server)